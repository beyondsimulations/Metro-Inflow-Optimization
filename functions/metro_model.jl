"""
    build_optimization_model(modelInstance)

Purpose: Creates a mathematical optimization model to determine optimal passenger inflow rates across the metro network.
Details: Constructs a JuMP optimization model that:
1. Sets variables X[o,p] representing inflow rates at each origin station o for each period p
2. Creates an objective function minimizing the difference between passenger demand and allowed inflow
3. Adds capacity constraints to ensure passenger flows don't exceed arc capacities
4. Configures the solver (HiGHS) with time limits and optimality gap parameters

Returns: Two items: the optimization model object and the decision variables X
"""
function build_optimization_model(modelInstance)

    # Initialize optimization model object
    im = Model(HiGHS.Optimizer)

    # Set attributes for presolving, time limit, and mip gap
    set_attribute(im, "presolve", "on")
    set_attribute(im, "time_limit", 120.0)
    set_attribute(im, "mip_rel_gap", 0.0)

    # Precompute total demand per origin per period (avoids repeated sum() calls)
    total_demand_od = [sum(modelInstance.demand_od_in_period[o, :, p]) for o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods]

    println("Preparing optimization model.")
    @variable(im,
        modelInstance.min_entry_origin .<= X[o=1:modelInstance.nr_nodes, p=1:modelInstance.nr_periods] .<= modelInstance.max_entry_origin * modelInstance.safety_factor
    )

    println("Preparing objective function.")
    @objective(im, Min,
        sum((sum(modelInstance.cum_demand_od_in_period[o, d, p] for d in 1:modelInstance.nr_nodes) - (X[o, p] * modelInstance.minutes_in_period)) for o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods)
    )

    println("Preparing capacity constraints.")
    @constraint(im, capacity[a in 1:modelInstance.nr_arcs, t in 1:modelInstance.nr_minutes, p_shifts in 0:ceil(Int, modelInstance.past_minutes / modelInstance.minutes_in_period); shift[a, t] != []],
        sum(X[o, p] * modelInstance.demand_od_in_period[o, d, max(1, p - p_shifts)] / total_demand_od[o, max(1, p - p_shifts)] for (o, d, p) in modelInstance.shift[a, t] if modelInstance.demand_od_in_period[o, d, max(1, p - p_shifts)] > 0 && total_demand_od[o, max(1, p - p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
    )

    return im, X
end

"""
    build_restricted_optimization_model(modelInstance, current_period, queue_period_age, inflow_raw)

Purpose: Creates a time-restricted optimization model focusing on a specific period and its temporal neighborhood.
Details: Similar to build_optimization_model() but with several key differences:
1. Considers only a slice of the time horizon centered around the current period
2. Handles cumulative demand specifically for the restricted time window
3. Supports different objective function types ("regularSqr" for squared differences or "linweight" for queue-age weighted linear difference)
4. Implements two different approaches to manage queues:
   - "shift_periods": Uses time-shifted demand based on fixed periods
   - "lag_periods": Adjusts demand based on how long passengers have been waiting

Returns: Two items: the restricted optimization model object and the decision variables X
"""
function build_restricted_optimization_model(modelInstance, current_period, queue_period_age, inflow_raw)
    # Initialize optimization model object
    im = Model(HiGHS.Optimizer)

    # Set attributes for presolving, time limit, and mip gap
    set_attribute(im, "presolve", "on")
    set_attribute(im, "time_limit", 120.0)
    set_attribute(im, "mip_rel_gap", 0.0)

    # Prepare the upper and lower bound of the observed time slice
    add_period = ceil(longest_path / modelInstance.minutes_in_period)
    lower_period::Int64 = max(1, current_period - ceil((longest_path + modelInstance.minutes_in_period) / modelInstance.minutes_in_period))
    upper_period::Int64 = min(current_period + add_period, modelInstance.nr_periods)
    lower_minute::Int64 = max(1, ceil(modelInstance.minutes_in_period * (current_period - 1) + 1))
    upper_minute::Int64 = min(ceil(modelInstance.minutes_in_period * (current_period + add_period)) + longest_path, modelInstance.nr_minutes)

    adjusted_cum_demand = copy(modelInstance.cum_demand_od_in_period)
    if current_period < upper_period
        for p in current_period+1:upper_period
            adjusted_cum_demand[:, :, p] .+= adjusted_cum_demand[:, :, p-1]
        end
    end

    # Precompute total demand per origin per period (avoids repeated sum() calls)
    println("Precomputing demand totals and ratios...")
    total_demand = zeros(Float64, modelInstance.nr_nodes, modelInstance.nr_periods)
    total_demand_od = zeros(Float64, modelInstance.nr_nodes, modelInstance.nr_periods)
    demand_ratio = zeros(Float64, modelInstance.nr_nodes, modelInstance.nr_nodes, modelInstance.nr_periods)
    @showprogress desc="Computing demand totals..." for o in 1:modelInstance.nr_nodes
        for p in 1:modelInstance.nr_periods
            total_demand[o, p] = sum(adjusted_cum_demand[o, :, p])
            total_demand_od[o, p] = sum(modelInstance.demand_od_in_period[o, :, p])
            # Precompute demand ratios: demand[o,d,p] / total_demand_od[o,p]
            if total_demand_od[o, p] > 0
                for d in 1:modelInstance.nr_nodes
                    demand_ratio[o, d, p] = modelInstance.demand_od_in_period[o, d, p] / total_demand_od[o, p]
                end
            end
        end
    end

    println("Preparing optimization model.")
    @variable(im,
        min(total_demand[o, p] / modelInstance.minutes_in_period, modelInstance.min_entry_origin) .<= X[o=1:modelInstance.nr_nodes, p=lower_period:upper_period] .<= modelInstance.max_entry_origin * modelInstance.safety_factor
    )

    if modelInstance.kind_opt == "regularSqr"
        println("Preparing regular objective function.")
        @objective(im, Min,
            sum((total_demand[o, p] - (X[o, p] * modelInstance.minutes_in_period))^2 for o in 1:modelInstance.nr_nodes, p in current_period:upper_period)
        )

    elseif modelInstance.kind_opt == "linweight"
        println("Preparing linear objective function.")
        @objective(im, Min,
            sum((total_demand[o, p] - X[o, p] * modelInstance.minutes_in_period) * queue_period_age[o, p] for o in 1:modelInstance.nr_nodes, p in current_period:upper_period)
        )

        println("Preparing constraint to prevent negative dispatch.")
        @constraint(
            im, queue[o in 1:modelInstance.nr_nodes, p in current_period:upper_period],
            total_demand[o, p] - X[o, p] * modelInstance.minutes_in_period >= 0
        )

    end

    if modelInstance.kind_queue == "shift_periods"
        # Fix the capacity inflow based on the fixed values of the past and the future (if min_entry_origin > 0)
        periods_to_consider = ceil(Int, modelInstance.past_minutes / modelInstance.minutes_in_period)
        adjusted_capacity_arcs = zeros(Float64, modelInstance.nr_minutes, modelInstance.nr_arcs)

        minute_range = lower_minute:upper_minute

        for t in minute_range
            adjusted_capacity_arcs[t, :] = modelInstance.capacity_arcs * modelInstance.safety_factor
        end

        # Precompute active (arc, minute) pairs to skip empty cells
        active_shifts = [(a, t) for a in 1:modelInstance.nr_arcs for t in minute_range if !isempty(shift[a, t])]
        println("Active (arc, minute) pairs: $(length(active_shifts)) / $(modelInstance.nr_arcs * length(minute_range))")

        println("Computing adjusted capacities (threads: $(Threads.nthreads()))...")
        progress = Progress(length(active_shifts), desc="Adjusting arc capacities...")
        min_entry = modelInstance.min_entry_origin
        Threads.@threads for idx in eachindex(active_shifts)
            (a, t) = active_shifts[idx]
            local_max = adjusted_capacity_arcs[t, a]
            for p_shifts in 0:periods_to_consider
                arcweight = 0.0
                @inbounds for (o, d, p) in shift[a, t]
                    p_idx = max(1, p - p_shifts)
                    ratio = demand_ratio[o, d, p_idx]
                    if ratio > 0
                        base = p < current_period ? inflow_raw[o, p] : min_entry
                        arcweight += base * ratio
                    end
                end
                local_max = max(local_max, arcweight)
            end
            adjusted_capacity_arcs[t, a] = local_max
            next!(progress)
        end
        println("Preparing capacity constraints for periodical demand.")
        @constraint(im, capacity_period[a in 1:modelInstance.nr_arcs, t in lower_minute:upper_minute, p_shifts in 0:ceil(Int, modelInstance.past_minutes / minutes_in_period); shift[a, t] != []],
            sum(X[o, p] * demand_ratio[o, d, max(1, p - p_shifts)] for (o, d, p) in modelInstance.shift[a, t] if demand_ratio[o, d, max(1, p - p_shifts)] > 0 && p <= upper_period) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
        )

    elseif modelInstance.kind_queue == "lag_periods"
        # Fix the capacity inflow based on the fixed values of the past and the future (if min_entry_origin > 0)
        periods_to_consider = ceil(Int, modelInstance.past_minutes / modelInstance.minutes_in_period)
        adjusted_capacity_arcs = zeros(Float64, modelInstance.nr_minutes, modelInstance.nr_arcs)

        minute_range = lower_minute:upper_minute

        for t in minute_range
            adjusted_capacity_arcs[t, :] = modelInstance.capacity_arcs * modelInstance.safety_factor
        end

        # Precompute active (arc, minute) pairs to skip empty cells
        active_shifts = [(a, t) for a in 1:modelInstance.nr_arcs for t in minute_range if !isempty(shift[a, t])]
        println("Active (arc, minute) pairs: $(length(active_shifts)) / $(modelInstance.nr_arcs * length(minute_range))")

        println("Computing adjusted capacities (threads: $(Threads.nthreads()))...")
        progress = Progress(length(active_shifts), desc="Adjusting arc capacities...")
        min_entry = modelInstance.min_entry_origin
        nr_periods = modelInstance.nr_periods
        Threads.@threads for idx in eachindex(active_shifts)
            (a, t) = active_shifts[idx]
            local_max = adjusted_capacity_arcs[t, a]
            for p_shifts in 0:periods_to_consider
                arcweight = 0.0
                @inbounds for (o, d, p) in shift[a, t]
                    p_idx = min(p, max(1, p - queue_period_age[o, p] + p_shifts), nr_periods)
                    ratio = demand_ratio[o, d, p_idx]
                    if ratio > 0
                        base = p < current_period ? inflow_raw[o, p] : min_entry
                        arcweight += base * ratio
                    end
                end
                local_max = max(local_max, arcweight)
            end
            adjusted_capacity_arcs[t, a] = local_max
            next!(progress)
        end

        println("Preparing capacity constraints for periodical demand shifted to the queue length.")
        @constraint(im, capacity_period[a in 1:modelInstance.nr_arcs, t in minute_range, p_shifts in 0:periods_to_consider; shift[a, t] != []],
            sum(X[o, p] * demand_ratio[o, d, min(p, max(1, p - queue_period_age[o, p] + p_shifts), modelInstance.nr_periods)] for (o, d, p) in modelInstance.shift[a, t] if demand_ratio[o, d, min(p, max(1, p - queue_period_age[o, p] + p_shifts), modelInstance.nr_periods)] > 0 && p <= upper_period) <= adjusted_capacity_arcs[t, a]
        )
    end

    return im, X
end

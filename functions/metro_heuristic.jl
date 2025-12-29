"""
    heuristic_adding_queues(im, config)

Purpose: Implements a time-based heuristic to optimize passenger inflow across metro stations while managing queues.
Details: Works period by period to:
1. Track remaining passenger queues at each origin station
2. Build and solve an optimization model for each time period
3. Determine how many passengers can enter the system from each station
4. Allocate available capacity to waiting passengers based on their destinations
5. Track queue lengths and waiting times
6. Calculate network utilization metrics

The function maintains a running queue of passengers who couldn't enter in previous periods and prioritizes them in subsequent periods.

Returns: Five items:
- results_queues: DataFrame of station-level statistics (allowed entries, moved passengers, queue length)
- results_arcs: DataFrame of connection-level utilization metrics
- optimization_duration: Time taken to solve each period's optimization model
- queue_period_age: Matrix tracking how long passengers have been waiting at each station
- infeasible_solutions: Count of periods where no feasible solution was found
"""
function heuristic_adding_queues(im, config)
    # Initialization of variables
    remaining_queue = copy(im.demand_od_in_period) # Initialize a new variable to store the remaining queue data
    inflow_raw = zeros(Float64, im.nr_nodes, im.nr_periods) .= im.min_entry_origin  # Initialize an array to store the values of the X variable
    optimization_duration = zeros(Float64, im.nr_periods) # Initialize optimization duration vector
    build_duration = zeros(Float64, im.nr_periods) # Initialize model build duration vector
    queue_period_age = zeros(Int64, im.nr_nodes, im.nr_periods) .= 1 # Computes the age of each queue per period (min 1 for optimization weighting)
    save_queue = zeros(Int64, nr_nodes, im.nr_periods) .= 1 # Saves the queue for the output

    stats = DataFrame(
        period=Int64[],
        new_demand=Float64[],
        in_queue=Float64[],
        moved_demand=Float64[],
    )

    infeasible_solutions = 0

    @showprogress for fix_period in 1:im.nr_periods

        println()
        println("Running period ", fix_period)
        im.cum_demand_od_in_period[:, :, fix_period] .= sum(remaining_queue[:, :, 1:fix_period], dims=3)

        # Prepare the upper and lower bound of the observed time slice
        lower_period::Int64 = max(1, fix_period - ceil((longest_path + im.minutes_in_period) / im.minutes_in_period))
        upper_period::Int64 = min(fix_period, im.nr_periods)

        println(sum(sum(remaining_queue[:, :, 1:fix_period], dims=3)))

        for o in 1:nr_nodes
            for queue_period in 1:fix_period
                if sum(remaining_queue[o, :, queue_period]) >= 1
                    queue_period_age[o, fix_period] = max(1, fix_period - queue_period)
                    break
                end
            end
        end

        println("Queue Age")
        println(queue_period_age[:, fix_period])
        println("Queue Length")
        println(round.(Int64, sum(im.cum_demand_od_in_period[:, :, fix_period], dims=2)))
        println("Total Queue Length")
        println(round.(Int64, sum(im.cum_demand_od_in_period[:, :, fix_period])))

        if im.closed_period[fix_period] == false
            # Only build and solve model during open periods
            build_duration[fix_period] = @elapsed begin
                model, X = build_restricted_optimization_model(im, fix_period, queue_period_age, inflow_raw)
            end

            for o in 1:im.nr_nodes
                if fix_period > 1
                    for p in lower_period:fix_period-1
                        if inflow_raw[o, p] > 0.0001
                            fix(X[o, p], inflow_raw[o, p]; force=true)
                        else
                            fix(X[o, p], 0; force=true)
                        end
                    end
                end
            end

            optimization_duration[fix_period] = @elapsed optimize!(model)

            if is_solved_and_feasible(model) == true
                for o in 1:im.nr_nodes
                    if value.(X)[o, fix_period] > 0.0001
                        inflow_raw[o, fix_period] = value.(X)[o, fix_period]
                    else
                        inflow_raw[o, fix_period] = 0.00
                    end
                end
            else
                infeasible_solutions += 1
                inflow_raw[:, fix_period] .= 0
            end

            # Release model memory to prevent RAM accumulation
            empty!(model)
            GC.gc()
        else
            # Closed period: skip model building entirely
            println("Skipping closed period $fix_period")
            inflow_raw[:, fix_period] .= 0
        end

        demand_fulfiled = inflow_raw[:, fix_period] .* im.minutes_in_period
        println("Dispatch")
        println(round.(demand_fulfiled))
        println("Total Dispatch")
        println(sum(demand_fulfiled))

        for o in 1:im.nr_nodes
            for p in 1:fix_period
                current_ratio = (remaining_queue[o, :, p] / sum(remaining_queue[o, :, p]))
                current_ratio .= ifelse.(isnan.(current_ratio), 0, current_ratio)

                while sum(remaining_queue[o, :, p]) >= 1 && demand_fulfiled[o] >= 1
                    remaining_queue[o, :, p] .-= im.minutes_in_period .* current_ratio
                    demand_fulfiled[o] -= im.minutes_in_period
                end

                if sum(remaining_queue[o, :, p]) <= 1
                    remaining_queue[o, :, p] .= 0
                end

                if demand_fulfiled[o] <= 1
                    demand_fulfiled[o] = 0
                end

                if demand_fulfiled[o] == 0
                    break
                end
            end
        end

        if im.closed_period[fix_period] == true
            remaining_queue[:, :, 1:fix_period] .= 0
        end
        for o in 1:im.nr_nodes
            save_queue[o, fix_period] = round(Int64, sum(remaining_queue[o, :, 1:fix_period]))
        end

        push!(stats, (
            period=fix_period,
            new_demand=sum(im.demand_od_in_period[:, :, fix_period]),
            in_queue=sum(save_queue[:, fix_period]),
            moved_demand=sum(inflow_raw[:, fix_period] .* minutes_in_period),
        )
        )
    end

    # Save optimization diagnostic plots
    plot_dir = "results/plots/$(config.name)_$(im.safety_factor)_$(im.minutes_in_period)_$(im.past_minutes)_$(im.max_entry_origin)_$(im.min_entry_origin)_$(im.scaling)_$(im.kind_opt)_$(im.kind_queue)"
    mkpath(plot_dir)

    p1 = plot(stats.period, stats.new_demand, label="new demand")
    plot!(p1, stats.period, stats.in_queue, label="in_queue")
    plot!(p1, stats.period, stats.moved_demand, label="moved", title="Optimization")
    savefig(p1, "$(plot_dir)/optimization_stats.png")

    p2 = plot(transpose(save_queue), title="Queue Length", label="")
    savefig(p2, "$(plot_dir)/optimization_queue_length.png")

    p3 = plot(transpose(queue_period_age), title="Waiting", label="")
    savefig(p3, "$(plot_dir)/optimization_waiting.png")

    println("Computing arc utilization metrics...")
    # Precompute sums per (origin, period) - avoids repeated sum() calls
    cum_sum_op = [@views sum(im.cum_demand_od_in_period[o, :, p]) for o in 1:im.nr_nodes, p in 1:im.nr_periods]
    demand_sum_op = [@views sum(im.demand_od_in_period[o, :, p]) for o in 1:im.nr_nodes, p in 1:im.nr_periods]

    utilization_aggregated = zeros(Float64, im.nr_arcs, im.nr_minutes)
    utilization_period = zeros(Float64, im.nr_arcs, im.nr_minutes)
    println("Using $(Threads.nthreads()) threads...")
    progress = Progress(im.nr_arcs, desc="Computing utilization...")
    Threads.@threads for a in 1:im.nr_arcs
        @inbounds for t in 1:im.nr_minutes
            for (o, d, p) in shift_original[a, t]
                cum_demand_odp = im.cum_demand_od_in_period[o, d, p]
                if cum_demand_odp > 0
                    cum_sum = cum_sum_op[o, p]
                    demand_sum = demand_sum_op[o, p]
                    if cum_sum > 0 && demand_sum > 0
                        utilization_aggregated[a, t] += inflow_raw[o, p] * cum_demand_odp / cum_sum
                        utilization_period[a, t] += inflow_raw[o, p] * im.demand_od_in_period[o, d, p] / demand_sum
                    end
                end
            end
        end
        next!(progress)
    end

    # Preallocate results_queues arrays
    n_queue_rows = nr_periods * nr_nodes
    rq_datetime = Vector{DateTime}(undef, n_queue_rows)
    rq_station = Vector{String}(undef, n_queue_rows)
    rq_allowed = Vector{Int64}(undef, n_queue_rows)
    rq_moved = Vector{Int64}(undef, n_queue_rows)
    rq_queued = Vector{Int64}(undef, n_queue_rows)

    @inbounds for p in 1:nr_periods
        base = (p - 1) * nr_nodes
        for o in 1:nr_nodes
            idx = base + o
            rq_datetime[idx] = periodrange[p]
            rq_station[idx] = nodes[o]
            rq_allowed[idx] = round(Int64, inflow_raw[o, p])
            rq_moved[idx] = round(Int64, inflow_raw[o, p])
            rq_queued[idx] = round(Int64, save_queue[o, p])
        end
    end
    results_queues = DataFrame(datetime=rq_datetime, station=rq_station, allowed=rq_allowed, moved=rq_moved, queued=rq_queued)

    # Preallocate results_arcs arrays
    timesteps = collect(start_time:Minute(1):end_time)
    n_arc_rows = length(timesteps) * nr_arcs
    ra_datetime = Vector{DateTime}(undef, n_arc_rows)
    ra_connection = Vector{Int64}(undef, n_arc_rows)
    ra_line = Vector{String}(undef, n_arc_rows)
    ra_util_agg = Vector{Float64}(undef, n_arc_rows)
    ra_util_period = Vector{Float64}(undef, n_arc_rows)

    @inbounds for t in 1:length(timesteps)
        base = (t - 1) * nr_arcs
        for a in 1:nr_arcs
            idx = base + a
            ra_datetime[idx] = timesteps[t]
            ra_connection[idx] = a
            ra_line[idx] = metroarcs[a].category
            ra_util_agg[idx] = utilization_aggregated[a, t] / metroarcs[a].capacity
            ra_util_period[idx] = utilization_period[a, t] / metroarcs[a].capacity
        end
    end
    results_arcs = DataFrame(datetime=ra_datetime, connection=ra_connection, line=ra_line, utilization_aggregated=ra_util_agg, utilization_period=ra_util_period)

    sort!(results_queues, [:datetime, :station])
    sort!(results_arcs, [:datetime, :line, :connection])

    return results_queues, results_arcs, optimization_duration, build_duration, queue_period_age, infeasible_solutions
end

"""
    create_od_queue!(im, real_od_queue, config)

Purpose: Initializes a tensor that tracks passenger demand by time, origin, and destination.
Details: Uses the SAME demand data as optimization (im.demand_od_in_period) to ensure consistency.
Distributes period-level demand across minutes within each period.
"""
function create_od_queue!(im, real_od_queue, config)
    interval = config.interval_minutes
    println("Creating OD queue tensor from optimization demand...")

    # Use demand_od_in_period (already scaled in framework) to ensure consistency
    # Distribute each period's demand evenly across its minutes
    for p in 1:im.nr_periods
        # Calculate minute range for this period
        period_start_minute = (p - 1) * im.minutes_in_period + 1
        period_end_minute = min(p * im.minutes_in_period, im.nr_minutes)

        for o in 1:im.nr_nodes
            for d in 1:im.nr_nodes
                period_demand = im.demand_od_in_period[o, d, p]
                if period_demand > 0
                    # Distribute demand evenly across minutes in this period
                    demand_per_minute = period_demand / im.minutes_in_period
                    for minute in period_start_minute:period_end_minute
                        real_od_queue[minute, o, d] = demand_per_minute
                    end
                end
            end
        end
    end
end

"""
    create_entry_list!(im, real_allowed_entry, queues)

Purpose: Populates a matrix with the number of passengers allowed to enter at each station during each minute.
Details: Converts period-level entry allowances into minute-level values, ensuring the allowed entries are consistent across each period's minutes.
"""
function create_entry_list!(im, real_allowed_entry, queues)
    timesteps = start_time:Minute(1):end_time
    timestep_id = Dict(timesteps[i] => i for i in eachindex(timesteps))
    for entry in eachrow(queues)
        real_allowed_entry[timestep_id[entry.datetime]:timestep_id[entry.datetime]+im.minutes_in_period-1, d_node_id[entry.station]] .= max(0, entry.allowed)
    end
end

"""
    simulate_metro(im, queues, opt_duration, build_duration, grapharcs, kind_sim, queue_period_age, infeasible_solutions, config)

Purpose: Simulates passenger flow through the metro network given entry constraints and demand patterns.
Details: Performs a minute-by-minute simulation that:
1. Processes new passenger demand
2. Implements station entry controls based on the optimization results
3. Traces passenger movements through the network along their routes
4. Tracks queue buildups at stations
5. Calculates utilization levels of each connection

The function supports three simulation modes:
- "bound": Limits entries based on optimization results
- "inflow": Limits entries based on maximum station capacity
- "unbound": Allows all waiting passengers to enter

The function also:
- Visualizes simulation results with plots
- Calculates performance metrics (queue lengths, utilization rates, etc.)
- Logs detailed statistics to CSV files for further analysis

Returns: Two DataFrames with simulation results - station queue information and arc utilization data
"""
function simulate_metro(im, queues, opt_duration, build_duration, grapharcs, kind_sim, queue_period_age, infeasible_solutions, config)
    real_arc_use = zeros(Float64, im.nr_minutes, im.nr_arcs)
    real_queue_use = zeros(Float64, im.nr_minutes, im.nr_nodes)
    real_od_queue = zeros(Float64, im.nr_minutes, im.nr_nodes, im.nr_nodes)

    exceeded = Vector{String}

    ## load all queues into a tensor
    create_od_queue!(im, real_od_queue, config)

    ## create the list that saves the number of allowed people
    real_allowed_entry = zeros(Float64, im.nr_minutes, im.nr_nodes)
    create_entry_list!(im, real_allowed_entry, queues)

    # === DEMAND VERIFICATION: Compare optimization vs simulation initial demand ===
    sim_total_demand = sum(real_od_queue)
    opt_total_demand = sum(im.demand_od_in_period)
    demand_diff_pct = opt_total_demand > 0 ? abs(sim_total_demand - opt_total_demand) / opt_total_demand * 100 : 0.0

    # Total allowed entries from optimization
    opt_total_allowed = sum(queues.allowed) * im.minutes_in_period  # allowed is per-minute rate
    sim_total_allowed = sum(real_allowed_entry)  # Already minute-level

    println("\n=== Demand & Capacity Verification ===")
    println("  Optimization total demand:  $(round(Int, opt_total_demand))")
    println("  Simulation total demand:    $(round(Int, sim_total_demand))")
    println("  Demand difference:          $(round(demand_diff_pct, digits=2))%")
    println("  ---")
    println("  Opt total allowed entries:  $(round(Int, opt_total_allowed))")
    println("  Sim total allowed entries:  $(round(Int, sim_total_allowed))")
    println("  Expected end queue (opt):   $(round(Int, max(0, opt_total_demand - opt_total_allowed)))")
    println("  Expected end queue (sim):   $(round(Int, max(0, sim_total_demand - sim_total_allowed)))")
    if demand_diff_pct > 0.1
        println("  WARNING: Demand sources differ! This will cause queue mismatch.")
    end
    println("======================================\n")

    stats = DataFrame(
        minute=Int64[],
        new_demand=Float64[],
        in_queue=Float64[],
        moved_demand=Float64[],
    )

    ## start the flow through the network
    println("Simulating passenger flow...")
    split_flow = zeros(Float64, im.nr_nodes)  # Preallocate once
    n_minutes = size(real_arc_use, 1)
    # Track earliest minute with non-zero queue per origin (avoids O(n²) iteration)
    earliest_nonempty = ones(Int, im.nr_nodes)

    # Queue age tracking for simulation-based metrics
    birth_weighted_sum = zeros(Float64, im.nr_nodes)   # Σ(queue_size × arrival_minute)
    queue_total_tracking = zeros(Float64, im.nr_nodes) # Σ(queue_size)

    # Accumulators for avg_queue_age
    total_age_contribution = 0.0   # Σ(weighted_age) across all open (minute, origin)
    total_queue_contribution = 0.0 # Σ(queue_total) across all open (minute, origin)

    # For end_queue_age - compute last_open_minute before simulation
    last_open_minute_sim = findlast(m -> any(real_allowed_entry[m, :] .> 0), 1:n_minutes)
    if last_open_minute_sim === nothing
        last_open_minute_sim = n_minutes
    end
    sim_end_queue_age = 0.0

    @showprogress desc="Simulating minutes..." for minute in 1:n_minutes

        new_demand = @views sum(real_od_queue[minute, :, :])
        moved_demand = 0.0

        # Track new arrivals for queue age calculation
        for o in 1:im.nr_nodes
            new_arrivals = @views sum(real_od_queue[minute, o, :])
            if new_arrivals > 0
                birth_weighted_sum[o] += new_arrivals * minute
                queue_total_tracking[o] += new_arrivals
            end
        end

        for origin in eachindex(nodes)

            ## determine the number of people allowed to enter in the minute at station
            queue_sum = @views sum(real_od_queue[earliest_nonempty[origin]:minute, origin, :])

            # Skip origins with no passengers waiting
            if queue_sum <= 0
                real_queue_use[minute, origin] = 0.0
                continue
            end

            if kind_sim == "bound"
                moved_minute = min(real_allowed_entry[minute, origin], queue_sum)
            elseif kind_sim == "inflow"
                moved_minute = min(queue_sum, im.max_entry_origin)
            else  # "unbound"
                moved_minute = queue_sum
            end

            ## dispatch the ratio according to each destination
            if moved_minute > 0
                for queue_minute in earliest_nonempty[origin]:minute
                    queue_length = @views sum(real_od_queue[queue_minute, origin, :])
                    # Advance pointer if this minute is empty
                    if queue_length <= 0
                        if queue_minute == earliest_nonempty[origin]
                            earliest_nonempty[origin] = queue_minute + 1
                        end
                        continue
                    end
                    all_inflow = min(moved_minute, queue_length)
                    @views split_flow .= (real_od_queue[queue_minute, origin, :] ./ queue_length) .* all_inflow
                    moved_demand += all_inflow
                    # Update age tracking - passengers from queue_minute are being served
                    birth_weighted_sum[origin] -= all_inflow * queue_minute
                    queue_total_tracking[origin] -= all_inflow
                    @inbounds for destination in 1:im.nr_nodes
                        if split_flow[destination] > 0
                            for movement in shift_start_end[origin, destination]
                                if minute + movement[2] - 1 <= n_minutes
                                    real_arc_use[minute+movement[2]-1, movement[1]] += split_flow[destination]
                                end
                            end
                        end
                    end
                    @views real_od_queue[queue_minute, origin, :] .-= split_flow
                    if @views sum(real_od_queue[queue_minute, origin, :]) < 1e-9
                        @views real_od_queue[queue_minute, origin, :] .= 0
                        # Advance pointer if we just emptied the earliest
                        if queue_minute == earliest_nonempty[origin]
                            earliest_nonempty[origin] = queue_minute + 1
                        end
                    end
                    moved_minute -= all_inflow
                    if moved_minute < 1e-9
                        break
                    end
                end
            end
            real_queue_use[minute, origin] = @views sum(real_od_queue[earliest_nonempty[origin]:minute, origin, :])
        end

        # Accumulate queue age data for open minutes
        if any(real_allowed_entry[minute, :] .> 0)
            for o in 1:im.nr_nodes
                if queue_total_tracking[o] > 0
                    weighted_age = queue_total_tracking[o] * minute - birth_weighted_sum[o]
                    total_age_contribution += weighted_age
                    total_queue_contribution += queue_total_tracking[o]
                end
            end
        end

        # Compute end_queue_age at last open minute
        if minute == last_open_minute_sim
            end_age_sum = 0.0
            end_queue_sum = 0.0
            for o in 1:im.nr_nodes
                if queue_total_tracking[o] > 0
                    weighted_age = queue_total_tracking[o] * minute - birth_weighted_sum[o]
                    end_age_sum += weighted_age
                    end_queue_sum += queue_total_tracking[o]
                end
            end
            sim_end_queue_age = end_queue_sum > 0 ? end_age_sum / end_queue_sum : 0.0
        end

        # Clear queues if metro is closed (all stations have 0 allowed entry)
        # Note: Removed infeasible_solutions check - queue should clear during closed hours
        # regardless of whether some optimization periods were infeasible
        if all(i -> i == 0, real_allowed_entry[minute, :])
            real_od_queue[1:minute, :, :] .= 0
            earliest_nonempty .= minute + 1  # Reset pointers after bulk clear
            # Reset age tracking
            birth_weighted_sum .= 0
            queue_total_tracking .= 0
        end

        push!(stats, (
            minute=minute,
            new_demand=new_demand,
            in_queue=sum(real_queue_use[minute, :]),
            moved_demand=moved_demand,
        )
        )
    end

    # Save simulation diagnostic plots
    plot_dir = "results/plots/$(config.name)_$(Date(start_time))_$(im.safety_factor)_$(im.minutes_in_period)_$(im.past_minutes)_$(im.max_entry_origin)_$(im.min_entry_origin)_$(im.scaling)_$(im.kind_opt)_$(kind_sim)_$(im.kind_queue)"
    mkpath(plot_dir)

    p1 = plot(stats.minute, stats.new_demand, label="new demand")
    plot!(p1, stats.minute, stats.moved_demand, label="moved")
    plot!(p1, stats.minute, stats.in_queue, label="queued", title="Simulation - $(kind_sim)")
    savefig(p1, "$(plot_dir)/simulation_stats_$(kind_sim).png")


    println("Computing arc utilization...")
    @showprogress desc="Computing utilization..." for arc in axes(grapharcs, 1)
        for minute in axes(real_arc_use, 1)
            if real_arc_use[minute, arc] > 0 && grapharcs.capacity[arc] > 0
                real_arc_use[minute, arc] = real_arc_use[minute, arc] / grapharcs.capacity[arc]
            end
        end
    end

    timesteps = collect(start_time:Minute(1):end_time)
    n_timesteps = length(timesteps)

    println("Building result DataFrames...")
    # Preallocate arrays for sim_queues
    n_queue_rows = n_timesteps * nr_nodes
    queue_datetimes = Vector{DateTime}(undef, n_queue_rows)
    queue_stations = Vector{String}(undef, n_queue_rows)
    queue_queued = Vector{Int64}(undef, n_queue_rows)
    queue_allowed = Vector{Int64}(undef, n_queue_rows)

    # Preallocate arrays for sim_arcs
    n_arc_rows = n_timesteps * nr_arcs
    arc_datetimes = Vector{DateTime}(undef, n_arc_rows)
    arc_connections = Vector{Int64}(undef, n_arc_rows)
    arc_lines = Vector{String}(undef, n_arc_rows)
    arc_utilizations = Vector{Float64}(undef, n_arc_rows)

    # Fill arrays
    @showprogress desc="Building results..." for t in 1:n_timesteps
        # Queue data
        queue_base = (t - 1) * nr_nodes
        @inbounds for o in 1:nr_nodes
            idx = queue_base + o
            queue_datetimes[idx] = timesteps[t]
            queue_stations[idx] = nodes[o]
            queue_queued[idx] = round(Int64, real_queue_use[t, o])
            queue_allowed[idx] = round(Int64, real_allowed_entry[t, o])
        end
        # Arc data
        arc_base = (t - 1) * nr_arcs
        @inbounds for a in 1:nr_arcs
            idx = arc_base + a
            arc_datetimes[idx] = timesteps[t]
            arc_connections[idx] = a
            arc_lines[idx] = metroarcs[a].category
            arc_utilizations[idx] = real_arc_use[t, a]
        end
    end

    # Create DataFrames from preallocated arrays
    sim_queues = DataFrame(datetime=queue_datetimes, station=queue_stations, queued=queue_queued, allowed=queue_allowed)
    sim_arcs = DataFrame(datetime=arc_datetimes, connection=arc_connections, line=arc_lines, utilization=arc_utilizations)

    logfile_name = if kind_sim == "unbound"
        "logfile_$(config.name)_$(start_time)_unbound.csv"
    else
        "logfile_$(config.name)_$(start_time)_$(im.minutes_in_period).csv"
    end
    if isfile(logfile_name)
        logfile = CSV.read(logfile_name, DataFrame)
    else
        logfile = DataFrame(
            timestamp=DateTime[],
            safety=Float64[],
            intervall=String[],
            period_length=Int64[],
            past_minutes=Int64[],
            max_enter=Float64[],
            min_enter=Float64[],
            scaling=Float64[],
            kind_optimization=String[],
            kind_simulation=String[],
            kind_queue=String[],
            infeasible=Int64[],
            avg_duration=Float64[],
            total_duration=Float64[],
            max_duration=Float64[],
            avg_build_time=Float64[],
            avg_queue=Float64[],
            end_queue=Float64[],
            avg_queue_age=Float64[],
            end_queue_age=Float64[],
            total_demand=Float64[],
            people_moved=Float64[],
            avg_utilization=Float64[],
            max_utilization=Float64[],
            exceeded_minutes=Float64[],
            exceeded_mean=Float64[],
            exceeded_median=Float64[],
            exceeded_080quant=Float64[],
            exceeded_090quant=Float64[],
            safety_minutes=Float64[],
            safety_mean=Float64[],
            safety_median=Float64[],
            safety_080quant=Float64[],
            safety_090quant=Float64[],)
    end

    # Filter to only open hours for accurate averages
    open_minutes = findall(m -> any(real_allowed_entry[m, :] .> 0), 1:n_minutes)
    n_open_minutes = length(open_minutes)

    # Calculate averages only over operating hours
    avg_queue_open = sum(real_queue_use[open_minutes, :]) / (n_open_minutes * nr_nodes)
    end_queue_total = sum(real_queue_use[last_open_minute_sim, :])

    # Compute queue ages - set to 0 if no queue exists (avoid misleading values)
    sim_avg_queue_age = (avg_queue_open > 0.01 && total_queue_contribution > 0) ? total_age_contribution / total_queue_contribution : 0.0
    sim_end_queue_age = end_queue_total > 0.01 ? sim_end_queue_age : 0.0
    open_arc_rows = filter(row -> hour(row.datetime) ∉ config.closed_hours, eachrow(sim_arcs))
    avg_util_open = isempty(open_arc_rows) ? 0.0 : sum(r.utilization for r in open_arc_rows) / length(open_arc_rows)

    ex = length(filter(x -> x > 1, sim_arcs.utilization)) == 0 ? [0] : filter(x -> x > 1, sim_arcs.utilization)
    sf = length(filter(x -> x > im.safety_factor, sim_arcs.utilization)) == 0 ? [0] : filter(x -> x > im.safety_factor, sim_arcs.utilization)

    # === QUEUE VERIFICATION: Compare optimization vs simulation queues ===
    # Find last OPEN period for comparison
    all_periods = sort(unique(queues.datetime))
    last_open_period_dt = nothing
    for dt in reverse(all_periods)
        if hour(dt) ∉ config.closed_hours
            last_open_period_dt = dt
            break
        end
    end

    if last_open_period_dt !== nothing
        opt_final_queues = filter(row -> row.datetime == last_open_period_dt, queues)
        opt_total_queue = sum(opt_final_queues.queued)
        sim_total_queue = sum(real_queue_use[last_open_minute_sim, :])
        queue_diff_pct = opt_total_queue > 0 ? abs(sim_total_queue - opt_total_queue) / opt_total_queue * 100 : 0.0

        if queue_diff_pct > 1.0
            println("\n  WARNING: Queue mismatch of $(round(queue_diff_pct, digits=1))% between optimization and simulation (likely due to rounding).\n")
        end
    end

    push!(logfile, (
        timestamp=now(),
        safety=im.safety_factor,
        intervall=string(start_time) * " to " * string(end_time),
        period_length=im.minutes_in_period,
        past_minutes=im.past_minutes,
        max_enter=im.max_entry_origin,
        min_enter=im.min_entry_origin,
        scaling=im.scaling,
        kind_optimization=im.kind_opt,
        kind_simulation=kind_sim,
        kind_queue=im.kind_queue,
        infeasible=infeasible_solutions,
        avg_duration=sum(opt_duration) / length(opt_duration),
        total_duration=sum(opt_duration),
        max_duration=maximum(opt_duration),
        avg_build_time=sum(build_duration) / length(build_duration),
        avg_queue=avg_queue_open,
        end_queue=end_queue_total,
        avg_queue_age=sim_avg_queue_age,
        end_queue_age=sim_end_queue_age,
        total_demand=sum(stats.new_demand),
        people_moved=sum(stats.moved_demand),
        avg_utilization=avg_util_open,
        max_utilization=maximum(sim_arcs.utilization),
        exceeded_minutes=length(filter(x -> x > 1, sim_arcs.utilization)),
        exceeded_mean=mean(filter(x -> x > 1, sim_arcs.utilization)),
        exceeded_median=quantile(ex, 0.5),
        exceeded_080quant=quantile(ex, 0.80),
        exceeded_090quant=quantile(ex, 0.90),
        safety_minutes=length(filter(x -> x > im.safety_factor, sim_arcs.utilization)),
        safety_mean=mean(filter(x -> x > im.safety_factor, sim_arcs.utilization)),
        safety_median=quantile(sf, 0.5),
        safety_080quant=quantile(sf, 0.80),
        safety_090quant=quantile(sf, 0.90),
    ))

    # Save logfile to both root folder and results folder
    CSV.write(logfile_name, logfile)
    CSV.write("results/$logfile_name", logfile)

    CSV.write("results/queues/sim_queues_$(config.name)_$(Date(start_time))_$(im.safety_factor)_$(im.minutes_in_period)_$(im.past_minutes)_$(im.max_entry_origin)_$(im.min_entry_origin)_$(im.scaling)_$(im.kind_opt)_$(kind_sim)_$(im.kind_queue).csv",
        sim_queues
    )

    CSV.write("results/arcs/sim_arcs_$(config.name)_$(Date(start_time))_$(im.safety_factor)_$(im.minutes_in_period)_$(im.past_minutes)_$(im.max_entry_origin)_$(im.min_entry_origin)_$(im.scaling)_$(im.kind_opt)_$(kind_sim)_$(im.kind_queue).csv",
        sim_arcs
    )
    return sim_queues, sim_arcs
end

function heuristic_adding_queues(im)
    # Initialization of variables
    remaining_queue = copy(im.demand_od_in_period) # Initialize a new variable to store the remaining queue data
    inflow_raw = zeros(Float64, im.nr_nodes, im.nr_periods) .= im.min_entry_origin  # Initialize an array to store the values of the X variable
    optimization_duration = zeros(Float64, im.nr_periods) # Initialize optimization duration vector
    queue_period_age = zeros(Int64,im.nr_nodes,im.nr_periods) .= 0 # Computes the age of each queue per period
    save_queue = zeros(Int64,nr_nodes,im.nr_periods) .= 1 # Saves the queue for the output

    stats = DataFrame(
        period = Int64[],
        new_demand = Float64[],
        in_queue = Float64[],
        moved_demand = Float64[],
    )

    infeasible_solutions = 0

    @showprogress for fix_period in 1:im.nr_periods
        
        println()
        println("Running period ", fix_period)
        im.cum_demand_od_in_period[:,:,fix_period] .= sum(remaining_queue[:,:,1:fix_period],dims=3)

        # Prepare the upper and lower bound of the observed time slice
        lower_period::Int64 = max(1,fix_period - floor((60+im.minutes_in_period)/im.minutes_in_period))
        upper_period::Int64 = min(fix_period,im.nr_periods)

        println(sum(sum(remaining_queue[:,:,1:fix_period],dims=3)))

        for o in 1:nr_nodes
            for queue_period in 1:fix_period
                if sum(remaining_queue[o,:,queue_period]) >= 1
                    queue_period_age[o,fix_period] = fix_period - queue_period
                    break
                end
            end
        end

        println("Queue Age")
        println(queue_period_age[:,fix_period])
        println("Queue Length")
        println(round.(Int64,sum(im.cum_demand_od_in_period[:,:,fix_period],dims=2)))
        println("Total Queue Length")
        println(round.(Int64,sum(im.cum_demand_od_in_period[:,:,fix_period])))

        model,X = build_restricted_optimization_model(im,fix_period,queue_period_age,inflow_raw)

        for o in 1:im.nr_nodes
            if fix_period > 1
                for p in lower_period:fix_period-1
                    if inflow_raw[o,p] > 0.0001
                        fix(X[o,p],inflow_raw[o,p]; force = true)
                    else
                        fix(X[o,p],0; force = true)
                    end
                end
            end
        end

        if im.closed_period[fix_period] == false

            optimization_duration[fix_period] = @elapsed optimize!(model)

            if is_solved_and_feasible(model) == false
                infeasible_solutions += 1
            end

            for o in 1:im.nr_nodes
                for p in fix_period:upper_period
                    if value.(X)[o,p] > 0.0001
                        inflow_raw[o,p] = value.(X)[o,p]
                    else
                        inflow_raw[o,p] = 0.00
                    end
                end
            end

        else

            inflow_raw[:,fix_period] .= 0
            
        end

        
        demand_fulfiled = inflow_raw[:,fix_period] .* im.minutes_in_period 
        println("Dispatch")
        println(round.(demand_fulfiled))
        println("Total Dispatch")
        println(sum(demand_fulfiled))

        for o in 1:im.nr_nodes 
            for p in 1:fix_period
                current_ratio = (remaining_queue[o,:,p]/sum(remaining_queue[o,:,p]))
                current_ratio .= ifelse.(isnan.(current_ratio), 0, current_ratio)

                while sum(remaining_queue[o,:,p]) >= 1 && demand_fulfiled[o] >= 1
                    remaining_queue[o,:,p] .-=  im.minutes_in_period .* current_ratio
                    demand_fulfiled[o] -= im.minutes_in_period
                end

                if sum(remaining_queue[o,:,p]) <= 1
                    remaining_queue[o,:,p] .= 0
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
            remaining_queue[:,:,1:fix_period] .= 0
        end
        for o in 1:im.nr_nodes
            save_queue[o,fix_period] = round(Int64,sum(remaining_queue[o,:,1:fix_period]))
        end


        push!(stats,(
            period = fix_period,
            new_demand = sum(im.demand_od_in_period[:,:,fix_period]),
            in_queue = sum(save_queue[:,fix_period]),
            moved_demand = sum(inflow_raw[:,fix_period] .* minutes_in_period),
            )
        )
    end

    plot(stats.period,stats.new_demand,label="new demand")
    plot!(stats.period,stats.in_queue,label="in_queue")
    display(plot!(stats.period,stats.moved_demand,label="moved",title="Optimization"))

    display(plot(transpose(save_queue),title="Queue Length", label=""))
    display(plot(transpose(queue_period_age),title="Waiting", label=""))
    
    results_queues = DataFrame(datetime = DateTime[],station=String[],allowed=Int64[],moved=Int64[],queued=Int64[])
    results_arcs   = DataFrame(datetime = DateTime[],connection=Int64[],line=String[],utilization_aggregated=Float64[],utilization_period=Float64[])

    utilization_aggregated = zeros(Float64,im.nr_arcs,im.nr_minutes)
    utilization_period = zeros(Float64,im.nr_arcs,im.nr_minutes)
    for a in 1:im.nr_arcs
        for t in 1:im.nr_minutes
            for (o,d,p) in shift_original[a,t]
                if im.cum_demand_od_in_period[o,d,p] > 0
                    utilization_aggregated[a,t] += sum(inflow_raw[o,p] * im.cum_demand_od_in_period[o,d,p]/sum(im.cum_demand_od_in_period[o,:,p]))
                    utilization_period[a,t] += sum(inflow_raw[o,p] * im.demand_od_in_period[o,d,p]/sum(im.demand_od_in_period[o,:,p]))
                end
            end
        end
    end

    for p in 1:nr_periods
        for o in 1:nr_nodes
            push!(results_queues,(
                datetime = periodrange[p],
                station=nodes[o],
                allowed=round(inflow_raw[o,p]),
                moved=round(inflow_raw[o,p]),
                queued=round(save_queue[o,p])))
        end
    end

    for t in 1:length(start_time:Minute(1):end_time)
        for a in 1:nr_arcs
            push!(results_arcs,(
                    datetime = (start_time:Minute(1):end_time)[t],
                    connection=a,
                    line=metroarcs[a].category,
                    utilization_aggregated=utilization_aggregated[a,t]/metroarcs[a].capacity,
                    utilization_period=utilization_period[a,t]/metroarcs[a].capacity
                )
            )
        end
    end

    sort!(results_queues,[:datetime,:station])
    sort!(results_arcs,[:datetime,:line,:connection])

    return results_queues, results_arcs, optimization_duration, queue_period_age, infeasible_solutions
end

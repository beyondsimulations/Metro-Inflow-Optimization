function heuristic_adding_queues()
    # Initialization of variables
    demand_od_heuristic = copy(demand_od)  # Initialize a new variable to store the aggregated demand data
    remaining_queue = copy(demand_od) .= 0  # Initialize a new variable to store the remaining queue data
    inflow_raw = zeros(Float64, nr_nodes, nr_periods)  # Initialize an array to store the values of the X variable
    optimization_duration = zeros(Float64, nr_periods) # Initialize optimization duration vector

    for fix_period in 1:nr_periods
        println("Running period ", fix_period)
        instance = build_model_instance(demand_od_heuristic,demand_od)

        model,X = build_restricted_optimization_model(instance,past_minutes,minutes_in_period,fix_period)

        for o in 1:nr_nodes
            if fix_period > 1
                for p in 1:fix_period-1
                    if inflow_raw[o,p] > 0.01
                        fix(X[o,p],inflow_raw[o,p]; force = true)
                    else
                        fix(X[o,p],0; force = true)
                    end
                end
            end
            
            if fix_period+ceil(Int,past_minutes/minutes_in_period) <= nr_periods
                for p in fix_period+ceil(Int,past_minutes/minutes_in_period):nr_periods
                    fix(X[o,p],0; force = true)
                end
            end
        end

        optimization_duration[fix_period] = @elapsed optimize!(model)

        for v in eachindex(value.(X))
            if value.(X)[v] > 0.1
                inflow_raw[v] = floor(value.(X)[v],digits=2)
            else
                inflow_raw[v] = 0.00
            end
        end
    
        demand_fulfiled = zeros(Float64,nr_nodes,nr_nodes)
        for o in 1:nr_nodes
            for d in 1:nr_nodes
                if demand_od_heuristic[o,d,fix_period] > 0
                    demand_fulfiled[o,d] = minutes_in_period * inflow_raw[o,fix_period] * (demand_od_heuristic[o,d,fix_period]/sum(demand_od_heuristic[o,:,fix_period]))
                end
            end
        end
        remaining_queue[:,:,fix_period] .= demand_od_heuristic[:,:,fix_period] .- demand_fulfiled

        if fix_period < nr_periods

            demand_od_heuristic[:,:,fix_period+1] .+= round.(remaining_queue[:,:,fix_period],digits=2)
            demand_od_heuristic[:,:,fix_period]   .= demand_fulfiled

        end
    end

    results_queues = DataFrame(datetime = DateTime[],station=String[],allowed=Int64[],moved=Int64[],queued=Int64[])
    results_arcs   = DataFrame(datetime = DateTime[],connection=Int64[],line=String[],utilization_aggregated=Float64[],utilization_period=Float64[])

    utilization_aggregated = zeros(Float64,nr_arcs,nr_minutes)
    utilization_period = zeros(Float64,nr_arcs,nr_minutes)
    for a in 1:nr_arcs
        for t in 1:nr_minutes
            for (o,d,p) in shift_original[a,t]
                if demand_od_heuristic[o,d,p] > 0
                    utilization_aggregated[a,t] += sum(inflow_raw[o,p] * demand_od_heuristic[o,d,p]/sum(demand_od_heuristic[o,:,p]))
                    utilization_period[a,t] += sum(inflow_raw[o,p] * demand_od[o,d,p]/sum(demand_od[o,:,p]))
                end
            end
        end
    end

    remaining_queue = sum(remaining_queue,dims=2)[:,1,:]

    for p in 1:nr_periods-1
        for o in 1:nr_nodes
            push!(results_queues,(
                datetime = periodrange[p],
                station=nodes[o],
                allowed=round(inflow_raw[o,p]),
                moved=round(inflow_raw[o,p]),
                queued=round(remaining_queue[o,p])))
        end
    end

    for t in 1:nr_minutes-1
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

    return results_queues, results_arcs, optimization_duration
end

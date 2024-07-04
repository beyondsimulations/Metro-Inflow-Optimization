function create_od_queue!(im,real_od_queue)
    demand = load_demand()
    timesteps = start_time:Minute(1):end_time
    timestep_id = Dict(timesteps[i] => i for i in eachindex(timesteps))
    for data in eachrow(demand)
        for minutestep in data.datetime:Minute(1):data.datetime+Minute(14)
            real_od_queue[timestep_id[minutestep],d_node_id[data.origin],d_node_id[data.destination]] = (data.value * im.scaling)/15
        end
    end

end

function create_entry_list!(im,real_allowed_entry,queues)
    timesteps = start_time:Minute(1):end_time
    timestep_id = Dict(timesteps[i] => i for i in eachindex(timesteps))
    for entry in eachrow(queues)
        real_allowed_entry[timestep_id[entry.datetime]:timestep_id[entry.datetime]+im.minutes_in_period-1,d_node_id[entry.station]] .= max(0,entry.allowed)
    end
end

function simulate_metro(im,queues,opt_duration,grapharcs,kind_sim,queue_period_age,infeasible_solutions)

    real_arc_use = zeros(Float64,im.nr_minutes,im.nr_arcs)
    real_queue_use = zeros(Float64,im.nr_minutes,im.nr_nodes)
    real_od_queue = zeros(Float64,im.nr_minutes,im.nr_nodes,im.nr_nodes)

    exceeded = Vector{String}

    ## load all queues into a tensor
    create_od_queue!(im,real_od_queue)

    ## create the list that saves the number of allowed people
    real_allowed_entry = zeros(Float64,im.nr_minutes,im.nr_nodes)
    create_entry_list!(im,real_allowed_entry,queues)

    stats = DataFrame(
        minute = Int64[],
        new_demand = Float64[],
        in_queue = Float64[],
        moved_demand = Float64[],
    )

    ## start the flow through the network
    @showprogress for minute in 1:size(real_arc_use,1)

        new_demand = sum(real_od_queue[minute,:,:])
        moved_demand = 0

        for origin in eachindex(nodes)

        ## determine the number of people allowed to enter in the minute at station
            if kind_sim == "bound"
                moved_minute = min(real_allowed_entry[minute,origin],sum(real_od_queue[1:minute,origin,:]))

            elseif kind_sim == "inflow"
                moved_minute = min(sum(real_od_queue[1:minute,origin,:]),im.max_entry_origin)

            else kind_sim == "unbound"
                moved_minute = sum(real_od_queue[1:minute,origin,:])

            end

            ## dispatch the ratio according to each destination
            if moved_minute > 0
                for queue_minute in 1:minute
                    if any(real_od_queue[queue_minute,origin,:] .> 0)
                        queue_length = sum(real_od_queue[queue_minute,origin,:])
                        all_inflow = min(moved_minute,queue_length)
                        split_flow = (real_od_queue[queue_minute,origin,:]/queue_length) .* all_inflow
                        moved_demand += all_inflow
                        for destination in 1:im.nr_nodes
                            if split_flow[destination] > 0
                                for movement in shift_start_end[origin,destination]
                                    if minute+movement[2]-1 <= size(real_arc_use,1)
                                        real_arc_use[minute+movement[2]-1,movement[1]] += split_flow[destination]
                                    end
                                end
                            end
                        end
                        real_od_queue[queue_minute,origin,:] .-= split_flow
                        if sum(real_od_queue[queue_minute,origin,:]) <= 0.1
                            real_od_queue[queue_minute,origin,:] .= 0
                        end
                        moved_minute -= all_inflow
                        if moved_minute <= 0.1
                            break
                        end
                    end
                end
            end
            real_queue_use[minute,origin] = sum(real_od_queue[1:minute,origin,:])

            if all(i->i==0,real_allowed_entry[minute,:]) && infeasible_solutions == 0
                real_od_queue[1:minute,:,:] .= 0
            end
        end

        push!(stats,(
            minute = minute,
            new_demand = new_demand,
            in_queue = sum(real_queue_use[minute,:]),
            moved_demand = moved_demand,
            )
        )
    end

    plot(stats.minute,stats.new_demand,label="new demand")
    display(plot!(stats.minute,stats.moved_demand,label="moved",title="Simulation - A"))
    display(plot!(stats.minute,stats.in_queue,label="queued",title="Simulation - B"))
  

    for arc in axes(grapharcs,1)
        for minute in axes(real_arc_use,1)
            if real_arc_use[minute,arc] > 0 && grapharcs.capacity[arc] > 0
                real_arc_use[minute,arc] = real_arc_use[minute,arc] / grapharcs.capacity[arc]
            end
        end
    end

    sim_queues = DataFrame(datetime = DateTime[],station=String[],queued=Int64[],allowed=Int64[])
    sim_arcs   = DataFrame(datetime = DateTime[],connection=Int64[],line=String[],utilization=Float64[])

    timesteps = start_time:Minute(1):end_time

    for t in 1:length(start_time:Minute(1):end_time)
        for o in 1:nr_nodes
            push!(sim_queues,(
                datetime = timesteps[t],
                station=nodes[o],
                queued=round(real_queue_use[t,o]),
                allowed=round(real_allowed_entry[t,o]))
            )
        end
        for a in 1:nr_arcs
            push!(sim_arcs,(
                    datetime = timesteps[t],
                    connection=a,
                    line=metroarcs[a].category,
                    utilization=real_arc_use[t,a],
                )
            )
        end
    end

    if isfile("logfile.csv")
        logfile = CSV.read("logfile.csv", DataFrame)
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
            safety_090quant=Float64[],
            
            )
    end

    ex = length(filter(x -> x > 1, sim_arcs.utilization)) == 0 ? [0] : filter(x -> x > 1, sim_arcs.utilization)
    sf = length(filter(x -> x > im.safety_factor, sim_arcs.utilization)) == 0 ? [0] : filter(x -> x > im.safety_factor, sim_arcs.utilization)

    push!(logfile, (
        timestamp = now(),
        safety = im.safety_factor,
        intervall = string(start_time) * " to " * string(end_time),
        period_length = im.minutes_in_period,
        past_minutes = im.past_minutes,
        max_enter= im.min_entry_origin,
        min_enter= im.max_entry_origin,
        scaling=im.scaling,
        kind_optimization=im.kind_opt,
        kind_simulation=kind_sim,
        kind_queue=im.kind_queue,
        infeasible=infeasible_solutions,
        avg_duration = sum(opt_duration)/length(opt_duration),
        total_duration =  sum(opt_duration),
        avg_queue = sum(sim_queues.queued)/nrow(sim_queues),
        end_queue = sum(real_queue_use[end,:]),
        avg_queue_age=(sum(queue_period_age)/(size(queue_period_age,1)*size(queue_period_age,2)))*im.minutes_in_period,
        end_queue_age=(sum(queue_period_age[end,:])/(size(queue_period_age,2)))*im.minutes_in_period,
        total_demand=sum(stats.new_demand),
        people_moved=sum(stats.moved_demand),
        avg_utilization = sum(sim_arcs.utilization)/nrow(sim_arcs),
        max_utilization = maximum(sim_arcs.utilization),
        exceeded_minutes = length(filter(x -> x > 1, sim_arcs.utilization)),
        exceeded_mean = mean(filter(x -> x > 1, sim_arcs.utilization)),
        exceeded_median = quantile(ex,0.5),
        exceeded_080quant = quantile(ex,0.80),
        exceeded_090quant = quantile(ex,0.90),
        safety_minutes = length(filter(x -> x > im.safety_factor, sim_arcs.utilization)),
        safety_mean = mean(filter(x -> x >  im.safety_factor, sim_arcs.utilization)),
        safety_median = quantile(sf,0.5),
        safety_080quant = quantile(sf,0.80),
        safety_090quant = quantile(sf,0.90),
    ))

    CSV.write("logfile.csv",logfile)
    return sim_queues,sim_arcs
end
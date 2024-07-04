# load packages
using Pkg
Pkg.activate("metroflow")

using DataFrames
using CSV
using DataStructures
using JuMP
using HiGHS
using Plots
using StatsPlots
using Dates
using Measures

# parameters for the modification
safety = 0.9            # safety factor substracted from the arc capacity
minutes_in_step = 15    # minutes in each timestep
horizon = 4             # "forecast" horizon (in steps)
max_enter = 250         # number of maximal entries per queue per station
max_change = 30         # maximal increase in the entry to the next period
rate_closed = 600       # removal rate per period during closed metro hours
closed = Hour.(3:5)     # closed hours of the metro
scaling = 1.0           # scaling of the metro queue

# daterange
daterange = Date("2022-11-27"):Date("2022-11-28")

# load the data
grapharcs = CSV.read("data_metro/metroarcs.csv", DataFrame)
grapharcs.traveltime = ceil.(Int64,grapharcs.traveltime)
grapharcs.utilization .= 0.0
sort!(grapharcs, :category)

# include the files with additional functions
include("functions/djikstra.jl")
include("functions/mathematical_model.jl")

# prepare the graph
nodes = sort(unique(vcat(grapharcs.origin,grapharcs.destination)))
nodeid = Dict([nodes[i] => i for i in eachindex(nodes)])
arcid = Dict((nodeid[grapharcs.origin[i]],nodeid[grapharcs.destination[i]]) => i for i in axes(grapharcs,1))
idarc = Dict(i => (nodeid[grapharcs.origin[i]],nodeid[grapharcs.destination[i]]) for i in axes(grapharcs,1))

distance, previous_nodes = djisktra(nodes,nodeid,grapharcs)

function create_graph(
    distance,
    previous_nodes,
    grapharcs,
    minutes_in_step
    )
    nr_nodes = size(previous_nodes,1)
    nr_arcs = size(grapharcs,1)
    on_path = zeros(Bool,nr_nodes,nr_nodes,nr_arcs,ceil(Int64,maximum(distance)+minutes_in_step*horizon),ceil(Int64,maximum(distance)))
    for origin in 1:nr_nodes
        for destination in 1:nr_nodes
            previous = previous_nodes[origin,destination]
            current_destination = destination
            destination_time = ceil(Int64,distance[origin,destination])
            while previous != 0
                start_time = ceil(Int64,destination_time - grapharcs.traveltime[arcid[(previous,current_destination)]])
                for startminute in 1:ceil(Int64,maximum(distance))
                    for time in start_time+startminute:destination_time+startminute-1
                        if time <= size(on_path,5)
                            on_path[origin,destination,arcid[(previous,current_destination)],time,startminute] = true
                        end
                    end
                end
                current_destination = previous
                previous = previous_nodes[origin,current_destination]
                destination_time = start_time
            end
        end
    end
    return on_path
end

function new_arrivals!(
    nodeid,
    demand,
    step,
    new_queue::Matrix{Float64}
    )
    for movement in axes(demand,1)
        if demand.datetime[movement] == step
            new_queue[nodeid[demand.origin[movement]],nodeid[demand.destination[movement]]] += demand.value[movement] * scaling
        end
    end
end

function new_arrivals!(
    nodeid,
    demand,
    step,
    dt,
    new_queue::Array{Float64,3}
    )
    for movement in axes(demand,1)
        if demand.datetime[movement] == dt
            new_queue[nodeid[demand.origin[movement]],nodeid[demand.destination[movement]],step] += demand.value[movement] * scaling
        end
    end
end

function dispatch_queues!(
    station_moved,
    queue
    )
    total_queue = sum(queue,dims=2)
    for o in axes(queue,1)
        for d in axes(queue,2)
            if queue[o,d] > 0
                queue[o,d] = max(0,round(queue[o,d] - station_moved[o] * (queue[o,d]/total_queue[o])))
            end
        end
    end
end

function arc_utilization!(
    grapharcs,
    arc_capacity,
    )
    grapharcs.utilization = (grapharcs.capacity .- vec(sum(arc_capacity[1:minutes_in_step,:],dims=1) ./ minutes_in_step)) ./ grapharcs.capacity
end

function arcs_lower_capacity!(station_moved,arc_capacity,queue,on_path)
    total_queue = sum(queue,dims=2)
    for o in axes(on_path,1)
        for d in axes(on_path,2)
            for a in axes(on_path,3)
                for t in 1:minutes_in_step
                    if on_path[o,d,a,t] == true
                        if queue[o,d] > 0
                            arc_capacity[t,a] -= station_moved[o]*(queue[o,d]/(total_queue[o]))
                        end
                    end
                end
            end
        end
    end
end

function arc_capacity!(arc_capacity,grapharcs,minutes_in_step)
    arc_capacity[1:size(arc_capacity,1)-minutes_in_step,:] .= arc_capacity[minutes_in_step+1:size(arc_capacity,1),:]
    for row in size(arc_capacity,1)-minutes_in_step+1:size(arc_capacity,1)
        arc_capacity[row,:] .= grapharcs.capacity
    end
end

function timestep_optimization(
    date,
    step,
    nodeid,
    nodes,
    demand,
    queue,
    grapharcs,
    on_path,
    minutes_in_step,
    safety,
    max_enter,
    max_change,
    arc_capacity,
    station_allowed,
    horizon
    )
    println("Starting timestep $step")
        new_arrivals!(nodeid,demand,step,queue)
        arc_capacity!(arc_capacity,grapharcs,minutes_in_step)
        if !(DateTime(date)+minimum(closed) <= step <= DateTime(date)+maximum(closed))
            inflow_model = create_model()
            X = build_mainmodel(
                inflow_model,
                nodes,on_path,
                max_enter,
                queue,
                minutes_in_step,
                horizon)
            update_metro_model(
                inflow_model,
                X,
                nodes,
                queue,
                on_path,
                safety,
                arc_capacity,
                grapharcs,
                max_change,
                station_allowed,
                )
            JuMP.optimize!(inflow_model)
            solution_summary(inflow_model)
            station_moved = value.(X[:,1])
            arcs_lower_capacity!(station_moved,arc_capacity,queue,on_path)
            station_allowed .= round.(Int64,station_moved)
        else
            station_moved = fill(rate_closed,length(nodes))
            station_allowed .= 0
        end

        arc_utilization!(grapharcs,arc_capacity)

        station_queue = sum(queue,dims=2)

        dispatch_queues!(station_moved*minutes_in_step,queue)

        return station_moved,station_queue
    end

function metro_inflow(
    daterange,
    nodeid,
    nodes,
    grapharcs,
    safety,
    minutes_in_step,
    max_enter,
    max_change,
    previous_nodes,
    distance,
    horizon)

    
    on_path = create_graph(distance,previous_nodes,grapharcs,minutes_in_step)
    queue = zeros(Float64,length(nodeid),length(nodeid))
    results_queues = DataFrame(datetime = DateTime[],station=String[],allowed=Int64[],moved=Int64[],queued=Int64[])
    results_arcs   = DataFrame(datetime = DateTime[],connection=Int64[],line=String[],utilization=Float64[])
    arc_capacity = zeros(Float64,ceil(Int64,maximum(distance)+minutes_in_step),size(grapharcs,1))
    station_allowed = zeros(Int64,length(nodes))

    for row in axes(arc_capacity,1)
        arc_capacity[row,:] .= grapharcs.capacity
    end

    date = minimum(daterange)

    demand = CSV.read("data_demand/OD_$(date)v3.csv", DataFrame)

    timesteps = DateTime(date):Minute(minutes_in_step):DateTime(date)+Day(1)-Minute(minutes_in_step)

    for step in timesteps
        
        station_moved,
        station_queue = timestep_optimization(
            date,
            step,
            nodeid,
            nodes,
            demand,
            queue,
            grapharcs,
            on_path,
            minutes_in_step,
            safety,
            max_enter,
            max_change,
            arc_capacity,
            station_allowed,
            horizon)

    end

    for date in daterange

        demand = CSV.read("data_demand/OD_$(date)v3.csv", DataFrame)

        timesteps = DateTime(date):Minute(minutes_in_step):DateTime(date)+Day(1)-Minute(minutes_in_step)

        for step in timesteps
            
            station_moved,
            station_queue = timestep_optimization(
                date,
                step,
                nodeid,
                nodes,
                demand,
                queue,
                grapharcs,
                on_path,
                minutes_in_step,
                safety,
                max_enter,
                max_change,
                arc_capacity,
                station_allowed,
                horizon)

                
            for station in eachindex(nodes)
                push!(results_queues,(
                    datetime = step,
                    station=nodes[station],
                    allowed=round(station_allowed[station]),
                    moved=round(station_moved[station]),
                    queued=round(station_queue[station])))
            end

            for arc in axes(grapharcs,1)
                push!(results_arcs,(
                    datetime = step,
                    connection=arc,
                    line=grapharcs.category[arc],
                    utilization=grapharcs.utilization[arc]))
            end
        end
        sort!(results_queues,[:datetime,:station])
        sort!(results_arcs,[:datetime,:line,:connection])
        CSV.write("results/qeues_$horizon",results_queues)
        CSV.write("results/arcs_$horizon",results_arcs)
    end
    return results_queues, results_arcs
end


# Optimization
results_queues, results_arcs = metro_inflow(
    daterange,
    nodeid,
    nodes,
    grapharcs,
    safety,
    minutes_in_step,
    max_enter,
    max_change,
    previous_nodes,
    distance,
    horizon
    )

# Plot Optimization results
function plot_optimization!(results_queues,results_arcs)
    Plots.scalefontsizes()
    Plots.scalefontsizes(2.2)
    steprange = unique(results_queues.datetime)
    max_queue = maximum(results_queues.queued)*1.1
    max_entry = maximum(results_queues.allowed)*1.1
    gp_results_queues = groupby(results_queues,:datetime)
    gp_results_arcs = groupby(results_arcs,:datetime)

    get_group(gdf, keys...) = gdf[(keys...,)]

    anim = @animate for dt in steprange
        begin
            new_plot = bar(
                get_group(gp_results_queues, dt).queued,
                size = (1920, 1080),
                title="Queue at timestep $dt",
                ylims=(0,max_queue*1.2),
                label="People in Queue",
                legend=:topright,
                colour=:steelblue2,
                margin=15mm,
            )
            bar!(
                get_group(gp_results_queues, dt).allowed .* minutes_in_step,
                label="People allowed to Enter",
                xticks=(1:length(nodes), 
                get_group(gp_results_queues, dt).station),
                xrotation=30,
                xtickfontsize=10,
                colour=:orange1,
            )
        end
    end
    gif(anim, "visuals/queues_fp04_$horizon.gif", fps = 4)

    anim = @animate for dt in steprange
        begin
            new_plot = bar(
                get_group(gp_results_queues, dt).allowed,
                label="People allowed to Enter per Minute",
                size = (1920, 1080),
                title="Entry per Minute at timestep $dt",
                ylims=(0,max_entry*1.2),
                legend=:topright,
                margin=15mm,
                xticks=(1:length(nodes),get_group(gp_results_queues, dt).station),
                xrotation=30,
                xtickfontsize=10,
                colour=:orange1,
            )
        end
    end
    gif(anim, "visuals/entry_fps04_$horizon.gif", fps = 4)

    anim = @animate for dt in steprange
        begin
            new_plot = bar(
                get_group(gp_results_arcs, dt).utilization,
                group = get_group(gp_results_arcs, dt).line,
                size = (1920, 1080),
                title="Theoretical Metroarc Usage at timestep $dt",
                #xticks=(1:length(idarc), [idarc[i] for i = 1:length(idarc)]),
                ylims=(0,1.5),
                ylabel = "Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=15mm,
            )
            hline!([safety],label="Restriction",linewidth=2,colour=:grey60,)
        end
    end
    gif(anim, "visuals/utilization_fps04_$horizon.gif", fps = 4)
end

plot_optimization!(results_queues,results_arcs)

function create_od_queue(nodeid,daterange,minutes_in_step)
    real_od_queue = zeros(Float64,length(nodeid),length(nodeid),length(daterange)*96)
    for datestep in eachindex(daterange)
        timesteps = DateTime(daterange[datestep]):Minute(minutes_in_step):DateTime(daterange[datestep])+Day(1)-Minute(minutes_in_step)
        demand = CSV.read("data_demand/OD_$(daterange[datestep])v3.csv", DataFrame)
        for step in eachindex(timesteps)
            all_step = (datestep-1)*length(timesteps) + step
            dt = timesteps[step]
            new_arrivals!(nodeid,demand,all_step,dt,real_od_queue)
        end
    end
    return real_od_queue
end

# Simulation
function simulate_metro(results_queues,nodeid,minutes_in_step,daterange,grapharcs,kind)
    
    real_arc_use = zeros(Float64,length(unique(results_queues.datetime))*minutes_in_step,length(arcid))
    real_queue_use = zeros(Float64,length(unique(results_queues.datetime))*minutes_in_step,length(nodes))

    real_od_queue = create_od_queue(nodeid,daterange,minutes_in_step)

    ## determine the final queue at each station
    real_station_queue = DataFrame(transpose(sum(real_od_queue,dims=2)[:,1,:]),nodes)
    r_s_q = DataFrame(transpose(sum(real_od_queue,dims=2)[:,1,:]),nodes)

    ## prepare the optimized entry per datetime
    optimized_entry = unstack(select(results_queues,[:datetime,:station,:allowed]),:station,:allowed)

    ## preallocate the state of each queue
    queue_state = fill(1,length(nodes))

    ## start the flow through the network
    for origin in eachindex(nodes)

        for minute in 1:size(real_arc_use,1)-1
            ## determine the number of people allowed to enter in the minute at station
            if r_s_q[queue_state[origin],origin] > 0

                if kind == "bound"
                    moved_minute = min(
                            optimized_entry[div(minute,minutes_in_step)+1,origin+1],
                            sum(real_station_queue[queue_state[origin]:div(minute,minutes_in_step)+1,origin])
                        )
                elseif kind == "unbound"
                    moved_minute = min(r_s_q[queue_state[origin],origin],max_enter)
                end

                
                if !(minimum(closed) <= Hour(optimized_entry[div(minute,minutes_in_step)+1,1]) <= maximum(closed))
                    ## dispatch the ratio according to each destination
                    for destination in eachindex(nodes)
                        previous = previous_nodes[origin,destination]
                        current_destination = destination
                        destination_time = minute + distance[origin,destination] + 1
                        if real_od_queue[origin,destination,queue_state[origin]] > 0
                            moving_to_destination = (real_od_queue[origin,destination,queue_state[origin]]/r_s_q[queue_state[origin],origin]) * moved_minute
                            while previous != 0
                                start_time = floor(Int64,destination_time - grapharcs.traveltime[arcid[(previous,current_destination)]])
                                if start_time <= size(real_arc_use,1)
                                    real_arc_use[start_time,arcid[(previous,current_destination)]] += moving_to_destination
                                end
                                current_destination = previous
                                previous = previous_nodes[origin,current_destination]
                                destination_time = start_time
                            end
                        end
                    end
                end

                ## remove the number of moved people from queue
                real_station_queue[queue_state[origin],origin] -= moved_minute
            end

            ## change to the next qeue, if queue of timestep is empty
            if real_station_queue[queue_state[origin],origin] <= 0 && queue_state[origin] < div(minute,minutes_in_step)+1
                queue_state[origin] += 1
                real_station_queue[queue_state[origin],origin] += real_station_queue[queue_state[origin]-1,origin]
                real_station_queue[queue_state[origin]-1,origin] = 0
            end
            if queue_state[origin] > div(minute,minutes_in_step)+1
                error("Queue state out of bounds for $origin")
            end
            real_queue_use[minute,origin] = sum(real_station_queue[queue_state[origin]:div(minute,minutes_in_step)+1,origin])
        end
    end

    for arc in axes(grapharcs,1)
        for minute in axes(real_arc_use,1)
            if real_arc_use[minute,arc]>0
                real_arc_use[minute,arc] = real_arc_use[minute,arc] / grapharcs.capacity[arc]
            end
        end
    end

    detailed_time = minimum(results_queues.datetime):Minute(1):maximum(results_queues.datetime)+Minute(minutes_in_step-1)
    anim = @animate for dt in 1:round(Int64,minutes_in_step/5):size(real_arc_use,1)-1
        begin
            new_plot = bar(
                vec(sum(real_arc_use[dt:(dt+round(Int64,minutes_in_step/5)-1),:],dims=1)./round(Int64,minutes_in_step/5)),
                group = grapharcs.category,
                size = (1920, 1080),
                title="Real Metroarc Utilization at timestep $(detailed_time[dt])",
                ylims=(0,1.5),
                ylabel = "Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=15mm,
            )
            hline!([safety],label="Restriction",linewidth=4,colour=:grey60,)
            hline!([1.0],label="Full Capacity",linewidth=4,colour=:firebrick,)
        end
    end
    gif(anim, "visuals/utilization_real_fps20_$(horizon)_$kind.gif", fps = 20)

    if kind == "unbound"
        anim = @animate for dt in 1:round(Int64,minutes_in_step/5):size(real_arc_use,1)-1
            begin
                new_plot = bar(
                    vec(sum(real_arc_use[dt:(dt+round(Int64,minutes_in_step/5)-1),:],dims=1)./round(Int64,minutes_in_step/5)),
                    group = grapharcs.category,
                    size = (1920, 1080),
                    title="Real Metroarc Utilization at timestep $(detailed_time[dt])",
                    ylabel = "Arc Utilization",
                    ylims=(0,maximum(real_arc_use)*1.2),
                    colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                    legend=:bottomright,
                    margin=15mm,
                )
                hline!([safety],label="Restriction",linewidth=4,colour=:grey60,)
                hline!([1.0],label="Full Capacity",linewidth=4,colour=:firebrick,)
            end
        end
        gif(anim, "visuals/utilization_real_fps20_$(horizon)_$(kind)_nolim.gif", fps = 20)
    end

    plot(real_arc_use[:,3])
end

simulate_metro(results_queues,nodeid,minutes_in_step,daterange,grapharcs,"bound")

simulate_metro(results_queues,nodeid,minutes_in_step,daterange,grapharcs,"unbound")

function compare_horizons()
    b = CSV.read("results/qeues_1",DataFrame)
    b.bench .= "Forecast Periods: 1"
    for x in 2:3
        b_new = CSV.read("results/qeues_$x",DataFrame)
        b_new.bench .= "Forecast Periods: $x"
        b = vcat(b,b_new)
    end

    b = groupby(b,[:bench,:datetime])
    b = combine(b, :allowed => sum, :queued => sum; renamecols = false)

    plot(b.datetime,b.queued,groups=b.bench,size = (1920, 1080),lw=2,title="Comparison of Periods")
    savefig("visuals/periods.pdf")
end

compare_horizons()
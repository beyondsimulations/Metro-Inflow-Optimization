# load packages
using Pkg
Pkg.activate("metroflow")

using CSV
using DataFrames
using DataStructures
using JuMP
using HiGHS
using Plots
using StatsPlots
using Dates
using Measures
using Graphs

# parameters for the actual model
safety = 0.90           # safety factor that limits the arc capacity
minutes_in_period = 30  # minutes in each period (in 15 minute intervals!)
max_enter = 200         # number of maximal entries per minute per station
scaling = 1.0           # scaling of the metro queue (to test lower or higher demand)
past_minutes = 60       # timeframe to consider from the past during the optimization

# define the start- and end time of the observed time horizon
# Make sure that the horizon contains only one shift!
start_time = DateTime("2022-11-27T05:00:00.00")
end_time = DateTime("2022-11-28T02:59:00.00")

# datetime preparation
daterange = Date(start_time):Date(end_time)
periodrange = start_time:Minute(minutes_in_period):end_time
nr_minutes = length(start_time:Minute(1):end_time)+1
nr_periods = ceil(Int64,nr_minutes/minutes_in_period)
@assert rem(minutes_in_period,15) == 0 "The length of each period has to be in 15 min intervalls."
@assert nr_minutes <= 1440 "The timeframe exceeds one day and does represent more than one shift."

# load the data
println("Loading data.")
grapharcs = CSV.read("data_metro/metroarcs.csv", DataFrame)
sort!(grapharcs, :category)

function load_demand()
    demand =[]
    for day in eachindex(daterange)
        if day == 1
            demand = CSV.read("data_demand/OD_$(daterange[day])v3.csv", DataFrame)
        else
            demand = vcat(demand, CSV.read("data_demand/OD_$(daterange[day])v3.csv", DataFrame))
        end
    end
    filter!(row -> row.datetime >= start_time && row.datetime <= end_time, demand)
    return demand::DataFrame
end

function aggregate_demand()
    demand = load_demand()
    period_mapping = start_time:Minute(minutes_in_period):end_time
    demand.period .= 1
    for period in eachindex(period_mapping)
        for row in eachrow(demand)
            if period == 1 && row.datetime <= period_mapping[period]
                row.period = period
            elseif period > 1 && row.datetime > period_mapping[period-1] && row.datetime <= period_mapping[period]
                row.period = period
            end
        end
    end
    demand = groupby(demand,[:origin,:destination,:period])
    demand = combine(demand, :value => sum => :value)
    return demand
end

println("Preparing demand data.")
demand = aggregate_demand()

struct Metroarc
    origin::String
    destination::String
    category::String
    capacity::Int64
    traveltime::Int64
end

metroarcs::Vector{Metroarc} = []
for arc in axes(grapharcs,1)
    push!(metroarcs, Metroarc(
        grapharcs.origin[arc],
        grapharcs.destination[arc],
        grapharcs.category[arc],
        grapharcs.capacity[arc],
        ceil(Int64,grapharcs.traveltime[arc])
        )
    )
end

# prepare the graph related hash tables
println("Preparing graph data.")
nodes = sort(unique(vcat(getproperty.(metroarcs, :origin),getproperty.(metroarcs, :destination))))
nr_nodes = length(nodes)
nr_arcs = size(grapharcs,1)
d_node_id = Dict([nodes[i] => i for i in eachindex(nodes)])
d_arc_id = Dict((d_node_id[grapharcs.origin[i]],d_node_id[grapharcs.destination[i]]) => i for i in axes(grapharcs,1))
d_id_arc = Dict(i => (d_node_id[grapharcs.origin[i]],d_node_id[grapharcs.destination[i]]) for i in axes(grapharcs,1))

# prepare the graph of the metro network
G = DiGraph(nr_nodes)
distance_matrix = zeros(Int64,nr_nodes,nr_nodes)
for arc in metroarcs
    add_edge!(G,d_node_id[arc.origin],d_node_id[arc.destination])
    distance_matrix[d_node_id[arc.origin],d_node_id[arc.destination]] = ceil(Int64,arc.traveltime)
end

println("Preparing shifting data.")
shift_original = [Tuple{Int64,Int64,Int64}[] for _ in 1:nr_arcs, _ in 1:nr_minutes]
shift_start_end = [Tuple{Int64,Int64}[] for _ in 1:nr_nodes, _ in 1:nr_nodes]
for origin in 1:nr_nodes # origin at which people are allowed into the metro
    state = dijkstra_shortest_paths(G, origin, distance_matrix; trackvertices=true)
    all_distances = state.dists
    all_paths = enumerate_paths(state) # holds all nodes on the path from origin to destination
    for path in all_paths
        if path != []
            stops_in_path = length(path)
            minute_at_arc_start = 1
            previous_stop = origin # initialize the first stop as origin of the path
            for next_node in path # we now iterate over all nodes on a path
                if next_node != origin # we skip the first node in the path, as we need the arc
                    push!(shift_start_end[origin,path[end]],(d_arc_id[(previous_stop,next_node)],minute_at_arc_start))
                    minute_at_arc_start += metroarcs[d_arc_id[(previous_stop,next_node)]].traveltime # add the traveltime to the arc
                end
                previous_stop = next_node # set the previous node to the current node
            end
        end
    end
    for period in 1:nr_periods # period in which people are allowed into the metro
        for destination in 1:nr_nodes # destination that the people allowed at the origin want to reach
            if origin != destination # skips cases where no travel is necessary
                stops_in_path = length(all_paths[destination]) # saves the number of stops
                previous_stop = origin # initialize the first stop as origin of the path
                minute_at_arc_start = (period-1) * minutes_in_period + 1 # first minute of the current period
                start_at_origin = minute_at_arc_start # used for the assertation of the path
                for next_node in all_paths[destination] # we now iterate over all nodes on a path
                    if next_node != origin # we skip the first node in the path, as we need the arc
                        for minute in minute_at_arc_start:(minute_at_arc_start+minutes_in_period-1) # access all minutes of that period
                            if minute <= nr_minutes # make sure that we don't exeed the overall timeframe
                                #shift[minute,d_arc_id[(previous_stop,next_node)],period,origin,destination] = true # allocate all corresponding periods
                                push!(shift_original[d_arc_id[(previous_stop,next_node)],minute],(origin,destination,period))
                            end
                        end
                        minute_at_arc_start += metroarcs[d_arc_id[(previous_stop,next_node)]].traveltime # add the traveltime to the arc
                        previous_stop = next_node # set the previous node to the current node
                        if next_node == destination
                            expected_distance = start_at_origin + all_distances[destination]
                            @assert expected_distance == minute_at_arc_start "Shortest Paths: Duration computed $minute_at_arc_start, expected $expected_distance"
                        end
                    end
                end
            end
        end
    end
end

shift = copy(shift_original)
for a in 1:nr_arcs
    for t in 1:nr_minutes
        current_element = shift[a,t]
        for tt in t+1:nr_minutes
            if shift[a,tt] == shift[a,t]
                shift[a,t] = []
            end
        end
    end
end

# prepare the demand data
demand_od = zeros(Float64,nr_nodes,nr_nodes,nr_periods)
for row in eachrow(demand)
    demand_od[d_node_id[row.origin],d_node_id[row.destination],row.period] = row.value
end
@assert sum(demand_od) == sum(demand.value) "Mismatch in demand matrix."

struct MetroInstance
    nr_nodes::Int64                       # Overall number of nodes
    nr_arcs::Int64                        # Overall number of arcs
    nr_periods::Int64                     # Overall number of periods
    nr_minutes::Int64                     # Overall number of minutes
    capacity_arcs::Vector{Int64}          # Capacity of each arc
    safety_factor::Float64                # Safety factor for the arc_capacity
    max_entry_origin::Int64               # Maximal number of people allowed to enter from outside
    cum_demand_od_in_period::Array{Float64,3} # Array with the cummulated demand from o to d for each period
    demand_od_in_period::Array{Float64,3} # Array with the demand from o to d for each period
    shift::Matrix{Vector{Tuple{Int64, Int64, Int64}}}
end

function build_model_instance(cum_demand_od,demand_od)
    println("Preparing model instance.")
    modelInstance = MetroInstance(
        nr_nodes,
        nr_arcs,
        nr_periods,
        nr_minutes,
        getproperty.(metroarcs, :capacity),
        safety,
        max_enter,
        cum_demand_od,
        demand_od,
        shift
    )
    return modelInstance
end

function build_optimization_model(modelInstance)
    im = Model(HiGHS.Optimizer)
    set_attribute(im, "presolve", "on")
    set_attribute(im, "time_limit", 120.0)
    set_attribute(im, "mip_rel_gap", 0.0)
    println("Preparing optimization model.")
    @variable(im, 
        0 .<= X[o=1:modelInstance.nr_nodes,p=1:modelInstance.nr_periods] .<= modelInstance.max_entry_origin * safety
    )
    println("Preparing objective function.")
    @objective(im, Min, 
        sum((sum(modelInstance.cum_demand_od_in_period[o,d,p] for d in 1:modelInstance.nr_nodes) - X[o,p] *  minutes_in_period)^2 for o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods)
    )
    println("Preparing capacity constraints.")
    @constraint(im, capacity[a in 1:modelInstance.nr_arcs,t in 1:modelInstance.nr_minutes, p_shifts in 0:ceil(Int,past_minutes/minutes_in_period); shift[a,t] != []],
        sum(X[o,p] * modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
    )
    return im,X
end

function build_restricted_optimization_model(modelInstance,past_minutes,minutes_in_period,current_period)
    im = Model(HiGHS.Optimizer)
    set_attribute(im, "presolve", "on")
    set_attribute(im, "time_limit", 120.0)
    set_attribute(im, "mip_rel_gap", 0.0)
    println("Preparing optimization model.")
    @variable(im, 
        0 .<= X[o=1:modelInstance.nr_nodes,p=1:modelInstance.nr_periods] .<= modelInstance.max_entry_origin * safety
    )
    lower_period::Int64 = max(1,current_period - floor(past_minutes/minutes_in_period))
    upper_period::Int64 = min(current_period + ceil(past_minutes/minutes_in_period),modelInstance.nr_periods)
    current_minute::Int64 = ceil(minutes_in_period * current_period)

    println("Preparing objective function.")
    @objective(im, Min, 
        sum((sum(modelInstance.cum_demand_od_in_period[o,d,p] for d in 1:modelInstance.nr_nodes) - X[o,p] *  minutes_in_period)^2 for o in 1:modelInstance.nr_nodes, p in lower_period:upper_period)
    )
    println("Preparing capacity constraints.")
    @constraint(im, capacity[a in 1:modelInstance.nr_arcs,t in max(1,current_minute-past_minutes):min(modelInstance.nr_minutes,current_minute+past_minutes), p_shifts in 0:ceil(Int,past_minutes/minutes_in_period); shift[a,t] != []],
        sum(X[o,p] * modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
    )
    return im,X
end

function heuristic_adding_queues()
    demand_od_heuristic = copy(demand_od)
    remaining_queue = copy(demand_od) .= 0
    inflow_raw = zeros(Float64,nr_nodes,nr_periods)

    optimization_duration = zeros(Float64,nr_periods)

    for fix_period in 1:nr_periods
        println("Running period ", fix_period)
        instance = build_model_instance(demand_od_heuristic,demand_od)
        #model,X = build_optimization_model(instance)
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

            demand_od_heuristic[:,:,fix_period+1] .+= remaining_queue[:,:,fix_period]
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

    if isfile("logfile.csv")
        log = CSV.read("logfile.csv", DataFrame)
    else
        log = DataFrame(
            timestamp=DateTime[],
            intervall=String[],
            period_length=Int64[],
            past_minutes=Int64[],
            avg_duration=Float64[],
            avg_utilization=Float64[],
            max_utilization=Float64[],
            exceeded_minutes=Int64[],
            queue_at_end=Float64[],
            avg_wait_time=Float64[]
            )
    end

    return inflow_raw, results_queues, results_arcs
end

# Start of the iterative allocation
optimized_inflow, queues, arcs = heuristic_adding_queues()

function plot_optimization!(results_arcs)
    Plots.scalefontsizes()
    Plots.scalefontsizes(2.2)
    timesteps = unique(results_arcs.datetime)
    gp_results_arcs = groupby(results_arcs,:datetime)
    get_group(gdf, keys...) = gdf[(keys...,)]

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_arcs, dt).utilization_aggregated,
                group = get_group(gp_results_arcs, dt).line,
                size = (1920, 1080),
                title="Real Metroarc Usage at timestep $dt",
                ylims=(0,1.5),
                ylabel = "Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=25mm,
            )
            hline!([safety],label="Restriction",linewidth=2,colour=:grey60,)
        end
    end
    gif(anim, "visuals/opt_utilization_fps30_aggregated_mip-$(minutes_in_period)_pam-$(past_minutes).gif", fps = 30)

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_arcs, dt).utilization_period,
                group = get_group(gp_results_arcs, dt).line,
                size = (1920, 1080),
                title="Real Metroarc Usage at timestep $dt",
                ylims=(0,1.5),
                ylabel = "Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=25mm,
            )
            hline!([safety],label="Restriction",linewidth=2,colour=:grey60,)
        end
    end
    gif(anim, "visuals/opt_utilization_fps30_period_mip-$(minutes_in_period)_pam-$(past_minutes).gif", fps = 30)
end

plot_optimization!(arcs)

CSV.write("results/queues_new_mip-$(minutes_in_period)_pam-$(past_minutes).csv", queues)
CSV.write("results/arcs_new_mip-$(minutes_in_period)_pam-$(past_minutes).csv", arcs)

################

function create_od_queue!(real_od_queue)
    demand = load_demand()
    timesteps = start_time:Minute(1):end_time
    timestep_id = Dict(timesteps[i] => i for i in eachindex(timesteps))
    for data in eachrow(demand)
        for minutestep in data.datetime:Minute(1):data.datetime+Minute(14)
            real_od_queue[timestep_id[minutestep],d_node_id[data.origin],d_node_id[data.destination]] = (data.value * scaling)/15
        end
    end

end

function create_entry_list!(real_allowed_entry,queues)
    timesteps = start_time:Minute(1):end_time
    timestep_id = Dict(timesteps[i] => i for i in eachindex(timesteps))
    for entry in eachrow(queues)
        for minutestep in entry.datetime:Minute(1):entry.datetime+Minute(minutes_in_period)-Minute(1)
            real_allowed_entry[timestep_id[minutestep],d_node_id[entry.station]] = max(1,entry.allowed)
        end
    end
end

# Simulation
function simulate_metro(queues,nr_minutes,grapharcs,kind)

    real_arc_use = zeros(Float64,nr_minutes,nr_arcs)
    real_queue_use = zeros(Float64,nr_minutes,nr_nodes)
    real_od_queue = zeros(Float64,nr_minutes,nr_nodes,nr_nodes)

    ## load all queues into a tensor
    create_od_queue!(real_od_queue)

    ## create the list that saves the number of allowed people
    real_allowed_entry = zeros(Float64,nr_minutes,nr_nodes)
    create_entry_list!(real_allowed_entry,queues)


    ## start the flow through the network
    for minute in 1:size(real_arc_use,1)
        for origin in eachindex(nodes)

        ## determine the number of people allowed to enter in the minute at station
            if kind == "bound"
                moved_minute = min(real_allowed_entry[minute,origin],real_queue_use[max(minute-1,1),origin])

            elseif kind == "unbound"
                moved_minute = max(real_queue_use[max(minute-1,1),origin],max_enter)
            end

            ## dispatch the ratio according to each destination
            if moved_minute > 0
                for queue_minute in 1:minute
                    if any(real_od_queue[queue_minute,origin,:] .> 0)
                        queue_length = sum(real_od_queue[queue_minute,origin,:])
                        all_inflow = min(moved_minute,queue_length)
                        split_flow = (real_od_queue[queue_minute,origin,:]/queue_length) .* all_inflow
                        for destination in 1:nr_nodes
                            if split_flow[destination] > 0
                                for movement in shift_start_end[origin,destination]
                                    if minute+movement[2] <= size(real_arc_use,1)
                                        real_arc_use[minute+movement[2],movement[1]] += split_flow[destination]
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
        end
    end

    for arc in axes(grapharcs,1)
        for minute in axes(real_arc_use,1)
            if real_arc_use[minute,arc] > 0
                real_arc_use[minute,arc] = real_arc_use[minute,arc] / grapharcs.capacity[arc]
            end
        end
    end

    sim_queues = DataFrame(datetime = DateTime[],station=String[],queued=Int64[],allowed=Int64[])
    sim_arcs   = DataFrame(datetime = DateTime[],connection=Int64[],line=String[],utilization=Float64[])

    timesteps = start_time:Minute(1):end_time

    for t in 1:nr_minutes-1
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

    gp_results_queues = groupby(sim_queues,:datetime)
    gp_results_arcs = groupby(sim_arcs,:datetime)
    get_group(gdf, keys...) = gdf[(keys...,)]

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_queues, dt).queued,
                size = (1920, 1080),
                title="Queue at timestep $dt",
                ylims=(0,100000),
                label="People in Queue",
                legend=:topright,
                colour=:steelblue2,
                margin=25mm,
                xticks=(1:length(nodes),get_group(gp_results_queues, dt).station),
                xrotation=30,
                xtickfontsize=10,
            )
        end
    end
    gif(anim, "visuals/sim_queues_fp30_$(kind)_mip-$(minutes_in_period)_pam-$(past_minutes).gif", fps = 30)

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_queues, dt).allowed,
                label="People allowed to Enter per Minute",
                size = (1920, 1080),
                title="Entry per Minute at timestep $dt",
                ylims=(0,max_enter*1.1),
                legend=:topright,
                margin=25mm,
                xticks=(1:length(nodes),get_group(gp_results_queues, dt).station),
                xrotation=30,
                xtickfontsize=10,
                colour=:orange1,
            )
        end
    end
    gif(anim, "visuals/sim_entry_fps30_$(kind)_mip-$(minutes_in_period)_pam-$(past_minutes).gif", fps = 30)

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_arcs, dt).utilization,
                group = get_group(gp_results_arcs, dt).line,
                size = (1920, 1080),
                title="Real Metroarc Usage at timestep $dt",
                ylims=(0,1.5),
                ylabel = "Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=25mm,
            )
            hline!([safety],label="Restriction",linewidth=2,colour=:grey60,)
        end
    end
    gif(anim, "visuals/sim_utilization_fps30_$(kind)_mip-$(minutes_in_period)_pam-$(past_minutes).gif", fps = 30)

    return sim_queues,sim_arcs
end

sim_queues,sim_arcs = simulate_metro(queues,nr_minutes,grapharcs,"bound")
CSV.write("results/sim_queues_new_bound_mip-$(minutes_in_period)_pam-$(past_minutes).csv", sim_queues)

sim_queues,sim_arcs = simulate_metro(queues,nr_minutes,grapharcs,"unbound")
CSV.write("results/sim_arcs_new_unbound_mip-$(minutes_in_period)_pam-$(past_minutes).csv", sim_arcs)
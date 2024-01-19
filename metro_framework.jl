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
using SparseArrays

# parameters for the actual model
safety = 0.7            # safety factor that limits the arc capacity
minutes_in_period = 60  # minutes in each period (in 15 minute intervals!)
max_enter = 100         # number of maximal entries per minute per station
rate_closed = 100       # taxi "queue removal rate" per minute during closed metro hours
closed = Hour.(3:5)     # hours during which the metro operation is closed
scaling = 1.0           # scaling of the metro queue (to test lower or higher demand)

# datetime preparation
start_time = DateTime("2022-11-27T05:00:00.00")
end_time = DateTime("2022-11-28T04:59:00.00")
daterange = Date(start_time):Date(end_time)
nr_minutes = length(start_time:Minute(1):end_time)
nr_periods = ceil(Int64,nr_minutes/minutes_in_period)
@assert rem(minutes_in_period,15) == 0 "The length of each period has to be in 15min intervalls."

# load the data
println("Loading data.")
grapharcs = CSV.read("data_metro/metroarcs.csv", DataFrame)
sort!(grapharcs, :category)

function aggregate_demand()
    demand =[]
    for day in eachindex(daterange)
        if day == 1
            demand = CSV.read("data_demand/OD_$(daterange[day])v3.csv", DataFrame)
        else
            demand = vcat(demand, CSV.read("data_demand/OD_$(daterange[day])v3.csv", DataFrame))
        end
    end
    filter!(row -> row.datetime >= start_time && row.datetime <= end_time, demand)
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

println("Preparing graph data.")
# prepare the graph related hash tables
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
shift = [Tuple{Int64,Int64,Int64}[] for _ in 1:nr_minutes, _ in 1:nr_arcs]
for origin in 1:nr_nodes # origin at which people are allowed into the metro
    state = dijkstra_shortest_paths(G, origin, distance_matrix; trackvertices=true)
    all_distances = state.dists
    all_paths = enumerate_paths(state) # holds all nodes on the path from origin to destination
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
                                push!(shift[minute,d_arc_id[(previous_stop,next_node)]],(origin,destination,period))
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

# prepare the demand data
demand_od = zeros(Int64,nr_nodes,nr_nodes,nr_periods)
for row in eachrow(demand)
    demand_od[d_node_id[row.origin],d_node_id[row.destination],row.period] = row.value
end
@assert sum(demand_od) == sum(demand.value) "Mismatch in demand matrix."

struct MetroInstance
    nr_nodes::Int64                     # Overall number of nodes
    nr_arcs::Int64                      # Overall number of arcs
    nr_periods::Int64                   # Overall number of periods
    nr_minutes::Int64                   # Overall number of minutes
    capacity_arcs::Vector{Int64}        # Capacity of each arc
    safety_factor::Float64              # Safety factor for the arc_capacity
    max_entry_origin::Int64             # Maximal number of people allowed to enter from outside
    demand_od_in_period::Array{Int64,3} # Array with the demand from o to d for each period
    shift::Matrix{Vector{Tuple{Int64, Int64, Int64}}}
end

println("Preparing model instance.")
modelInstance = MetroInstance(
    nr_nodes,
    nr_arcs,
    nr_periods,
    nr_minutes,
    getproperty.(metroarcs, :capacity),
    safety,
    max_enter,
    demand_od,
    shift
)

im = Model(HiGHS.Optimizer)
set_attribute(im, "presolve", "on")
set_attribute(im, "time_limit", 120.0)
set_attribute(im, "mip_rel_gap", 0.0)

println("Preparing optimization model.")
@variable(im, 
    X[o=1:modelInstance.nr_nodes,p=1:modelInstance.nr_periods] .>= 0
)

println("Preparing objective function.")
@objective(im, Min, 
    sum(sum(sum(modelInstance.demand_od_in_period[o,d,p] for d in 1:modelInstance.nr_nodes) - X[o,p] *  minutes_in_period)^2 for o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods)
)

println("Preparing capacity constraints.")
for t in 1:minutes_in_period:modelInstance.nr_minutes 
    for a in 1:modelInstance.nr_arcs
        @constraint(im,
            sum(X[o,p] * modelInstance.demand_od_in_period[o,d,p]/sum(modelInstance.demand_od_in_period[o,:,p]) for (o,d,p) in modelInstance.shift[t,a] if modelInstance.demand_od_in_period[o,d,p] > 0) <=  modelInstance.capacity_arcs[a] * modelInstance.safety_factor
        )
    end
end

println("Preparing capacity inflow constraints.")
@constraint(
    im, capacity_station[o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods],
    X[o,p] <= modelInstance.max_entry_origin
)

println("Starting optimization.")
optimize!(im)

rem_queue_steps = zeros(Float64, nr_nodes, nr_periods)
queue_steps = sum(modelInstance.demand_od_in_period[:,:,:], dims=2)[:,1,:]
add_queue_steps = copy(queue_steps)
inflow = value.(X) .* minutes_in_period
for p in 1:nr_periods
    rem_queue_steps[:,p] = queue_steps[:,p] .- inflow[:,p]
    if p > 1
        rem_queue_steps[:,p] +=  rem_queue_steps[:,p-1]
        add_queue_steps[:,p] +=  add_queue_steps[:,p-1]
    end
end

#=
println("Preparing optimization model.")
@variable(im, 
    X[o=1:modelInstance.nr_nodes,p=1:modelInstance.nr_periods,u=max(1,p-max_shift):p] .>= 0
)

println("Preparing objective function.")
@objective(im, Min, 
    sum((sum(modelInstance.demand_od_in_period[o,d,u] for d in 1:modelInstance.nr_nodes, u in 1:p) - sum(X[o,pp,u] for pp in 1:p, u in 1:p))^2 for o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods)
)

println("Preparing capacity constraints.")
for t in 1:minutes_in_period:modelInstance.nr_minutes 
    for a in 1:modelInstance.nr_arcs
        @constraint(im,
            sum(X[o,p,u] * modelInstance.demand_od_in_period[o,d,u]/sum(modelInstance.demand_od_in_period[o,dd,u] for dd in 1:modelInstance.nr_nodes) for (o,d,u) in modelInstance.shift[t,a], p in 1:modelInstance.nr_periods if modelInstance.demand_od_in_period[o,d,u] > 0 && u <= p && u >= p - max_shift) <=  modelInstance.capacity_arcs[a] * modelInstance.safety_factor
        )
    end
end

println("Preparing capacity inflow constraints.")
@constraint(
    im, capacity_station[o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods],
    sum(X[o,p,u] for u in 1:modelInstance.nr_periods) <= modelInstance.max_entry_origin
)
println("Starting optimization.")
optimize!(im)

=#
which_node = 1
display(plot(inflow[which_node,:], label = "Inflow"))
display(plot!(queue_steps[which_node,:], label = "New Queue"))
display(plot!(rem_queue_steps[which_node,:], label = "Remaining Queue"))
display(plot!(add_queue_steps[which_node,:], label = "Cummulated Queue"))


display(plot(sum(inflow,dims=1)[:], label = "Inflow"))
display(plot!(sum(queue_steps,dims=1)[:], label = "New Queue"))
display(plot!(sum(rem_queue_steps,dims=1)[:], label = "Remaining Queue"))
display(plot!(sum(add_queue_steps,dims=1)[:], label = "Cummulated Queue"))

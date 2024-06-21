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
using Statistics

include("metro_functions.jl")
include("metro_model.jl")
include("metro_heuristic.jl")
include("metro_simulation.jl")
include("metro_visuals.jl")

# Parameters for the actual model
safety = 0.90           # safety factor that limits the arc capacity
minutes_in_period = 15  # minutes in each period (in 15 minute intervals!)
max_enter = 180         # number of maximal entries per minute per station
scaling = 1.0           # scaling of the metro queue (to test lower or higher demand)
past_minutes = 300      # timeframe to consider from the past during the optimization

# Define the start- and end time of the observed time horizon
# Make sure that the horizon contains only one shift!
start_time = DateTime("2022-11-28T05:00:00.00")
end_time = DateTime("2022-11-29T02:59:00.00")

struct MetroInstance
    nr_nodes::Int64                                     # Overall number of nodes
    nr_arcs::Int64                                      # Overall number of arcs
    nr_periods::Int64                                   # Overall number of periods
    nr_minutes::Int64                                   # Overall number of minutes
    capacity_arcs::Vector{Int64}                        # Capacity of each arc
    safety_factor::Float64                              # Safety factor for the arc_capacity
    max_entry_origin::Int64                             # Maximal number of people allowed to enter from outside
    cum_demand_od_in_period::Array{Float64,3}           # Array with the cummulated demand from o to d for each period
    demand_od_in_period::Array{Float64,3}               # Array with the demand from o to d for each period
    shift::Matrix{Vector{Tuple{Int64, Int64, Int64}}}   # Array with the timelag of entry and utilization
end

# Input validation and datetime preparation
daterange = Date(start_time):Date(end_time)
periodrange = start_time:Minute(minutes_in_period):end_time
nr_minutes = length(start_time:Minute(1):end_time)+1
nr_periods = ceil(Int64,nr_minutes/minutes_in_period)
@assert rem(minutes_in_period,15) == 0 "The length of each period has to be in 15 min intervalls."
@assert nr_minutes <= 1440 "The timeframe exceeds one day and does represent more than one shift."

# Load data of the metroarcs
println("Loading data.")
grapharcs = CSV.read("data_metro/metroarcs.csv", DataFrame)
sort!(grapharcs, :category)

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

# Prepare the graph related hash tables
println("Preparing graph data.")
nodes = sort(unique(vcat(getproperty.(metroarcs, :origin),getproperty.(metroarcs, :destination))))
nr_nodes = length(nodes)
nr_arcs = size(grapharcs,1)
d_node_id = Dict([nodes[i] => i for i in eachindex(nodes)])
d_arc_id = Dict((d_node_id[grapharcs.origin[i]],d_node_id[grapharcs.destination[i]]) => i for i in axes(grapharcs,1))
d_id_arc = Dict(i => (d_node_id[grapharcs.origin[i]],d_node_id[grapharcs.destination[i]]) for i in axes(grapharcs,1))

# Load and prepare the demand data
println("Preparing demand data.")
demand = aggregate_demand()
demand_od = zeros(Float64,nr_nodes,nr_nodes,nr_periods)
for row in eachrow(demand)
    demand_od[d_node_id[row.origin],d_node_id[row.destination],row.period] = row.value
end
@assert sum(demand_od) == sum(demand.value) "Mismatch in demand matrix."

# Prepare the time-shift in the metro data
println("Preparing shifting data.")
shift, shift_original, shift_start_end = compute_shift()

# Start of the iterative allocation
queues, arcs, opt_duration = heuristic_adding_queues()
plot_optimization!(arcs)

# Save the results from the heuristic
CSV.write("results/queues_new_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).csv", queues)
CSV.write("results/arcs_new_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).csv", arcs)

# Start the simulation unbound
println("Start simulation bound.")
sim_queues,sim_arcs = simulate_metro(queues,nr_minutes,grapharcs,"bound")
CSV.write("results/sim_queues_new_bound_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).csv", sim_queues)
plot_simulation!(sim_queues,sim_arcs,"bound")

# Start the simulation in unbound mode
#println("Start simulation inflowbound.")
#sim_queues,sim_arcs = simulate_metro(queues,nr_minutes,grapharcs,"inflow")
#CSV.write("results/sim_arcs_new_inflowbound_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).csv", sim_arcs)
#plot_simulation!(sim_queues,sim_arcs,"inflow")

# Start the simulation in unbound mode
#println("Start simulation unbound.")
#sim_queues,sim_arcs = simulate_metro(queues,nr_minutes,grapharcs,"unbound")
#CSV.write("results/sim_arcs_new_unbound_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).csv", sim_arcs)
#plot_simulation!(sim_queues,sim_arcs,"unbound")
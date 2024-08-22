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
using Gurobi
using ProgressMeter

include("metro_functions.jl")
include("metro_model.jl")
include("metro_heuristic.jl")
include("metro_simulation.jl")
include("metro_visuals.jl")

# Parameters for the actual model
set_safety = [0.9]                      # safety factor that limits the arc capacity
set_max_enter = [120]                           # number of maximal entries per minute per station
set_min_enter = [10]                              # min number of people allowed to enter
set_scaling = [1.0]                                 # scaling of the metro queue (to test lower or higher demand)
set_past_periods = [4]                      # timeframe to consider from the past during the optimization
set_kind_opt = ["regularSqr"]                        # "regularSqr","linweight"
set_kind_queue = ["lag_periods"]                    # "shift_periods","lag_periods"

# Define static simulation data
kind_sim = "bound"                  # "bound","inflow","unbound"
minutes_in_period = 30              # minutes in each period (in 15 minute intervals!)

# Define the start- and end time of the observed time horizon
# Make sure that the horizon contains only one shift!
start_time = DateTime("2022-11-29T05:00:00.00")
end_time = DateTime("2022-11-30T04:59:00.00")

struct MetroInstance
    kind_opt::String
    kind_queue::String
    nr_nodes::Int64                                     # Overall number of nodes
    nr_arcs::Int64                                      # Overall number of arcs
    nr_periods::Int64                                   # Overall number of periods
    nr_minutes::Int64                                   # Overall number of minutes
    past_minutes::Int64
    minutes_in_period::Int64
    capacity_arcs::Vector{Int64}                        # Capacity of each arc
    safety_factor::Float64                              # Safety factor for the arc_capacity
    min_entry_origin::Int64                             # Minimal number of people allowed to enter from outside
    max_entry_origin::Int64                             # Maximal number of people allowed to enter from outside
    cum_demand_od_in_period::Array{Float64,3}           # Array with the cummulated demand from o to d for each period
    demand_od_in_period::Array{Float64,3}               # Array with the demand from o to d for each period
    shift::Matrix{Vector{Tuple{Int64, Int64, Int64}}}   # Array with the timelag of entry and utilization
    scaling::Float64                                    # scaling of the metro queue (to test lower or higher demand)
    closed_period::Vector{Bool}                         # defines the periods with a closure of the metro
end

# Input validation and datetime preparation
set_past_minutes = set_past_periods .* minutes_in_period
daterange = Date(start_time):Date(end_time)
periodrange = start_time:Minute(minutes_in_period):end_time
nr_minutes = length(start_time:Minute(1):end_time)+ 120
nr_periods = ceil(Int64,(nr_minutes-120)/minutes_in_period)
closed_period = zeros(Bool,nr_periods)
for period in eachindex(periodrange)
    if hour(periodrange[period]) == 3 || hour(periodrange[period]) == 4 || hour(periodrange[period]) == 5
        closed_period[period] = true
    end
end
@assert rem(minutes_in_period,15) == 0 "The length of each period has to be in 15 min intervalls."

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
shift, shift_original, shift_start_end, longest_path = compute_shift()
println("Longest path is $longest_path.")

for min_enter in set_min_enter
    for safety in set_safety  
        for max_enter in set_max_enter
            for scaling in set_scaling
                for past_minutes in set_past_minutes
                    for kind_opt in set_kind_opt
                        for kind_queue in set_kind_queue
                            
                            println("Start: safety $safety, mip $minutes_in_period, pm $past_minutes, $kind_opt, $kind_sim, $kind_queue")

                            # Create a MetroInstance object with the following parameters:
                            # * kind_opt: 
                            # * kind_queue: 
                            # * nr_nodes: number of nodes in the network (i.e., metro stations)
                            # * nr_arcs: number of arcs between nodes in the network
                            # * nr_periods: number of time periods to consider in the problem
                            # * nr_minutes: number of minutes per period
                            # * past_minutes: timeframe of previous minutes considered during the optimization
                            # * getproperty.(metroarcs, :capacity): list of capacities for each arc 
                            # * safety: ratio of max. arc utilization
                            # * min_enter: minimum number of people that can enter a node at any minute
                            # * max_enter: maximum number of people that can enter a node at any minute
                            # * cum_demand_od: cummulated demand from origin to destination per period
                            # * demand_od: demand from origin to destination per period
                            # * shift: shift data that maps minutes and periods and arcs
                            # * scaling: test higher and lower demands by multiplication

                            modelInstance = MetroInstance(
                                kind_opt,
                                kind_queue,
                                nr_nodes,
                                nr_arcs,
                                nr_periods,
                                nr_minutes,
                                past_minutes,
                                minutes_in_period,
                                getproperty.(metroarcs, :capacity),
                                safety,
                                min_enter,
                                max_enter,
                                copy(demand_od),
                                copy(demand_od),
                                shift,
                                scaling,
                                closed_period,
                            )

                            # Start of the iterative allocation
                            queues, arcs, opt_duration, queue_period_age, infeasible_solutions = heuristic_adding_queues(modelInstance)
                            #plot_optimization!(arcs)

                            # Save the results from the heuristic
                            #CSV.write("results/queues_new_$(kind_opt)_$(kind_sim)_$(kind_queue)_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).csv", queues)
                            #CSV.write("results/arcs_new_$(kind_opt)_$(kind_sim)_$(kind_queue)_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).csv", arcs)

                            # Start the simulation
                            println("Start simulation $kind_sim.")
                            sim_queues,sim_arcs = simulate_metro(modelInstance,queues,opt_duration,grapharcs,kind_sim,queue_period_age,infeasible_solutions)

                            # Plot results
                            #plot_simulation!(sim_queues,sim_arcs,kind_sim)
                        end
                    end
                end
            end
        end
    end
end
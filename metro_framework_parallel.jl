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

include("functions/metro_functions.jl")
include("functions/metro_model.jl")
include("functions/metro_heuristic.jl")
include("functions/metro_simulation.jl")
include("functions/metro_visuals.jl")

# Parse command line arguments
function parse_args()
    if length(ARGS) != 3
        println("Usage: julia metro_framework_parallel.jl <minutes_in_period> <start_time> <end_time>")
        println("Example: julia metro_framework_parallel.jl 15 '2022-11-29T05:00:00.00' '2022-11-30T04:59:00.00'")
        exit(1)
    end

    minutes_in_period = parse(Int, ARGS[1])
    start_time = DateTime(ARGS[2])
    end_time = DateTime(ARGS[3])

    # Validate minutes_in_period
    if !(minutes_in_period in [15, 30, 45, 60])
        println("Error: minutes_in_period must be one of: 15, 30, 45, 60")
        exit(1)
    end

    # Validate that minutes_in_period is divisible by 15
    if rem(minutes_in_period, 15) != 0
        println("Error: The length of each period has to be in 15 min intervals.")
        exit(1)
    end

    return minutes_in_period, start_time, end_time
end

# Parse arguments
minutes_in_period, start_time, end_time = parse_args()

println("Running with parameters:")
println("  minutes_in_period: $minutes_in_period")
println("  start_time: $start_time")
println("  end_time: $end_time")

# Parameters for the actual model
set_safety = [0.6,0.7,0.8,0.9,1.0]        # safety factor that limits the arc capacity
set_max_enter = [80,160,240]              # number of maximal entries per minute per station
set_min_enter = [0,10]                    # min number of people allowed to enter
set_scaling = [0.8,1.0,1.2]               # scaling of the metro queue (to test lower or higher demand)
set_past_periods = [1,2,3,4,5,6]          # timeframe to consider from the past during the optimization
set_kind_opt = ["linweight"]              # "regularSqr","linweight"
set_kind_queue = ["lag_periods"]          # "shift_periods","lag_periods"

# Define static simulation data
kind_sim = "bound"                        # "bound","inflow","unbound"

struct MetroInstance
    kind_opt::String
    kind_queue::String
    nr_nodes::Int64                                     # Overall number of nodes
    nr_arcs::Int64                                      # Overall number of arcs
    nr_periods::Int64                                   # Overall number of periods
    nr_minutes::Int64                                   # Overall number of minutes
    past_minutes::Int64
    minutes_in_period::Int64                            # Number of minutes in each period
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

# Create results directories with parameter-specific subdirectories
result_suffix = "mip_$(minutes_in_period)_$(Dates.format(start_time, "yyyymmdd_HHMM"))_$(Dates.format(end_time, "yyyymmdd_HHMM"))"
if !isdir("results")
    mkdir("results")
end
if !isdir("results/queues")
    mkdir("results/queues")
end
if !isdir("results/arcs")
    mkdir("results/arcs")
end
if !isdir("results/visuals")
    mkdir("results/visuals")
end

# Create parameter-specific subdirectories
for subdir in ["queues", "arcs", "visuals"]
    param_dir = "results/$subdir/$result_suffix"
    if !isdir(param_dir)
        mkpath(param_dir)
    end
end

# Load data of the metroarcs
println("Loading data.")
grapharcs = CSV.read("data_metro/metroarcs.csv", DataFrame)
sort!(grapharcs, :category)

# Create Metroarc struct
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

# Create a log file for this specific run
log_file = open("results/log_$result_suffix.txt", "w")
println(log_file, "Metro Framework Run Log")
println(log_file, "Started at: $(now())")
println(log_file, "Parameters:")
println(log_file, "  minutes_in_period: $minutes_in_period")
println(log_file, "  start_time: $start_time")
println(log_file, "  end_time: $end_time")
println(log_file, "  nr_periods: $nr_periods")
println(log_file, "  nr_minutes: $nr_minutes")
println(log_file, "")
flush(log_file)

total_combinations = length(set_min_enter) * length(set_safety) * length(set_max_enter) *
                    length(set_scaling) * length(set_past_minutes) * length(set_kind_opt) *
                    length(set_kind_queue)

println("Total parameter combinations to process: $total_combinations")
println(log_file, "Total parameter combinations to process: $total_combinations")
flush(log_file)

global combination_count = 0

for min_enter in set_min_enter
    for safety in set_safety
        for max_enter in set_max_enter
            for scaling in set_scaling
                for past_minutes in set_past_minutes
                    for kind_opt in set_kind_opt
                        for kind_queue in set_kind_queue
                            global combination_count += 1

                            println("[$combination_count/$total_combinations] Start: safety $safety, mip $minutes_in_period, pm $past_minutes, $kind_opt, $kind_sim, $kind_queue")
                            println(log_file, "[$combination_count/$total_combinations] Processing: safety=$safety, max_enter=$max_enter, min_enter=$min_enter, scaling=$scaling, past_minutes=$past_minutes, kind_opt=$kind_opt, kind_queue=$kind_queue")
                            flush(log_file)

                            # Create a MetroInstance object with the following parameters:
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

                            try
                                # Start of the iterative allocation
                                queues, arcs, opt_duration, queue_period_age, infeasible_solutions = heuristic_adding_queues(modelInstance)

                                # Start the simulation
                                println("Start simulation $kind_sim.")
                                sim_queues,sim_arcs = simulate_metro(modelInstance,queues,opt_duration,grapharcs,kind_sim,queue_period_age,infeasible_solutions)

                                # Save results with parameter-specific naming
                                result_name = "safety_$(safety)_maxenter_$(max_enter)_minenter_$(min_enter)_scaling_$(scaling)_past_$(past_minutes)_$(kind_opt)_$(kind_queue)"

                                # You can add specific result saving here if needed
                                println(log_file, "  ✓ Completed successfully")

                                # Plot results (uncomment if needed)
                                # plot_optimization!(arcs)
                                # plot_simulation!(sim_queues,sim_arcs,kind_sim)

                            catch e
                                println("Error in combination $combination_count: $e")
                                println(log_file, "  ✗ Error: $e")
                            end

                            flush(log_file)
                        end
                    end
                end
            end
        end
    end
end

println("Completed all combinations for minutes_in_period=$minutes_in_period")
println(log_file, "")
println(log_file, "Completed at: $(now())")
close(log_file)

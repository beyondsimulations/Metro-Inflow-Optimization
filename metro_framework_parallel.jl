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
using TOML

include("functions/config.jl")
include("functions/metro_functions.jl")
include("functions/metro_model.jl")
include("functions/metro_heuristic.jl")
include("functions/metro_simulation.jl")
include("functions/metro_visuals.jl")

# Parse command line arguments
function parse_args()
    # Check for flags (--config and --unbound can appear in any order before positional args)
    config_path = "config/doha.toml"  # Default configuration
    unbound_mode = false
    positional_args = String[]

    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--config" && i + 1 <= length(ARGS)
            config_path = ARGS[i + 1]
            i += 2
        elseif ARGS[i] == "--unbound"
            unbound_mode = true
            i += 1
        else
            push!(positional_args, ARGS[i])
            i += 1
        end
    end

    # Check positional arguments
    if length(positional_args) != 3
        println("Usage: julia metro_framework_parallel.jl [--config <config_file>] [--unbound] <minutes_in_period> <start_time> <end_time>")
        println("")
        println("Options:")
        println("  --config <file>  Configuration file (default: config/doha.toml)")
        println("  --unbound        Run baseline simulation without optimization (only varies scaling)")
        println("")
        println("Examples:")
        println("  # Optimized run (default)")
        println("  julia metro_framework_parallel.jl --config config/shanghai.toml 60 '2017-05-15T05:00:00' '2017-05-16T04:59:00'")
        println("")
        println("  # Baseline run (no optimization, just simulation)")
        println("  julia metro_framework_parallel.jl --config config/shanghai.toml --unbound 60 '2017-05-15T05:00:00' '2017-05-16T04:59:00'")
        exit(1)
    end

    # Load configuration
    config = load_config(config_path)

    minutes_in_period = parse(Int, positional_args[1])
    start_time = DateTime(positional_args[2])
    end_time = DateTime(positional_args[3])

    # Validate minutes_in_period
    if minutes_in_period <= 0
        println("Error: minutes_in_period must be positive, got $minutes_in_period")
        exit(1)
    end
    base_interval = config.interval_minutes
    if rem(minutes_in_period, base_interval) != 0
        println("Error: minutes_in_period ($minutes_in_period) must be a multiple of the base interval ($base_interval) for $(config.display_name)")
        exit(1)
    end

    return config, minutes_in_period, start_time, end_time, unbound_mode
end

# Parse arguments
config, minutes_in_period, start_time, end_time, unbound_mode = parse_args()

println("Running with parameters:")
println("  region: $(config.display_name)")
println("  minutes_in_period: $minutes_in_period")
println("  start_time: $start_time")
println("  end_time: $end_time")
println("  mode: $(unbound_mode ? "UNBOUND (baseline)" : "BOUND (optimized)")")

# Parameters from config file (reproducible per region)
set_safety = config.safety_factors        # safety factor that limits the arc capacity
set_max_enter = config.max_enter          # number of maximal entries per minute per station
set_min_enter = config.min_enter          # min number of people allowed to enter
set_scaling = config.scaling_factors      # scaling of the metro queue (to test lower or higher demand)
set_past_periods = config.past_periods    # timeframe to consider from the past during the optimization
set_kind_opt = config.kind_opt            # "regularSqr","linweight"
set_kind_queue = config.kind_queue        # "shift_periods","lag_periods"

if unbound_mode
    println("\nUNBOUND MODE: Skipping optimization, only varying scaling factors")
    println("  scaling_factors: $set_scaling")
else
    println("\nOptimization parameters from config:")
    println("  safety_factors: $set_safety")
    println("  min_enter: $set_min_enter")
    println("  max_enter: $set_max_enter")
    println("  scaling_factors: $set_scaling")
    println("  past_periods: $set_past_periods")
    println("  kind_opt: $set_kind_opt")
    println("  kind_queue: $set_kind_queue")
end

# Define simulation mode based on --unbound flag
kind_sim = unbound_mode ? "unbound" : "bound"

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
    if hour(periodrange[period]) in config.closed_hours
        closed_period[period] = true
    end
end

# Create results directories with parameter-specific subdirectories
result_suffix = "$(config.name)_mip_$(minutes_in_period)_$(Dates.format(start_time, "yyyymmdd_HHMM"))_$(Dates.format(end_time, "yyyymmdd_HHMM"))"
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
mkpath("results/plots")

# Load data of the metroarcs
println("Loading data for $(config.display_name).")
grapharcs = CSV.read(get_metroarcs_path(config), DataFrame)
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
# Demand data preparation (prints progress inside aggregate_demand)
demand = aggregate_demand(config)
demand_od = zeros(Float64,nr_nodes,nr_nodes,nr_periods)
for row in eachrow(demand)
    demand_od[d_node_id[row.origin],d_node_id[row.destination],row.period] = row.value
end
@assert sum(demand_od) == sum(demand.value) "Mismatch in demand matrix."

# Prepare the time-shift in the metro data
shift, shift_original, shift_start_end, longest_path = compute_shift(demand_od)
println("Longest path is $longest_path.")

# Create a log file for this specific run
log_file = open("results/log_$result_suffix.txt", "w")
println(log_file, "Metro Framework Run Log")
println(log_file, "Started at: $(now())")
println(log_file, "Parameters:")
println(log_file, "  region: $(config.display_name)")
println(log_file, "  config_file: $(config.name)")
println(log_file, "  minutes_in_period: $minutes_in_period")
println(log_file, "  start_time: $start_time")
println(log_file, "  end_time: $end_time")
println(log_file, "  nr_periods: $nr_periods")
println(log_file, "  nr_minutes: $nr_minutes")
println(log_file, "")
flush(log_file)

if unbound_mode
    # UNBOUND MODE: Only vary scaling, skip optimization entirely
    total_combinations = length(set_scaling)
    println("Total scaling scenarios to process: $total_combinations")
    println(log_file, "Mode: UNBOUND (baseline simulation)")
    println(log_file, "Total scaling scenarios to process: $total_combinations")
    flush(log_file)

    global combination_count = 0
    progress = Progress(total_combinations, desc="Scaling scenarios: ", showspeed=true)

    for scaling in set_scaling
        global combination_count += 1

        println("/n[$combination_count/$total_combinations] UNBOUND: scaling=$scaling")
        println(log_file, "[$combination_count/$total_combinations] Processing: scaling=$scaling (unbound baseline)")
        flush(log_file)

        # Use default/dummy values for optimization parameters (not used in unbound)
        safety = 1.0
        min_enter = 0
        max_enter = set_max_enter[1]
        past_minutes = set_past_minutes[1]
        kind_opt = "none"
        kind_queue = "none"

        # Create a minimal MetroInstance (optimization params don't matter for unbound)
        # Apply scaling to demand here so optimization and simulation use same values
        scaled_demand = copy(demand_od) .* scaling
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
            scaled_demand,
            copy(scaled_demand),
            shift,
            scaling,
            closed_period,
        )

        try
            # Create dummy optimization outputs (not used in unbound simulation)
            # queues DataFrame: needs datetime, station, allowed columns
            periodrange_local = start_time:Minute(minutes_in_period):end_time
            queues = DataFrame(
                datetime = repeat(collect(periodrange_local), inner=nr_nodes),
                station = repeat(nodes, outer=length(periodrange_local)),
                allowed = zeros(Float64, nr_nodes * length(periodrange_local)),
                queued = zeros(Float64, nr_nodes * length(periodrange_local))
            )
            # Set allowed = max_enter for open periods, 0 for closed (respects metro hours)
            for row in eachrow(queues)
                period_hour = hour(row.datetime)
                if period_hour in config.closed_hours
                    row.allowed = 0
                else
                    row.allowed = max_enter  # Large value, but doesn't matter for unbound
                end
            end

            opt_duration = zeros(Float64, nr_periods)
            build_duration = zeros(Float64, nr_periods)
            queue_period_age = zeros(Float64, nr_periods, nr_nodes)
            infeasible_solutions = 0

            # Run simulation directly (no optimization)
            println("Start simulation $kind_sim (baseline - no optimization).")
            sim_queues, sim_arcs = simulate_metro(modelInstance, queues, opt_duration, build_duration, grapharcs, kind_sim, queue_period_age, infeasible_solutions, config)

            println(log_file, "  ✓ Completed successfully")

            # Generate visualization GIFs (uncomment to enable - adds ~15min per run)
            # println("Generating visualization GIFs...")
            # plot_simulation!(sim_queues, sim_arcs, kind_sim, config, kind_opt, kind_queue, minutes_in_period, past_minutes, start_time, safety, max_enter, nodes)

        catch e
            println("Error in unbound scenario $combination_count: $e")
            println(log_file, "  ✗ Error: $e")
            showerror(stdout, e, catch_backtrace())
        end

        flush(log_file)
        next!(progress)
    end

else
    # BOUND MODE: Full optimization with all parameter combinations (existing behavior)
    total_combinations = length(set_min_enter) * length(set_safety) * length(set_max_enter) *
                        length(set_scaling) * length(set_past_minutes) * length(set_kind_opt) *
                        length(set_kind_queue)

    println("Total parameter combinations to process: $total_combinations")
    println(log_file, "Mode: BOUND (optimized)")
    println(log_file, "Total parameter combinations to process: $total_combinations")
    flush(log_file)

    global combination_count = 0
    progress = Progress(total_combinations, desc="Parameter combinations: ", showspeed=true)

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
                                # Apply scaling to demand here so optimization and simulation use same values
                                scaled_demand = copy(demand_od) .* scaling
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
                                    scaled_demand,
                                    copy(scaled_demand),
                                    shift,
                                    scaling,
                                    closed_period,
                                )

                                try
                                    # Start of the iterative allocation
                                    queues, arcs, opt_duration, build_duration, queue_period_age, infeasible_solutions = heuristic_adding_queues(modelInstance, config)

                                    # Start the simulation
                                    println("Start simulation $kind_sim.")
                                    sim_queues,sim_arcs = simulate_metro(modelInstance,queues,opt_duration,build_duration,grapharcs,kind_sim,queue_period_age,infeasible_solutions,config)

                                    # Save results with parameter-specific naming
                                    result_name = "safety_$(safety)_maxenter_$(max_enter)_minenter_$(min_enter)_scaling_$(scaling)_past_$(past_minutes)_$(kind_opt)_$(kind_queue)"

                                    # You can add specific result saving here if needed
                                    println(log_file, "  ✓ Completed successfully")

                                    # Generate visualization GIFs (uncomment to enable - adds ~15min per run)
                                    # println("Generating visualization GIFs...")
                                    # plot_optimization!(arcs, config, kind_opt, kind_queue, minutes_in_period, past_minutes, start_time, safety)
                                    # plot_simulation!(sim_queues, sim_arcs, kind_sim, config, kind_opt, kind_queue, minutes_in_period, past_minutes, start_time, safety, max_enter, nodes)

                                catch e
                                    println("Error in combination $combination_count: $e")
                                    println(log_file, "  ✗ Error: $e")
                                end

                                flush(log_file)
                                next!(progress)
                            end
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

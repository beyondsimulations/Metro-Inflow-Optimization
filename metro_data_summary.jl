# load packages
using Pkg
Pkg.activate("metroflow")

using CSV
using DataFrames
using Dates
using Plots
using StatsPlots
using Measures
using ColorSchemes
using Statistics
using Colors
using TOML
using ProgressMeter

include("functions/config.jl")

# =============================================================================
# CONFIGURATION - Supports command line --config argument
# =============================================================================
# Parse --config argument (same pattern as metro_framework_parallel.jl)
function parse_summary_args()
    config_path = "config/doha.toml"  # Default
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--config" && i + 1 <= length(ARGS)
            config_path = ARGS[i + 1]
            i += 2
        else
            i += 1
        end
    end
    return config_path
end

CONFIG_FILE = parse_summary_args()

# Load region configuration
config = load_config(CONFIG_FILE)
println("Analyzing data for: $(config.display_name)")

# Auto-select example station based on region
ANALYSIS_STATION = if config.name == "shanghai"
    "Metro_LongcaoRoad"
elseif config.name == "doha"
    "Metro_Lusail"
else
    "Metro_Central"  # Fallback
end
println("Analysis station: $ANALYSIS_STATION")

# =============================================================================
# Helper Functions
# =============================================================================

# Helper function to compute changesort (hours after midnight for sorting)
function compute_changesort(hour_val, closed_hours)
    max_closed = maximum(closed_hours)
    if hour_val <= max_closed
        return hour_val + 24
    else
        return hour_val
    end
end

# Compute which changesort values to filter (closed hours)
function get_closed_changesort(closed_hours)
    return [h <= 5 ? h + 24 : h for h in closed_hours]
end

# Assign daycode based on operating day (spans midnight)
function assign_daycode(date, hour_val, unique_dates, first_open_hour, last_closed_hour)
    for i in 1:length(unique_dates)-1
        current_date = unique_dates[i]
        next_date = unique_dates[i+1]
        if (date == current_date && hour_val >= first_open_hour) ||
           (date == next_date && hour_val <= last_closed_hour)
            return string(current_date)
        end
    end
    return "irrelevant"
end

# Parse all parameters from filename
# Format: sim_arcs_region_date_safety_period_past_maxenter_minenter_scaling_kindopt_kindsim_kindqueue.csv
# Example: sim_arcs_doha_2022-11-27_0.9_60_240_160_10_1.0_linweight_bound_lag_periods.csv
struct FileParams
    safety::Union{Nothing, Float64}
    period_length::Union{Nothing, Int}
    past_minutes::Union{Nothing, Int}
    max_enter::Union{Nothing, Int}
    min_enter::Union{Nothing, Int}
    scaling::Union{Nothing, Float64}
    is_baseline::Bool
end

function parse_file_params(filename::String)::FileParams
    is_baseline = occursin("unbound", filename)

    # Use regex to extract parameters
    # Pattern: _safety_period_past_maxenter_minenter_scaling_kindopt_
    m = match(r"_(\d+\.?\d*)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+\.?\d*)_(?:linweight|none)_", filename)

    if m !== nothing
        return FileParams(
            tryparse(Float64, m.captures[1]),  # safety
            tryparse(Int, m.captures[2]),       # period_length
            tryparse(Int, m.captures[3]),       # past_minutes
            tryparse(Int, m.captures[4]),       # max_enter
            tryparse(Int, m.captures[5]),       # min_enter
            tryparse(Float64, m.captures[6]),   # scaling
            is_baseline
        )
    end

    return FileParams(nothing, nothing, nothing, nothing, nothing, nothing, is_baseline)
end

# Check if file matches plotting filters from config
function matches_plot_filters(params::FileParams, config)::Bool
    # Always include baseline files (they don't have optimization parameters)
    if params.is_baseline
        # Only check scaling for baseline files
        if params.scaling !== nothing && !(params.scaling in config.scaling_factors)
            return false
        end
        return true
    end

    # For optimized files, check all plotting filters
    if config.plot_safety !== nothing && params.safety != config.plot_safety
        return false
    end
    if config.plot_period_length !== nothing && params.period_length != config.plot_period_length
        return false
    end
    # Check past_periods (computed from past_minutes / period_length)
    if config.plot_past_periods !== nothing && params.period_length !== nothing && params.past_minutes !== nothing
        computed_past_periods = div(params.past_minutes, params.period_length)
        if computed_past_periods != config.plot_past_periods
            return false
        end
    end
    if config.plot_max_enter !== nothing && params.max_enter != config.plot_max_enter
        return false
    end
    if config.plot_min_enter !== nothing && params.min_enter != config.plot_min_enter
        return false
    end
    # Also check scaling is in allowed list
    if params.scaling !== nothing && !(params.scaling in config.scaling_factors)
        return false
    end
    return true
end

# Extract just the scaling factor (for grouping)
function extract_scaling_from_filename(filename::String)
    m = match(r"_(\d+\.?\d*)_(?:linweight|none)", filename)
    if m !== nothing
        return tryparse(Float64, m.captures[1])
    end
    return nothing
end

# Extended color palette
const PLOT_COLORS = [:firebrick, :cyan4, :seashell4, :purple, :orange, :green, :blue,
                     :brown, :pink, :teal, :gold, :coral, :navy, :olive]
const PLOT_LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

# =============================================================================
# Load demand data from configured region
# =============================================================================
println("\n" * "="^60)
println("Loading demand data...")
println("="^60)

od_files = filter(f -> startswith(f, "OD_") && endswith(f, ".csv"), readdir(config.base_dir))
od_dfs = DataFrame[]
@showprogress "Loading OD files: " for file in od_files
    push!(od_dfs, CSV.read(joinpath(config.base_dir, file), DataFrame))
end
data = reduce(vcat, od_dfs)
println("Loaded $(nrow(data)) demand records from $(length(od_files)) files")

data.hour = hour.(data.datetime)

# Dynamically assign daycode based on actual dates in data
unique_dates = sort(unique(data.date))
println("Found dates: ", unique_dates)

# Determine first open hour (hour after last closed hour)
first_open_hour = maximum(config.closed_hours) + 1
last_closed_hour = minimum(config.closed_hours)

# Vectorized daycode assignment
data.daycode = assign_daycode.(data.date, data.hour, Ref(unique_dates), first_open_hour, last_closed_hour)

# Overall demand of each interval
gpd_demand = groupby(data, :daycode)
gpd_demand = combine(gpd_demand, :value => sum)

# Peak for all stations
gpd_peak = groupby(data, [:daycode,:datetime,:hour])
gpd_peak = combine(gpd_peak, :value => sum)
gpd_peak.changesort = compute_changesort.(gpd_peak.hour, Ref(config.closed_hours))
sort!(gpd_peak, [:changesort, :datetime])
gpd_peak.stringday .= string.(Time.(gpd_peak.datetime))
filter!(row -> !(row.daycode == "irrelevant"), gpd_peak)
gpd_peak.day .= Date.(gpd_peak.daycode)

# Dynamic styles based on number of unique days
n_days = length(unique(gpd_peak.day))
day_linestyles = repeat(PLOT_LINESTYLES, ceil(Int, n_days / length(PLOT_LINESTYLES)))[1:n_days]
day_colors = repeat(PLOT_COLORS, ceil(Int, n_days / length(PLOT_COLORS)))[1:n_days]
day_markers = repeat([:rect, :star5, :xcross, :circle, :diamond, :utriangle], ceil(Int, n_days / 6))[1:n_days]

@df gpd_peak plot(
    :stringday,
    :value_sum,
    group=:day,
    xlabel="time",
    ylabel="demand",
    linestyle=permutedims(day_linestyles),
    markershape=permutedims(day_markers),
    markersize = 2,
    linewidth = 0.8,
    markerstrokewidth = 0,
    xrotation = 20,
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    color=permutedims(day_colors),
)
savefig("visuals/transportdemand_$(config.name).pdf")
println("Saved: visuals/transportdemand_$(config.name).pdf")


# Peak for individual stations
station = ANALYSIS_STATION
gpd_peak_station = filter(row -> row.origin == station, data)
if nrow(gpd_peak_station) > 0
    gpd_peak_station = groupby(gpd_peak_station, [:daycode,:datetime,:hour])
    gpd_peak_station = combine(gpd_peak_station, :value => sum)
    gpd_peak_station.changesort = compute_changesort.(gpd_peak_station.hour, Ref(config.closed_hours))
    sort!(gpd_peak_station, [:changesort, :datetime])
    gpd_peak_station.stringday .= string.(Time.(gpd_peak_station.datetime))
    filter!(row -> !(row.daycode == "irrelevant"), gpd_peak_station)
    gpd_peak_station.day .= Date.(gpd_peak_station.daycode)

    n_days_station = length(unique(gpd_peak_station.day))
    station_linestyles = repeat(PLOT_LINESTYLES, ceil(Int, n_days_station / length(PLOT_LINESTYLES)))[1:n_days_station]
    station_colors = repeat(PLOT_COLORS, ceil(Int, n_days_station / length(PLOT_COLORS)))[1:n_days_station]
    station_markers = repeat([:rect, :star5, :xcross, :circle, :diamond], ceil(Int, n_days_station / 5))[1:n_days_station]

    @df gpd_peak_station plot(
        :stringday,
        :value_sum,
        group=:day,
        xlabel="time",
        ylabel="demand",
        linestyle=permutedims(station_linestyles),
        markershape=permutedims(station_markers),
        markersize = 2,
        linewidth = 0.8,
        markerstrokewidth = 0,
        xrotation = 20,
        size=(700,350),
        margin=5mm,
        fontfamily="Computer Modern",
        color=permutedims(station_colors),
    )
    savefig("visuals/transportdemand_$(config.name)_$station.pdf")
    println("Saved: visuals/transportdemand_$(config.name)_$station.pdf")
else
    println("Warning: No data found for station $station")
end

# =============================================================================
# Aggregate analysis from logfiles (adapted from data_paper analysis)
# =============================================================================
println("\n" * "="^60)
println("Aggregate Analysis from Logfiles")
println("="^60)

# Find all logfiles for this region
logfile_pattern = "logfile_$(config.name)_"
logfiles = filter(f -> startswith(f, logfile_pattern) && endswith(f, ".csv"), readdir("results"))
println("Found $(length(logfiles)) logfiles for $(config.name)")

if length(logfiles) > 0
    log_dfs = DataFrame[]
    @showprogress "Loading logfiles: " for file in logfiles
        push!(log_dfs, CSV.read("results/$file", DataFrame))
    end
    logdata = reduce(vcat, log_dfs)

    # Filter to config parameters only (allows keeping old results while analyzing current config)
    if nrow(logdata) > 0
        valid_safety = Set(config.safety_factors)
        valid_min = Set(Float64.(config.min_enter))
        valid_max = Set(Float64.(config.max_enter))
        valid_scaling = Set(config.scaling_factors)

        before_filter = nrow(logdata)
        filter!(row ->
            row.safety in valid_safety &&
            row.min_enter in valid_min &&
            row.max_enter in valid_max &&
            row.scaling in valid_scaling,
            logdata
        )
        println("Filtered to config params: $(nrow(logdata)) rows (from $before_filter)")
    end

    if nrow(logdata) > 0 && hasproperty(logdata, :kind_simulation)
        # Analysis by simulation parameters
        ana_util_entry_day = groupby(logdata, [
            :kind_simulation,
            :min_enter,
            :max_enter,
            :safety,
        ])

        ana_util_entry_day = combine(ana_util_entry_day,
            :safety_minutes => mean,
            :exceeded_minutes => mean,
            :max_utilization => mean,
            :avg_utilization => mean,
            :avg_queue => mean,
            :people_moved => mean,
        )
        sort!(ana_util_entry_day, :max_enter)
        CSV.write("results/analysis_util_entry_$(config.name).csv", ana_util_entry_day)
        println("Saved: results/analysis_util_entry_$(config.name).csv")

        # Analysis by period length
        if hasproperty(logdata, :past_minutes) && hasproperty(logdata, :period_length)
            logdata.past_periods .= ceil.(Int, logdata.past_minutes ./ logdata.period_length)
            length_periods_day = groupby(logdata, [
                :kind_simulation,
                :period_length,
                :past_periods,
            ])

            length_periods_day = combine(length_periods_day,
                :avg_duration => mean,
                :safety_minutes => mean,
                :exceeded_minutes => mean,
                :max_utilization => mean,
                :avg_utilization => mean,
                :avg_queue => mean,
                :people_moved => mean,
            )
            CSV.write("results/length_periods_$(config.name).csv", length_periods_day)
            println("Saved: results/length_periods_$(config.name).csv")
        end
    else
        println("No valid logdata found or missing kind_simulation column")
    end
else
    println("No logfiles found - skipping aggregate analysis")
end

# =============================================================================
# Detailed Plots arc utilization - PER SCALING FACTOR
# =============================================================================
println("\n" * "="^60)
println("Arc Utilization Plots (per scaling factor)")
println("="^60)

dir = "results/arcs"

# Only process files for the current region
all_region_files = filter(f -> startswith(f, "sim_arcs_$(config.name)") && endswith(f, ".csv"), readdir(dir))
println("Found $(length(all_region_files)) total arc files for $(config.name)")

# Print active plotting filters
if any(x -> x !== nothing, [config.plot_safety, config.plot_period_length, config.plot_past_periods, config.plot_max_enter, config.plot_min_enter])
    println("Active plotting filters:")
    config.plot_safety !== nothing && println("  - safety: $(config.plot_safety)")
    config.plot_period_length !== nothing && println("  - period_length: $(config.plot_period_length)")
    config.plot_past_periods !== nothing && println("  - past_periods: $(config.plot_past_periods)")
    config.plot_max_enter !== nothing && println("  - max_enter: $(config.plot_max_enter)")
    config.plot_min_enter !== nothing && println("  - min_enter: $(config.plot_min_enter)")
end

# Filter files by config parameters (scaling + optional plotting filters)
region_files = String[]
file_scaling_map = Dict{String, Float64}()
file_params_map = Dict{String, FileParams}()

println("Filtering files by config parameters...")
@showprogress "Parsing filenames: " for file in all_region_files
    params = parse_file_params(file)
    if matches_plot_filters(params, config) && params.scaling !== nothing
        push!(region_files, file)
        file_scaling_map[file] = params.scaling
        file_params_map[file] = params
    end
end
println("Filtered to $(length(region_files)) files matching config parameters")

# Colors for the 3 days: red, blue, grey
day_colors = [:firebrick, :teal, :grey]

# Process each scaling factor separately
for target_scaling in sort(collect(config.scaling_factors))
    println("\n--- Processing scaling factor: $target_scaling ---")

    scaling_files = filter(f -> get(file_scaling_map, f, -1.0) == target_scaling, region_files)
    println("  Files for this scaling: $(length(scaling_files))")

    if isempty(scaling_files)
        println("  No files found for scaling $target_scaling, skipping...")
        continue
    end

    arc_dfs = DataFrame[]
    @showprogress "  Loading arc files (scaling=$target_scaling): " for file in scaling_files
        new_data = CSV.read("$dir/$(file)", DataFrame)
        new_data.kind .= occursin("unbound", file) ? "baseline" : "optimized"
        new_data.scaling .= target_scaling
        new_data.date .= hasproperty(new_data, :datetime) && nrow(new_data) > 0 ? string(Date(new_data.datetime[1])) : "unknown"
        push!(arc_dfs, new_data)
    end
    data_scaling = isempty(arc_dfs) ? DataFrame() : reduce(vcat, arc_dfs)

    if nrow(data_scaling) == 0
        println("  No data loaded for scaling $target_scaling")
        continue
    end

    arcutil = groupby(data_scaling, [:kind,:date,:datetime])
    arcutil = combine(arcutil, :utilization => maximum)
    arcutil.combination .= arcutil.date .* " " .* arcutil.kind
    arcutil.changesort = compute_changesort.(hour.(arcutil.datetime), Ref(config.closed_hours))
    sort!(arcutil, [:changesort, :datetime])
    arcutil.stringday .= string.(Time.(arcutil.datetime))

    # Filter out closed hours using config
    closed_changesort = get_closed_changesort(config.closed_hours)
    filter!(row -> !(row.changesort in closed_changesort), arcutil)

    if nrow(arcutil) == 0
        println("  No data after filtering for scaling $target_scaling")
        continue
    end

    # Build color and linestyle arrays based on combinations (date + kind)
    local unique_dates = sort(unique(arcutil.date))
    local unique_combinations = sort(unique(arcutil.combination))

    local comb_colors = Symbol[]
    local comb_linestyles = Symbol[]
    for comb in unique_combinations
        local day_idx = findfirst(d -> occursin(d, comb), unique_dates)
        day_idx = day_idx === nothing ? 1 : day_idx
        push!(comb_colors, day_colors[mod1(day_idx, length(day_colors))])
        if occursin("baseline", comb)
            push!(comb_linestyles, :solid)
        else
            push!(comb_linestyles, :dot)
        end
    end

    # Combined plot for this scaling factor
    local scaling_label = target_scaling == 1.0 ? "observed" : "+$(round(Int, (target_scaling - 1) * 100))%"

    @df arcutil plot(
        :stringday,
        :utilization_maximum,
        group=:combination,
        xlabel="time",
        ylabel="maximal arc utilization",
        linestyle=permutedims(comb_linestyles),
        linewidth = 1,
        xrotation = 20,
        size=(700,350),
        margin=5mm,
        fontfamily="Computer Modern",
        color=permutedims(comb_colors),
        legend = :topleft,
    )
    local scaling_suffix = replace(string(target_scaling), "." => "_")
    savefig("visuals/arcutil_$(config.name)_scaling_$(scaling_suffix).pdf")
    println("  Saved: visuals/arcutil_$(config.name)_scaling_$(scaling_suffix).pdf")

    # Per-day plots for this scaling factor
    for (y, x) in enumerate(unique_dates)
        colorset = day_colors[mod1(y, length(day_colors))]
        day = filter(row -> row.date == x, arcutil)

        kind_linestyles = [k == "baseline" ? :solid : :dot for k in sort(unique(day.kind))]

        @df day plot(
            :stringday,
            :utilization_maximum,
            group=:kind,
            xlabel="time",
            ylabel="maximal arc utilization",
            linestyle=permutedims(kind_linestyles),
            linewidth = 1.2,
            xrotation = 20,
            size=(700,280),
            margin=5mm,
            fontfamily="Computer Modern",
            legend = :topleft,
            color = colorset,
            ylim = (0, max(1.5, maximum(arcutil.utilization_maximum))),
        )
        savefig("visuals/arcutil_$(config.name)_$(x)_scaling_$(scaling_suffix).pdf")
        println("  Saved: visuals/arcutil_$(config.name)_$(x)_scaling_$(scaling_suffix).pdf")
    end
end


# =============================================================================
# Detailed Plots waiting time - PER SCALING FACTOR
# =============================================================================
println("\n" * "="^60)
println("Queue Plots (per scaling factor)")
println("="^60)

dir = "results/queues"

# Only process files for the current region
all_queue_files = filter(f -> startswith(f, "sim_queues_$(config.name)") && endswith(f, ".csv"), readdir(dir))
println("Found $(length(all_queue_files)) total queue files for $(config.name)")

# Filter files by config parameters (scaling + optional plotting filters)
queue_files = String[]
queue_scaling_map = Dict{String, Float64}()

println("Filtering files by config parameters...")
@showprogress "Parsing queue filenames: " for file in all_queue_files
    params = parse_file_params(file)
    if matches_plot_filters(params, config) && params.scaling !== nothing
        push!(queue_files, file)
        queue_scaling_map[file] = params.scaling
    end
end
println("Filtered to $(length(queue_files)) queue files matching config parameters")

# Process each scaling factor separately
for target_scaling in sort(collect(config.scaling_factors))
    println("\n--- Processing queues for scaling factor: $target_scaling ---")

    scaling_files = filter(f -> get(queue_scaling_map, f, -1.0) == target_scaling, queue_files)
    println("  Files for this scaling: $(length(scaling_files))")

    if isempty(scaling_files)
        println("  No queue files found for scaling $target_scaling, skipping...")
        continue
    end

    queue_dfs = DataFrame[]
    @showprogress "  Loading queue files (scaling=$target_scaling): " for file in scaling_files
        new_data = CSV.read("$dir/$(file)", DataFrame)
        new_data.kind .= occursin("unbound", file) ? "baseline" : "optimized"
        new_data.scaling .= target_scaling
        new_data.date .= hasproperty(new_data, :datetime) && nrow(new_data) > 0 ? string(Date(new_data.datetime[1])) : "unknown"
        push!(queue_dfs, new_data)
    end
    data_scaling = isempty(queue_dfs) ? DataFrame() : reduce(vcat, queue_dfs)

    if nrow(data_scaling) == 0
        println("  No queue data loaded for scaling $target_scaling")
        continue
    end

    queues = groupby(data_scaling, [:kind,:date,:datetime])
    queues = combine(queues, :queued => mean)
    queues.combination .= queues.date .* " " .* queues.kind
    queues.changesort = compute_changesort.(hour.(queues.datetime), Ref(config.closed_hours))
    sort!(queues, [:changesort, :datetime])
    queues.stringday .= string.(Time.(queues.datetime))

    # Filter out closed hours
    closed_changesort = get_closed_changesort(config.closed_hours)
    filter!(row -> !(row.changesort in closed_changesort), queues)

    if nrow(queues) == 0
        println("  No queue data after filtering for scaling $target_scaling")
        continue
    end

    # Dynamic styles
    local unique_dates = sort(unique(queues.date))
    local unique_combinations = sort(unique(queues.combination))

    local comb_colors = Symbol[]
    local comb_linestyles = Symbol[]
    for comb in unique_combinations
        local day_idx = findfirst(d -> occursin(d, comb), unique_dates)
        day_idx = day_idx === nothing ? 1 : day_idx
        push!(comb_colors, day_colors[mod1(day_idx, length(day_colors))])
        if occursin("baseline", comb)
            push!(comb_linestyles, :solid)
        else
            push!(comb_linestyles, :dot)
        end
    end

    local scaling_label = target_scaling == 1.0 ? "observed" : "+$(round(Int, (target_scaling - 1) * 100))%"
    local scaling_suffix = replace(string(target_scaling), "." => "_")

    @df queues plot(
        :stringday,
        :queued_mean,
        group=:combination,
        xlabel="time",
        ylabel="mean queue",
        linestyle=permutedims(comb_linestyles),
        linewidth = 1,
        xrotation = 20,
        size=(700,350),
        margin=5mm,
        fontfamily="Computer Modern",
        color=permutedims(comb_colors),
        legend = :topleft,
    )
    savefig("visuals/queues_$(config.name)_scaling_$(scaling_suffix).pdf")
    println("  Saved: visuals/queues_$(config.name)_scaling_$(scaling_suffix).pdf")

    # Per-day queue plots
    for (idx, x) in enumerate(unique_dates)
        day = filter(row -> row.date == x, queues)
        colorset = day_colors[mod1(idx, length(day_colors))]

        kind_linestyles = [k == "baseline" ? :solid : :dot for k in sort(unique(day.kind))]

        @df day plot(
            :stringday,
            :queued_mean,
            group=:kind,
            xlabel="time",
            ylabel="mean queue",
            linestyle=permutedims(kind_linestyles),
            linewidth = 1,
            xrotation = 20,
            size=(700,350),
            margin=5mm,
            fontfamily="Computer Modern",
            color=colorset,
            legend = :topleft,
        )
        savefig("visuals/queues_$(config.name)_$(x)_scaling_$(scaling_suffix).pdf")
        println("  Saved: visuals/queues_$(config.name)_$(x)_scaling_$(scaling_suffix).pdf")
    end
end

println("\n" * "="^60)
println("Analysis complete for $(config.display_name)")
println("="^60)

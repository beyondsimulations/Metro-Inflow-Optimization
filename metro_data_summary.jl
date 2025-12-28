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

# Helper function to compute changesort (hours after midnight for sorting)
function compute_changesort(hour_val, closed_hours)
    # Hours in closed_hours that are early morning get shifted to 24+
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

# Extended color palette
const PLOT_COLORS = [:firebrick, :cyan4, :seashell4, :purple, :orange, :green, :blue,
                     :brown, :pink, :teal, :gold, :coral, :navy, :olive]
const PLOT_LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

# =============================================================================
# Load demand data from configured region
# =============================================================================
data = DataFrame()
od_files = filter(f -> startswith(f, "OD_") && endswith(f, ".csv"), readdir(config.base_dir))
for file in od_files
    new_data = CSV.read(joinpath(config.base_dir, file), DataFrame)
    global data = vcat(data, new_data)
end
println("Loaded $(nrow(data)) demand records from $(length(od_files)) files")

data.hour = hour.(data.datetime)
data.daycode .= "irrelevant"

# Dynamically assign daycode based on actual dates in data
# A "day" runs from first open hour to last closed hour
unique_dates = sort(unique(data.date))
println("Found dates: ", unique_dates)

# Determine first open hour (hour after last closed hour)
first_open_hour = maximum(config.closed_hours) + 1
last_closed_hour = minimum(config.closed_hours)

for row in eachrow(data)
    for i in 1:length(unique_dates)-1
        current_date = unique_dates[i]
        next_date = unique_dates[i+1]
        # Day runs from first_open_hour on current_date to last_closed_hour on next_date
        if (row.date == current_date && row.hour >= first_open_hour) ||
           (row.date == next_date && row.hour <= last_closed_hour)
            row.daycode = string(current_date)
            break
        end
    end
end

# Overall demand of each interval
gpd_demand = groupby(data, :daycode)
gpd_demand = combine(gpd_demand, :value => sum)

# Peak for all stations
gpd_peak = groupby(data, [:daycode,:datetime,:hour])
gpd_peak = combine(gpd_peak, :value => sum)
gpd_peak.changesort .= 0
for row in eachrow(gpd_peak)
    row.changesort = compute_changesort(row.hour, config.closed_hours)
end
sort!(gpd_peak,:datetime)
sort!(gpd_peak,:changesort)
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
    gpd_peak_station.changesort .= 0
    for row in eachrow(gpd_peak_station)
        row.changesort = compute_changesort(row.hour, config.closed_hours)
    end
    sort!(gpd_peak_station,:datetime)
    sort!(gpd_peak_station,:changesort)
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
println("\n--- Aggregate Analysis from Logfiles ---")

# Find all logfiles for this region
logfile_pattern = "logfile_$(config.name)_"
logfiles = filter(f -> startswith(f, logfile_pattern) && endswith(f, ".csv"), readdir("results"))
println("Found $(length(logfiles)) logfiles for $(config.name)")

if length(logfiles) > 0
    logdata = DataFrame()
    for file in logfiles
        new_data = CSV.read("results/$file", DataFrame)
        global logdata = vcat(logdata, new_data)
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
# Detailed Plots arc utilization
# =============================================================================
println("\n--- Arc Utilization Plots ---")

data = DataFrame()
dir = "results/arcs"

# Only process files for the current region
region_files = filter(f -> startswith(f, "sim_arcs_$(config.name)") && endswith(f, ".csv"), readdir(dir))
println("Found $(length(region_files)) arc result files for $(config.name)")

for file in region_files
    new_data = CSV.read("$dir/$(file)", DataFrame)
    new_data.kind .= "optimized"
    if occursin("unbound", file)
        new_data.kind .= "baseline"
    end
    # Extract date from datetime column
    if hasproperty(new_data, :datetime) && nrow(new_data) > 0
        new_data.date .= string(Date(new_data.datetime[1]))
    else
        new_data.date .= "unknown"
    end
    global data = vcat(data, new_data)
end

if nrow(data) > 0
    arcutil = groupby(data, [:kind,:date,:datetime])
    arcutil = combine(arcutil, :utilization => maximum)
    arcutil.combination .= arcutil.date .* " " .* arcutil.kind

    arcutil.changesort .= 0
    for row in eachrow(arcutil)
        row.changesort = compute_changesort(hour(row.datetime), config.closed_hours)
    end
    sort!(arcutil, :datetime)
    sort!(arcutil, :changesort)
    arcutil.stringday .= string.(Time.(arcutil.datetime))

    # Filter out closed hours using config
    closed_changesort = get_closed_changesort(config.closed_hours)
    filter!(row -> !(row.changesort in closed_changesort), arcutil)

    # Dynamic styles based on combinations
    n_combinations = length(unique(arcutil.combination))
    comb_linestyles = repeat(PLOT_LINESTYLES, ceil(Int, n_combinations / length(PLOT_LINESTYLES)))[1:n_combinations]
    comb_colors = repeat(PLOT_COLORS, ceil(Int, n_combinations / length(PLOT_COLORS)))[1:n_combinations]

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
    )
    savefig("visuals/arcutil_$(config.name).pdf")
    println("Saved: visuals/arcutil_$(config.name).pdf")

    # Plot per-day arc utilization
    unique_arc_dates = unique(arcutil.date)

    for (y, x) in enumerate(unique_arc_dates)
        colorset = PLOT_COLORS[mod1(y, length(PLOT_COLORS))]
        day = filter(row -> row.date == x, arcutil)

        n_kinds = length(unique(day.kind))
        kind_linestyles = repeat([:solid, :dash], ceil(Int, n_kinds / 2))[1:n_kinds]

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
            ylim = (0, maximum(arcutil.utilization_maximum)),
        )
        hline!([0.9], lw=0.5, label="restriction", linestyle=:dash, color=:black)
        savefig("visuals/arcutil_$(config.name)_$x.pdf")
        println("Saved: visuals/arcutil_$(config.name)_$x.pdf")
    end
else
    println("No arc data found for $(config.name)")
end


# =============================================================================
# Detailed Plots waiting time
# =============================================================================
println("\n--- Queue Plots ---")

data = DataFrame()
dir = "results/queues"

# Only process files for the current region
region_queue_files = filter(f -> startswith(f, "sim_queues_$(config.name)") && endswith(f, ".csv"), readdir(dir))
println("Found $(length(region_queue_files)) queue result files for $(config.name)")

for file in region_queue_files
    new_data = CSV.read("$dir/$(file)", DataFrame)
    new_data.kind .= "optimized"
    if occursin("unbound", file)
        new_data.kind .= "baseline"
    end
    # Extract date from datetime column
    if hasproperty(new_data, :datetime) && nrow(new_data) > 0
        new_data.date .= string(Date(new_data.datetime[1]))
    else
        new_data.date .= "unknown"
    end
    global data = vcat(data, new_data)
end

if nrow(data) > 0
    queues = groupby(data, [:kind,:date,:datetime])
    queues = combine(queues, :queued => mean)
    queues.combination .= queues.date .* " " .* queues.kind

    queues.changesort .= 0
    for row in eachrow(queues)
        row.changesort = compute_changesort(hour(row.datetime), config.closed_hours)
    end
    sort!(queues, :datetime)
    sort!(queues, :changesort)
    queues.stringday .= string.(Time.(queues.datetime))

    # Filter out closed hours using config
    closed_changesort = get_closed_changesort(config.closed_hours)
    filter!(row -> !(row.changesort in closed_changesort), queues)

    # Dynamic styles based on combinations
    n_combinations = length(unique(queues.combination))
    comb_linestyles = repeat(PLOT_LINESTYLES, ceil(Int, n_combinations / length(PLOT_LINESTYLES)))[1:n_combinations]
    comb_colors = repeat(PLOT_COLORS, ceil(Int, n_combinations / length(PLOT_COLORS)))[1:n_combinations]

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
    )
    savefig("visuals/queues_$(config.name).pdf")
    println("Saved: visuals/queues_$(config.name).pdf")

    # Per-day queue plots
    for (idx, x) in enumerate(unique(queues.date))
        day = filter(row -> row.date == x, queues)
        colorset = PLOT_COLORS[mod1(idx, length(PLOT_COLORS))]

        n_kinds = length(unique(day.kind))
        kind_linestyles = repeat([:solid, :dash], ceil(Int, n_kinds / 2))[1:n_kinds]

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
        )
        savefig("visuals/queues_$(config.name)_$x.pdf")
        println("Saved: visuals/queues_$(config.name)_$x.pdf")
    end
else
    println("No queue data found for $(config.name)")
end

println("\n=== Analysis complete for $(config.display_name) ===")

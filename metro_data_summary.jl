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
# CONFIGURATION - Edit these settings for your analysis
# =============================================================================
CONFIG_FILE = "config/doha.toml"      # Options: "config/doha.toml", "config/shanghai.toml"

# Load region configuration
config = load_config(CONFIG_FILE)
println("Analyzing data for: $(config.display_name)")

# Example station for individual analysis (region-specific)
# Doha: "Metro_Lusail", "Metro_Msheireb", "Metro_AlWakra"
# Shanghai: "Metro_LongcaoRoad", "Metro_JinganTemple"
ANALYSIS_STATION = "Metro_Lusail"

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
# A "day" runs from 6:00 to next day 3:00 (before metro closure)
unique_dates = sort(unique(data.date))
println("Found dates: ", unique_dates)

for row in eachrow(data)
    for i in 1:length(unique_dates)-1
        current_date = unique_dates[i]
        next_date = unique_dates[i+1]
        # Day runs from 6:00 on current_date to 3:00 on next_date
        if (row.date == current_date && row.hour >= 6) || (row.date == next_date && row.hour < 3)
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
gpd_peak = combine(gpd_peak,  :value => sum)
gpd_peak.changesort .= 0
for row in eachrow(gpd_peak)
    if row.hour <= 5
        row.changesort = row.hour + 24
    else
        row.changesort = row.hour
    end
end
sort!(gpd_peak,:datetime)
sort!(gpd_peak,:changesort)
gpd_peak.stringday .= string.(Time.(gpd_peak.datetime))
filter!(row -> !(row.daycode == "irrelevant") , gpd_peak)
gpd_peak.day .= Date.(gpd_peak.daycode)

@df gpd_peak plot(
    :stringday, 
    :value_sum, 
    group=:day, 
    xlabel="time", 
    ylabel="demand",
    linestyle=[:dash :solid :dot],
    markershape=[:rect :star5 :xcross],
    markersize = 2,
    linewidth = 0.8,
    markerstrokewidth = 0,
    xrotation = 20, 
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    color = [:firebrick :cyan4 :seashell4],
    )
savefig("visuals/transportdemand_$(config.name).pdf")


# Peak for individual stations
station = ANALYSIS_STATION
gpd_peak = filter(row -> row.origin == station, data)
gpd_peak = groupby(gpd_peak, [:daycode,:datetime,:hour])
gpd_peak = combine(gpd_peak,  :value => sum)
gpd_peak.changesort .= 0
for row in eachrow(gpd_peak)
    if row.hour <= 5
        row.changesort = row.hour + 24
    else
        row.changesort = row.hour
    end
end
sort!(gpd_peak,:datetime)
sort!(gpd_peak,:changesort)
gpd_peak.stringday .= string.(Time.(gpd_peak.datetime))
filter!(row -> !(row.daycode == "irrelevant") , gpd_peak)
gpd_peak.day .= Date.(gpd_peak.daycode)

@df gpd_peak plot(
    :stringday, 
    :value_sum, 
    group=:day, 
    xlabel="time", 
    ylabel="demand",
    linestyle=[:dash :solid :dot],
    markershape=[:rect :star5 :xcross],
    markersize = 2,
    linewidth = 0.8,
    markerstrokewidth = 0,
    xrotation = 20, 
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    color = [:firebrick :cyan4 :seashell4],
    )
savefig("visuals/transportdemand_$(config.name)_$station.pdf")

# evaluate the results of the analysis

data = DataFrame()
dir = "data_paper"
for file in readdir(dir)
    new_data = CSV.read("$dir/$(file)", DataFrame)
    global data = vcat(data,new_data)
end

# Overall demand of each interval
ana_util_entry_day = groupby(data, [
    :kind_simulation,
    :min_enter, 
    :max_enter, 
    :safety,
    ]
    )
    
ana_util_entry_day = combine(ana_util_entry_day,
    :safety_minutes => mean,
    :exceeded_minutes => mean,
    :max_utilization => mean,
    :avg_utilization => mean,
    :avg_queue => mean,
    :people_moved => mean,
    )

    sort!(ana_util_entry_day,:max_enter)

    CSV.write("results/analysis_util_entry_day.csv", ana_util_entry_day)

# Overall demand of each interval
data.past_periods .= ceil.(Int, data.past_minutes ./ data.period_length)
length_periods_day = groupby(data, [
    :kind_simulation,
    :period_length, 
    :past_periods,
    ]
    )
    
length_periods_day = combine(length_periods_day,
    :avg_duration => mean,
    :safety_minutes => mean,
    :exceeded_minutes => mean,      
    :max_utilization => mean,
    :avg_utilization => mean,
    :avg_queue => mean,
    :people_moved => mean,
    )

CSV.write("results/length_periods_day.csv", length_periods_day)

# Detailed Plots arc utilization
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

arcutil = groupby(data,[:kind,:date,:datetime])
arcutil = combine(arcutil, :utilization => maximum)
arcutil.combination .= arcutil.date .* " " .* arcutil.kind

arcutil.changesort .= 0
for row in eachrow(arcutil)
    if hour(row.datetime) <= 5
        row.changesort = hour(row.datetime) + 24
    else
        row.changesort = hour(row.datetime)
    end
end
sort!(arcutil,:datetime)
sort!(arcutil,:changesort)
arcutil.stringday .= string.(Time.(arcutil.datetime))
filter!(row -> !(row.changesort in [28,29]), arcutil)

@df arcutil plot(
    :stringday, 
    :utilization_maximum, 
    group=:combination, 
    xlabel="time", 
    ylabel="maximal arc utilization",
    linestyle=[:dot :solid :dot :solid :dot :solid],
    #markershape=[:rect :star5 :xcross],
    #markersize = 2,
    linewidth = 1,
    #markerstrokewidth = 0,
    xrotation = 20, 
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    color = [:firebrick :firebrick :cyan4 :cyan4 :seashell4 :seashell4],
    )
savefig("visuals/arcutil_$(config.name).pdf")

# Plot per-day arc utilization
unique_arc_dates = unique(arcutil.date)
colors = [:firebrick, :cyan4, :seashell4, :purple, :orange, :green, :blue]

for (y, x) in enumerate(unique_arc_dates)
    colorset = colors[mod1(y, length(colors))]
    day = filter(row -> row.date == x, arcutil)

    @df day plot(
        :stringday,
        :utilization_maximum,
        group=:kind,
        xlabel="time",
        ylabel="maximal arc utilization",
        linestyle=[:dot :solid :dot :solid :dot :solid],
        linewidth = 1.2,
        xrotation = 20,
        size=(700,280),
        margin=5mm,
        fontfamily="Computer Modern",
        legend = :topleft,
        color = colorset,
        ylim = (0,maximum(arcutil.utilization_maximum)),
        )
        hline!([0.9],lw=0.5,label = "restriction", linestyle = :dash, color = :black)
    savefig("visuals/arcutil_$(config.name)_$x.pdf")
end


# Detailed Plots waiting time
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

queues = groupby(data,[:kind,:date,:datetime])
queues = combine(queues, :queued => mean)
queues.combination .= queues.date .* " " .* queues.kind

queues.changesort .= 0
for row in eachrow(queues)
    if hour(row.datetime) <= 5
        row.changesort = hour(row.datetime) + 24
    else
        row.changesort = hour(row.datetime)
    end
end
sort!(queues,:datetime)
sort!(queues,:changesort)
queues.stringday .= string.(Time.(queues.datetime))
filter!(row -> !(row.changesort in [28,29]), queues)

@df queues plot(
        :stringday, 
        :queued_mean, 
        group=:combination, 
        xlabel="time", 
        ylabel="maximal queue",
        linestyle=[:dot :solid :dot :solid :dot :solid],
        #markershape=[:rect :star5 :xcross],
        #markersize = 2,
        linewidth = 1,
        #markerstrokewidth = 0,
        xrotation = 20, 
        size=(700,350),
        margin=5mm,
        fontfamily="Computer Modern",
        color = [:firebrick :firebrick :cyan4 :cyan4 :seashell4 :seashell4],
        )
    savefig("visuals/queues_$(config.name).pdf")

for x in unique(queues.date)

    day = filter(row -> row.date == x, queues)

    @df day plot(
        :stringday,
        :queued_mean,
        group=:kind,
        xlabel="time",
        ylabel="maximal queue",
        linestyle=[:dot :solid :dot :solid :dot :solid],
        linewidth = 1,
        xrotation = 20,
        size=(700,350),
        margin=5mm,
        fontfamily="Computer Modern",
        color = [:firebrick :firebrick :cyan4 :cyan4 :seashell4 :seashell4],
        )
    savefig("visuals/queues_$(config.name)_$x.pdf")
end

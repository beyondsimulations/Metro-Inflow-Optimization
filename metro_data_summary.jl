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

# evaluate input demand data
data = DataFrame()
for file in readdir("data_demand")
    new_data = CSV.read("data_demand/$(file)", DataFrame)
    global data = vcat(data,new_data)
end

data.hour = hour.(data.datetime)
data.hour = hour.(data.datetime)
data.daycode .= "irrelevant"
for row in eachrow(data)
    if (row.date == Date("2022-11-27") && row.hour >= 6) || (row.date == Date("2022-11-28") && row.hour < 3)
        row.daycode = "2022-11-27"
    elseif (row.date == Date("2022-11-28") && row.hour >= 6) || (row.date == Date("2022-11-29") && row.hour < 3)
        row.daycode = "2022-11-28"
    elseif (row.date == Date("2022-11-29") && row.hour >= 6) || (row.date == Date("2022-11-30") && row.hour < 3)
        row.daycode = "2022-11-29"
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
    linewidth = 0.5,
    markerstrokewidth = 0,
    xrotation = 20, 
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    palette = palette([:black, :grey], 3)
    )
savefig("visuals/transportdemand.pdf")


# Peak for individual stations
station = "Metro_Lusail"
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
    linewidth = 0.5,
    markerstrokewidth = 0,
    xrotation = 20, 
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    palette = palette([:black, :grey], 3)
    )
savefig("visuals/transportdemand$station.pdf")

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

    CSV.write("results_paper/ana_util_entry_day.csv", ana_util_entry_day)

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

CSV.write("results_paper/length_periods_day.csv", length_periods_day)


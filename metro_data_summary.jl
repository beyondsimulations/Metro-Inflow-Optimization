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
    linewidth = 0.8,
    markerstrokewidth = 0,
    xrotation = 20, 
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    color = [:firebrick :cyan4 :seashell4],
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
    linewidth = 0.8,
    markerstrokewidth = 0,
    xrotation = 20, 
    size=(700,350),
    margin=5mm,
    fontfamily="Computer Modern",
    color = [:firebrick :cyan4 :seashell4],
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

# Detailed Plots arc utilization
data = DataFrame()
dir = "results_arcs"
for file in readdir(dir)
    new_data = CSV.read("$dir/$(file)", DataFrame)
    new_data.kind .= "optimized"
    if occursin("unbound",file)
        new_data.kind .= "baseline"
    end
    new_data.date .= "empty"
    if occursin("2022-11-27T05:00:00_",file)
        new_data.date  .= "2022-11-27"
    end
    if occursin("2022-11-28T05:00:00_",file)
        new_data.date  .= "2022-11-28"
    end
    if occursin("2022-11-29T05:00:00_",file)
        new_data.date  .= "2022-11-29"
    end
    global data = vcat(data,new_data)
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
savefig("visuals/arcutil.pdf")

for y in 1:3

    x = unique(arcutil.date)[y]
    colorset = [:firebrick,:cyan4,:seashell4][y]

    day = filter(row -> row.date == x, arcutil)

    @df day plot(
        :stringday, 
        :utilization_maximum, 
        group=:kind, 
        xlabel="time", 
        ylabel="maximal arc utilization",
        linestyle=[:dot :solid :dot :solid :dot :solid],
        #markershape=[:rect :star5 :xcross],
        #markersize = 2,
        linewidth = 1.2,
        #markerstrokewidth = 0,
        xrotation = 20, 
        size=(700,280),
        margin=5mm,
        fontfamily="Computer Modern",
        legend = :topleft,
        color = colorset,
        ylim = (0,maximum(arcutil.utilization_maximum)),
        )
        hline!([0.9],lw=0.5,label = "restriction", linestyle = :dash, color = :black)
    savefig("visuals/arcutil_$x.pdf")
end


# Detailed Plots waiting time
data = DataFrame()
dir = "results_queues"
for file in readdir(dir)
    new_data = CSV.read("$dir/$(file)", DataFrame)
    new_data.kind .= "optimized"
    if occursin("unbound",file)
        new_data.kind .= "baseline"
    end
    new_data.date .= "empty"
    if occursin("2022-11-27T05:00:00_",file)
        new_data.date  .= "2022-11-27"
    end
    if occursin("2022-11-28T05:00:00_",file)
        new_data.date  .= "2022-11-28"
    end
    if occursin("2022-11-29T05:00:00_",file)
        new_data.date  .= "2022-11-29"
    end
    global data = vcat(data,new_data)
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
    savefig("visuals/queues.pdf")

for x in unique(queues.date)

    day = filter(row -> row.date == x, queues)

    @df day plot(
        :stringday, 
        :queued_mean, 
        group=:kind, 
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
    savefig("visuals/queues_$x.pdf")
end

# load packages
using Pkg
Pkg.activate("metroflow")

using CSV
using DataFrames
using Dates

data = DataFrame()
for file in readdir("data_demand")
    new_data = CSV.read("data_demand/$(file)", DataFrame)
    data = vcat(data,new_data)
end

data.hour = hour.(data.datetime)
data.daycode .= 0
for row in eachrow(data)
    if (row.date == Date("2022-11-27") && row.hour >= 6) || 
        (row.date == Date("2023-11-28") && row.hour < 6)
        row.daycode = 1
    end
    if (row.date == Date("2022-11-28") && row.hour >= 6) || 
        (row.date == Date("2023-11-29") && row.hour < 6)
        row.daycode = 2
    end
    if (row.date == Date("2022-11-29") && row.hour >= 6) || 
        (row.date == Date("2023-11-30") && row.hour < 6)
        row.daycode = 3
    end
end

gpd = groupby(data, :daycode)
gpd = combine(gpd, :value => sum)

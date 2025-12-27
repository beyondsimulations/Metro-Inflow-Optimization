"""
Shanghai OD Data Transformation Script

Transforms Shanghai metro OD flow data into the Doha-compatible format.
Run from the Shanghai data folder.

Usage:
    cd data_public/Shanghai
    julia --project=../../metroflow transform_od.jl 2017-05-10 2017-05-14   # Process Wed-Sun
    julia --project=../../metroflow transform_od.jl 2017-05-10 2017-05-16   # Process full week (Wed-Tue)
    julia --project=../../metroflow transform_od.jl 2017-05-01 2017-05-31   # Process entire month

Note: May 8-9, 2017 have incomplete OD data (<2% of usual ridership) in the source.
      Use May 10+ for complete data.

Input:
    - metroData_ODFlow.csv (11GB, in same folder)
    - stationInfo.csv (station ID → name mapping)

Output:
    - OD_2017-05-10.csv, OD_2017-05-11.csv, ... (one file per day)
    - Format: date,datetime,origin,destination,value
"""

using Pkg
Pkg.activate("../../metroflow")

using CSV
using DataFrames
using Dates
using ProgressMeter

# Parse command line arguments
function parse_args()
    if length(ARGS) != 2
        println("Usage: julia transform_od.jl <start_date> <end_date>")
        println("Example: julia transform_od.jl 2017-05-08 2017-05-14")
        exit(1)
    end

    start_date = Date(ARGS[1])
    end_date = Date(ARGS[2])

    if end_date < start_date
        println("Error: end_date must be >= start_date")
        exit(1)
    end

    return start_date, end_date
end

"""
Format station name to match metro arcs convention: Metro_StationName
Matches the format used in build_shanghai_metroarcs.jl
"""
function format_station_name(name::AbstractString)
    formatted = replace(name, " " => "")
    formatted = replace(formatted, "'" => "")
    formatted = replace(formatted, "-" => "")
    formatted = replace(formatted, "/" => "_")
    formatted = replace(formatted, "(" => "")
    formatted = replace(formatted, ")" => "")
    return "Metro_" * formatted
end

"""
Build station ID to formatted name mapping from stationInfo.csv
"""
function build_station_map(station_file::String)
    println("Loading station information...")
    stations = CSV.read(station_file, DataFrame)

    station_map = Dict{Int, String}()
    for row in eachrow(stations)
        station_id = row.stationID
        station_name = format_station_name(row.name)
        station_map[station_id] = station_name
    end

    println("Loaded $(length(station_map)) station mappings")
    return station_map
end

"""
Convert date from YYYYMMDD integer to Date
"""
function parse_date_int(date_int::Int)::Date
    year = div(date_int, 10000)
    month = div(mod(date_int, 10000), 100)
    day = mod(date_int, 100)
    return Date(year, month, day)
end

"""
Convert time from HHMMSS integer to Time
"""
function parse_time_int(time_int::Int)::Time
    hour = div(time_int, 10000)
    minute = div(mod(time_int, 10000), 100)
    second = mod(time_int, 100)
    return Time(hour, minute, second)
end

"""
Convert date integer and start time integer to DateTime string in Doha format
Format: 2017-05-08T06:00:00.0
"""
function format_datetime(date_int::Int, start_time_int::Int)::String
    date = parse_date_int(date_int)
    time = parse_time_int(start_time_int)
    dt = DateTime(date, time)
    return Dates.format(dt, "yyyy-mm-ddTHH:MM:SS.0")
end

"""
Format date for output filename and date column
Format: 2017-05-08
"""
function format_date_str(date::Date)::String
    return Dates.format(date, "yyyy-mm-dd")
end

"""
Process OD flow data and write daily output files.
Reads entire file into memory (requires ~20GB RAM for 11GB file).
"""
function process_od_data(od_file::String, station_map::Dict{Int, String},
                         start_date::Date, end_date::Date)

    # Convert dates to integer format for comparison
    start_date_int = year(start_date) * 10000 + month(start_date) * 100 + day(start_date)
    end_date_int = year(end_date) * 10000 + month(end_date) * 100 + day(end_date)

    println("Processing OD data from $(format_date_str(start_date)) to $(format_date_str(end_date))...")
    println("Loading full file into memory...")

    # Read entire file with proper types (normalizenames handles spaces in headers)
    df = CSV.read(od_file, DataFrame; normalizenames=true,
        types=Dict(:date => Int, :timeslot => Int, :startTime => Int, :endTime => Int,
                   :originStation => Int, :destinationStation => Int, :Flow => Int))

    total_rows_read = nrow(df)
    println("Loaded $(total_rows_read) rows")

    # Filter to date range
    println("Filtering to date range...")
    df = filter(row -> row.date >= start_date_int && row.date <= end_date_int, df)
    rows_in_range = nrow(df)
    println("Found $(rows_in_range) rows in date range")

    # Accumulate data by day
    daily_data = Dict{Date, DataFrame}()
    for d in start_date:Day(1):end_date
        daily_data[d] = DataFrame(
            date = String[],
            datetime = String[],
            origin = String[],
            destination = String[],
            value = Int[]
        )
    end

    # Track statistics
    rows_skipped_missing_station = 0

    # Process filtered rows
    println("Processing rows...")
    @showprogress for row in eachrow(df)
        # Look up station names
        origin_id = row.originStation
        dest_id = row.destinationStation

        if !haskey(station_map, origin_id)
            rows_skipped_missing_station += 1
            continue
        end
        if !haskey(station_map, dest_id)
            rows_skipped_missing_station += 1
            continue
        end

        origin_name = station_map[origin_id]
        dest_name = station_map[dest_id]

        # Parse date and format datetime
        date = parse_date_int(row.date)
        date_str = format_date_str(date)
        datetime_str = format_datetime(row.date, row.startTime)

        # Use Flow column as value
        flow_value = row.Flow

        # Add to daily data
        push!(daily_data[date], (date_str, datetime_str, origin_name, dest_name, flow_value))
    end

    println("\nProcessing complete:")
    println("  Total rows read: $(total_rows_read)")
    println("  Rows in date range: $(rows_in_range)")
    println("  Rows skipped (missing station): $(rows_skipped_missing_station)")

    # Write output files
    println("\nWriting output files...")
    for (date, df) in daily_data
        if nrow(df) > 0
            output_file = "OD_$(format_date_str(date)).csv"
            CSV.write(output_file, df)
            println("  Wrote $output_file ($(nrow(df)) rows)")
        else
            println("  Skipping $(format_date_str(date)) - no data")
        end
    end

    return daily_data
end

# Main execution
function main()
    start_date, end_date = parse_args()

    println("=== Shanghai OD Data Transformation ===")
    println("Date range: $(format_date_str(start_date)) to $(format_date_str(end_date))")
    println("")

    # Build station mapping
    station_map = build_station_map("stationInfo.csv")

    # Process OD data
    daily_data = process_od_data("metroData_ODFlow.csv", station_map, start_date, end_date)

    # Summary
    total_rows = sum(nrow(df) for df in values(daily_data))
    println("\n=== Summary ===")
    println("Total OD records written: $total_rows")
    println("Files created: $(count(df -> nrow(df) > 0, values(daily_data)))")
    println("Done!")
end

main()

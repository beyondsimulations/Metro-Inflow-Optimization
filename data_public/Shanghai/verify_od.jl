"""
Shanghai OD Data Verification Script

Compares generated OD files against metroData_InOutFlow.csv to verify transformation accuracy.
Run from the Shanghai data folder.

Usage:
    cd data_public/Shanghai
    julia --project=../../metroflow verify_od.jl 2017-05-10 2017-05-14

Note: May 8-9, 2017 have incomplete OD data (<2% of usual ridership) in the source.
      Use May 10-14 (Wed-Sun) for a complete sample week, or May 10-16 for Mon-Sun.

Expected: ~2% difference between OD and InOutFlow is normal. OD only contains complete
origin-destination pairs, while InOutFlow counts all taps (including orphaned entries/exits).
"""

using Pkg
Pkg.activate("../../metroflow")

using CSV
using DataFrames
using Dates
using Printf

# Parse command line arguments
function parse_args()
    if length(ARGS) != 2
        println("Usage: julia verify_od.jl <start_date> <end_date>")
        println("Example: julia verify_od.jl 2017-05-08 2017-05-14")
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
Same as in transform_od.jl
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
    stations = CSV.read(station_file, DataFrame)
    station_map = Dict{Int, String}()
    for row in eachrow(stations)
        station_id = row.stationID
        station_name = format_station_name(row.name)
        station_map[station_id] = station_name
    end
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
Load and aggregate InOutFlow data for verification
"""
function load_inoutflow(inoutflow_file::String, station_map::Dict{Int, String},
                        start_date::Date, end_date::Date)
    println("Loading InOutFlow data...")

    # Convert dates to integer format
    start_date_int = year(start_date) * 10000 + month(start_date) * 100 + day(start_date)
    end_date_int = year(end_date) * 10000 + month(end_date) * 100 + day(end_date)

    # Read full file
    df = CSV.read(inoutflow_file, DataFrame; normalizenames=true,
        types=Dict(:date => Int, :timeslot => Int, :startTime => Int, :endTime => Int,
                   :station => Int, :inFlow => Int, :outFlow => Int))

    # Filter to date range
    df = filter(row -> row.date >= start_date_int && row.date <= end_date_int, df)

    # Filter to stations we know about
    df = filter(row -> haskey(station_map, row.station), df)

    # Convert station IDs to names
    df.station_name = [station_map[s] for s in df.station]
    df.date_parsed = [parse_date_int(d) for d in df.date]

    println("  Loaded $(nrow(df)) InOutFlow records")
    return df
end

"""
Load generated OD file for a specific date
"""
function load_od_file(date::Date)
    filename = "OD_$(Dates.format(date, "yyyy-mm-dd")).csv"
    if !isfile(filename)
        return nothing
    end

    df = CSV.read(filename, DataFrame)
    return df
end

"""
Perform verification for a single date
"""
function verify_date(date::Date, inoutflow_df::DataFrame, station_map::Dict{Int, String})
    date_str = Dates.format(date, "yyyy-mm-dd")

    # Load OD file
    od_df = load_od_file(date)
    if od_df === nothing
        println("\n  $date_str: OD file not found")
        return Dict("status" => "missing", "date" => date)
    end

    # Filter InOutFlow to this date
    daily_inout = filter(row -> row.date_parsed == date, inoutflow_df)

    if nrow(daily_inout) == 0
        println("\n  $date_str: No InOutFlow data")
        return Dict("status" => "no_inout", "date" => date)
    end

    # Calculate totals from InOutFlow
    total_inflow = sum(daily_inout.inFlow)
    total_outflow = sum(daily_inout.outFlow)

    # Calculate totals from OD
    total_od = sum(od_df.value)

    # Per-station aggregations from OD
    od_by_origin = combine(groupby(od_df, :origin), :value => sum => :total)
    od_by_dest = combine(groupby(od_df, :destination), :value => sum => :total)

    # Per-station aggregations from InOutFlow
    inout_by_station = combine(groupby(daily_inout, :station_name),
        :inFlow => sum => :total_in,
        :outFlow => sum => :total_out)

    # Compare totals
    od_origins_total = sum(od_by_origin.total)
    od_dests_total = sum(od_by_dest.total)

    # Calculate differences
    outflow_diff = abs(od_origins_total - total_outflow)
    inflow_diff = abs(od_dests_total - total_inflow)
    outflow_pct = total_outflow > 0 ? 100 * outflow_diff / total_outflow : 0
    inflow_pct = total_inflow > 0 ? 100 * inflow_diff / total_inflow : 0

    # Station count comparison
    od_stations = union(Set(od_df.origin), Set(od_df.destination))
    inout_stations = Set(daily_inout.station_name)
    stations_in_od_only = setdiff(od_stations, inout_stations)
    stations_in_inout_only = setdiff(inout_stations, od_stations)

    # Print results
    println("\n  $date_str:")
    println("     InOutFlow: outFlow=$(format_number(total_outflow)), inFlow=$(format_number(total_inflow))")
    println("     OD file:   origins=$(format_number(od_origins_total)), dests=$(format_number(od_dests_total))")
    println("     Difference: outFlow=$(format_number(outflow_diff)) ($(@sprintf("%.2f", outflow_pct))%), inFlow=$(format_number(inflow_diff)) ($(@sprintf("%.2f", inflow_pct))%)")
    println("     Stations:  OD=$(length(od_stations)), InOut=$(length(inout_stations))")

    if length(stations_in_od_only) > 0 && length(stations_in_od_only) <= 5
        println("     Stations in OD only: $(join(collect(stations_in_od_only), ", "))")
    elseif length(stations_in_od_only) > 5
        println("     Stations in OD only: $(length(stations_in_od_only)) stations")
    end

    if length(stations_in_inout_only) > 0 && length(stations_in_inout_only) <= 5
        println("     Stations in InOut only: $(join(collect(stations_in_inout_only), ", "))")
    elseif length(stations_in_inout_only) > 5
        println("     Stations in InOut only: $(length(stations_in_inout_only)) stations")
    end

    # Status
    if outflow_pct < 1 && inflow_pct < 1
        println("     PASSED")
        status = "passed"
    elseif outflow_pct < 5 && inflow_pct < 5
        println("     MINOR DIFFERENCES")
        status = "minor_diff"
    else
        println("     SIGNIFICANT DIFFERENCES")
        status = "failed"
    end

    return Dict(
        "status" => status,
        "date" => date,
        "total_outflow" => total_outflow,
        "total_inflow" => total_inflow,
        "od_origins" => od_origins_total,
        "od_dests" => od_dests_total,
        "outflow_diff_pct" => outflow_pct,
        "inflow_diff_pct" => inflow_pct
    )
end

"""
Format large numbers with commas
"""
function format_number(n)
    return replace(string(round(Int, n)), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

# Main execution
function main()
    start_date, end_date = parse_args()

    println("=== Shanghai OD Verification ===")
    println("Date range: $(Dates.format(start_date, "yyyy-mm-dd")) to $(Dates.format(end_date, "yyyy-mm-dd"))")

    # Build station mapping
    station_map = build_station_map("stationInfo.csv")
    println("Loaded $(length(station_map)) station mappings")

    # Load InOutFlow data
    inoutflow_df = load_inoutflow("metroData_InOutFlow.csv", station_map, start_date, end_date)

    # Verify each date
    results = []
    println("\n--- Verification Results ---")

    for date in start_date:Day(1):end_date
        result = verify_date(date, inoutflow_df, station_map)
        push!(results, result)
    end

    # Summary
    passed = count(r -> r["status"] == "passed", results)
    minor = count(r -> r["status"] == "minor_diff", results)
    failed = count(r -> r["status"] == "failed", results)
    missing_files = count(r -> r["status"] == "missing", results)

    println("\n=== Summary ===")
    println("  Passed: $passed")
    println("  Minor differences: $minor")
    println("  Failed: $failed")
    println("  Missing OD files: $missing_files")

    # Overall assessment
    if failed == 0 && missing_files == 0
        println("\nVerification PASSED")
    elseif failed > 0
        println("\nVerification FAILED - check discrepancies above")
    else
        println("\nVerification incomplete - some files missing")
    end
end

main()

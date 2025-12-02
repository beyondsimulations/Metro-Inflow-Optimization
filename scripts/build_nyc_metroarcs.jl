using HTTP
using ZipFile
using CSV
using DataFrames
using Statistics

const GTFS_URL = "http://web.mta.info/developers/data/nyct/subway/google_transit.zip"
const DATA_DIR = "data_public/New-York/gtfs"

"""
Download and extract MTA GTFS data.
"""
function download_gtfs()
    mkpath(DATA_DIR)
    zip_path = joinpath(DATA_DIR, "google_transit.zip")

    println("Downloading GTFS data from MTA...")
    response = HTTP.get(GTFS_URL)
    open(zip_path, "w") do f
        write(f, response.body)
    end
    println("Downloaded $(round(filesize(zip_path) / 1024 / 1024, digits=2)) MB")

    println("Extracting...")
    reader = ZipFile.Reader(zip_path)
    for file in reader.files
        filepath = joinpath(DATA_DIR, file.name)
        mkpath(dirname(filepath))
        open(filepath, "w") do f
            write(f, read(file))
        end
    end
    close(reader)
    println("Extracted to $DATA_DIR")
end

"""
Load GTFS files into DataFrames.
"""
function load_gtfs()
    stops = CSV.read(joinpath(DATA_DIR, "stops.txt"), DataFrame)
    stop_times = CSV.read(joinpath(DATA_DIR, "stop_times.txt"), DataFrame)
    trips = CSV.read(joinpath(DATA_DIR, "trips.txt"), DataFrame)
    routes = CSV.read(joinpath(DATA_DIR, "routes.txt"), DataFrame)
    return stops, stop_times, trips, routes
end

"""
Parse GTFS time string (HH:MM:SS) to seconds. Handles times > 24:00.
"""
function parse_time(t::AbstractString)
    parts = split(t, ":")
    h, m, s = parse(Int, parts[1]), parse(Int, parts[2]), parse(Int, parts[3])
    return h * 3600 + m * 60 + s
end

"""
Build metro arcs from GTFS data.
"""
function build_metroarcs(; output_file="data_public/New-York/metroarcs_nyc.csv")
    println("Loading GTFS data...")
    stops, stop_times, trips, routes = load_gtfs()

    # Join stop_times with trips to get route_id
    stop_times = leftjoin(stop_times, trips[:, [:trip_id, :route_id, :direction_id]], on=:trip_id)

    # Sort by trip and stop sequence
    sort!(stop_times, [:trip_id, :stop_sequence])

    println("Extracting station connections...")

    # Build a dictionary of stop_id -> stop_name
    stop_names = Dict(row.stop_id => row.stop_name for row in eachrow(stops))

    # Build route_id -> route_short_name (line name)
    route_names = Dict(row.route_id => row.route_short_name for row in eachrow(routes))

    # Extract consecutive stop pairs with travel times
    arcs = Dict{Tuple{String,String,String}, Vector{Float64}}()  # (origin, dest, line) -> travel times

    prev_row = nothing
    for row in eachrow(stop_times)
        if !isnothing(prev_row) && prev_row.trip_id == row.trip_id
            # Same trip, consecutive stops
            origin = prev_row.stop_id
            dest = row.stop_id
            route = row.route_id

            # Calculate travel time in minutes
            arr_time = parse_time(row.arrival_time)
            dep_time = parse_time(prev_row.departure_time)
            travel_time = (arr_time - dep_time) / 60.0

            if travel_time > 0 && travel_time < 30  # Filter unreasonable times
                key = (string(origin), string(dest), string(route))
                if !haskey(arcs, key)
                    arcs[key] = Float64[]
                end
                push!(arcs[key], travel_time)
            end
        end
        prev_row = row
    end

    println("Found $(length(arcs)) unique station-to-station arcs")

    # Aggregate to median travel time per arc
    arc_data = DataFrame(
        origin = String[],
        destination = String[],
        capacity = Int[],
        traveltime = Int[],
        category = String[]
    )

    # NYC subway car capacity ~150 passengers, 8-10 cars per train
    # Peak frequency varies by line: express ~12 trains/hr, local ~8 trains/hr
    # We'll use a simplified capacity model
    line_capacities = Dict(
        # Express lines (higher frequency)
        "1" => 300, "2" => 300, "3" => 300,
        "4" => 300, "5" => 300, "6" => 300,
        "7" => 300,
        "A" => 300, "C" => 250, "E" => 300,
        "B" => 250, "D" => 300, "F" => 300, "M" => 250,
        "G" => 200,
        "J" => 200, "Z" => 200,
        "L" => 250,
        "N" => 300, "Q" => 300, "R" => 250, "W" => 250,
        "S" => 150,  # Shuttles
        "SI" => 150  # Staten Island Railway
    )
    default_capacity = 200

    for ((origin, dest, route), times) in arcs
        # Use parent stop ID (remove direction suffix like N/S)
        origin_base = length(origin) > 1 && origin[end] in ['N', 'S'] ? origin[1:end-1] : origin
        dest_base = length(dest) > 1 && dest[end] in ['N', 'S'] ? dest[1:end-1] : dest

        origin_name = get(stop_names, origin, get(stop_names, origin_base, origin))
        dest_name = get(stop_names, dest, get(stop_names, dest_base, dest))

        line_name = get(route_names, route, route)
        capacity = get(line_capacities, line_name, default_capacity)

        # Format station names like Doha: Metro_StationName
        origin_fmt = "Metro_" * replace(origin_name, " " => "", "/" => "_", "-" => "", "'" => "")
        dest_fmt = "Metro_" * replace(dest_name, " " => "", "/" => "_", "-" => "", "'" => "")

        push!(arc_data, (
            origin = origin_fmt,
            destination = dest_fmt,
            capacity = capacity,
            traveltime = max(1, round(Int, median(times))),
            category = line_name * " Line"
        ))
    end

    # Remove duplicates (same origin/dest/line)
    unique!(arc_data, [:origin, :destination, :category])

    println("Adding transfer connections...")

    # Find stations served by multiple lines (transfer points)
    # Group stops by parent station (stop_id without direction suffix)
    stop_routes = Dict{String, Set{String}}()
    for row in eachrow(stop_times)
        stop_base = string(row.stop_id)
        stop_base = length(stop_base) > 1 && stop_base[end] in ['N', 'S'] ? stop_base[1:end-1] : stop_base

        if !haskey(stop_routes, stop_base)
            stop_routes[stop_base] = Set{String}()
        end
        push!(stop_routes[stop_base], string(row.route_id))
    end

    # Add transfer arcs (bidirectional) at multi-line stations
    transfer_time = 3  # minutes
    transfer_capacity = 500  # pedestrian capacity

    for (stop_base, routes_set) in stop_routes
        if length(routes_set) > 1
            stop_name = get(stop_names, stop_base, get(stop_names, stop_base * "N", stop_base))
            stop_fmt = "Metro_" * replace(stop_name, " " => "", "/" => "_", "-" => "", "'" => "")

            routes_list = collect(routes_set)
            for i in 1:length(routes_list)
                for j in (i+1):length(routes_list)
                    line_i = get(route_names, routes_list[i], routes_list[i])
                    line_j = get(route_names, routes_list[j], routes_list[j])

                    origin_platform = stop_fmt * "_" * line_i
                    dest_platform = stop_fmt * "_" * line_j

                    # Bidirectional transfer
                    push!(arc_data, (
                        origin = origin_platform,
                        destination = dest_platform,
                        capacity = transfer_capacity,
                        traveltime = transfer_time,
                        category = "Pedestrian Path"
                    ))
                    push!(arc_data, (
                        origin = dest_platform,
                        destination = origin_platform,
                        capacity = transfer_capacity,
                        traveltime = transfer_time,
                        category = "Pedestrian Path"
                    ))
                end
            end
        end
    end

    # Sort by category, then origin
    sort!(arc_data, [:category, :origin])

    # Save
    mkpath(dirname(output_file))
    CSV.write(output_file, arc_data)

    println("\nDone! Saved $(nrow(arc_data)) arcs to $output_file")
    println("  - $(count(r -> r.category != "Pedestrian Path", eachrow(arc_data))) train arcs")
    println("  - $(count(r -> r.category == "Pedestrian Path", eachrow(arc_data))) transfer arcs")

    return arc_data
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("NYC Metro Arcs Builder")
    println("=" ^ 50)
    println()

    if !isfile(joinpath(DATA_DIR, "stops.txt"))
        download_gtfs()
    end

    build_metroarcs()
end

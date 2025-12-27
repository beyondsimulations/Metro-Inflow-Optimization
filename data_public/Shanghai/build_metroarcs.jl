"""
Shanghai Metro Arcs Builder

Generates metroarcs_shanghai.csv from station info and OD flow data.
Follows methodology from the paper:
- Travel time: 80 km/h + 20s dwell, rounded to 1-min steps (min 1 min)
- Capacity: Derived from peak observed hourly OD flow / 60 (passengers per minute)
- Consistency: All arcs on same line have identical capacity (line maximum)

Usage:
    cd data_public/Shanghai
    julia --project=../../metroflow build_metroarcs.jl

Input:
    - stationInfo.csv (302 stations with coordinates and neighbor relationships)
    - metroData_ODFlow.csv (11GB OD flow data for capacity estimation)

Output:
    - metroarcs_shanghai.csv (origin, destination, capacity, traveltime, category)
"""

using Pkg
Pkg.activate("../../metroflow")

using CSV
using DataFrames
using Statistics
using Graphs
using ProgressMeter
using Dates
using Printf

# CONSTANTS

const TRAIN_SPEED_KMH = 80.0    # km/h
const DWELL_TIME_MIN = 2/6      # 20 seconds

"""
Haversine formula to calculate distance between two points in kilometers.
"""
function haversine(lat1, lon1, lat2, lon2)
    R = 6371.0  # Earth radius in km
    φ1 = deg2rad(lat1)
    φ2 = deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)

    a = sin(Δφ/2)^2 + cos(φ1) * cos(φ2) * sin(Δλ/2)^2
    c = 2 * atan(sqrt(a), sqrt(1-a))

    return R * c
end

"""
Estimate travel time in minutes from distance.
Uses 80 km/h train speed + 30s dwell time, rounded to 1-minute steps.
Minimum 1 minute.
"""
function estimate_travel_time(distance_km)::Int
    travel = (distance_km / TRAIN_SPEED_KMH) * 60  # minutes
    return max(1, round(Int, travel + DWELL_TIME_MIN))
end

"""
Format station name to match convention: Metro_StationName
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
Parse the neighbor string from stationInfo.csv.
Format: "[312, 314, 2012, 2038]"
"""
function parse_neighbors(neighbor_str::AbstractString)
    cleaned = strip(neighbor_str, ['[', ']', ' '])
    if isempty(cleaned)
        return Int[]
    end
    parts = split(cleaned, ",")
    return [parse(Int, strip(p)) for p in parts if !isempty(strip(p))]
end

"""
Load station line assignments from station_lines.csv.
Returns: Dict{String, Vector{String}} mapping station name → [lines]
"""
function load_station_lines(station_lines_file::String="station_lines.csv")
    println("Loading station line assignments from $station_lines_file...")
    df = CSV.read(station_lines_file, DataFrame)

    station_lines = Dict{String, Vector{String}}()
    for row in eachrow(df)
        name = row[Symbol("Station Name")]
        lines_str = row[Symbol("Metro Line(s)")]
        # Parse "Line 1, Line 3, Line 4" → ["Line 1", "Line 3", "Line 4"]
        lines = [strip(String(l)) for l in split(lines_str, ",")]
        station_lines[name] = lines
    end

    println("  Loaded $(length(station_lines)) station line assignments")
    return station_lines
end

"""
Normalize station name for matching.
"""
function normalize_for_match(name::String)
    n = lowercase(name)
    n = replace(n, " " => "")
    n = replace(n, "'" => "")
    n = replace(n, "-" => "")
    n = replace(n, "·" => "")  # Handle special middle dot
    # Don't remove "road" or "lu" - too aggressive for substring matching
    return n
end

"""
Build mapping from station ID to line(s) by matching names.
Uses 3-level matching: exact → normalized → substring/fuzzy.
"""
function build_station_to_lines(station_names::Dict{Int, String},
                                 station_lines::Dict{String, Vector{String}})
    station_to_lines = Dict{Int, Vector{String}}()
    unmatched = String[]
    matched_fuzzy = Dict{String, String}()  # Track fuzzy matches for logging

    # Build normalized lookup
    normalized_lookup = Dict{String, Tuple{String, Vector{String}}}()
    for (ref_name, lines) in station_lines
        normalized_lookup[normalize_for_match(ref_name)] = (ref_name, lines)
    end

    for (id, name) in station_names
        # Level 1: Exact match
        if haskey(station_lines, name)
            station_to_lines[id] = station_lines[name]
            continue
        end

        # Level 2: Normalized match
        normalized = normalize_for_match(name)
        if haskey(normalized_lookup, normalized)
            station_to_lines[id] = normalized_lookup[normalized][2]
            continue
        end

        # Level 3: Substring match (e.g., "Xintiandi" in "Site of...Xintiandi")
        # Only match if substring is at least 6 characters to avoid false positives
        found = false
        if length(normalized) >= 6
            for (ref_name, lines) in station_lines
                ref_norm = normalize_for_match(ref_name)
                if length(ref_norm) >= 6 && (occursin(normalized, ref_norm) || occursin(ref_norm, normalized))
                    station_to_lines[id] = lines
                    matched_fuzzy[name] = ref_name
                    found = true
                    break
                end
            end
        end

        if !found
            push!(unmatched, name)
        end
    end

    # Log fuzzy matches
    if !isempty(matched_fuzzy)
        println("  Fuzzy matches:")
        for (orig, matched) in matched_fuzzy
            println("    '$orig' → '$matched'")
        end
    end

    if !isempty(unmatched)
        println("  Warning: $(length(unmatched)) stations not matched (will use 'Metro'):")
        for n in unmatched[1:min(10, length(unmatched))]
            println("    - $n")
        end
    end

    return station_to_lines
end

"""
Get line name for an arc between two stations.
Handles multi-line stations by finding the common line.
"""
function get_arc_line(station1_id::Int, station2_id::Int,
                       station_to_lines::Dict{Int, Vector{String}})
    lines1 = get(station_to_lines, station1_id, String[])
    lines2 = get(station_to_lines, station2_id, String[])

    if isempty(lines1) && isempty(lines2)
        return "Metro"
    elseif isempty(lines1)
        return lines2[1]  # Use first line of destination
    elseif isempty(lines2)
        return lines1[1]  # Use first line of origin
    else
        # Find common line
        common = intersect(lines1, lines2)
        if !isempty(common)
            return common[1]  # Use first common line
        else
            return lines1[1]  # At transfer stations, use origin's primary line
        end
    end
end

# ============================================================================
# EXPANDED NETWORK MODEL
# Transfer stations have multiple nodes: central + one per line
# ============================================================================

const PEDESTRIAN_CAPACITY = 9999  # Effectively unlimited
const TRANSFER_TIME = 1           # 1 minute for transfers

"""
Extract line number from line name (e.g., "Line 1" → 1, "Line 16" → 16).
"""
function get_line_number(line::String)::Int
    m = match(r"(\d+)", line)
    return isnothing(m) ? 0 : parse(Int, m.captures[1])
end

"""
Format line suffix with zero-padding (e.g., 1 → "L01", 16 → "L16").
"""
function format_line_suffix(line::String)::String
    num = get_line_number(line)
    return @sprintf("L%02d", num)
end

"""
Check if a station is a transfer station (serves 2+ lines).
"""
function is_transfer_station(station_id::Int, station_to_lines::Dict{Int, Vector{String}})::Bool
    lines = get(station_to_lines, station_id, String[])
    return length(lines) > 1
end

"""
Detect triangles in the network and return edges to skip (longest edge in each triangle).

When three stations A, B, C on the same line form a triangle (all pairs are neighbors),
the longest edge is likely an incorrect "shortcut" that should be skipped.

Returns: Set of (station_id1, station_id2, line) tuples to skip
"""
function detect_triangle_shortcuts(station_neighbors::Dict{Int, Vector{Int}},
                                    station_to_lines::Dict{Int, Vector{String}},
                                    station_coords::Dict{Int, Tuple{Float64, Float64}},
                                    station_names::Dict{Int, String})
    shortcuts_to_skip = Set{Tuple{Int, Int, String}}()

    # Build per-line adjacency: line -> set of (station_id, station_id) edges
    line_edges = Dict{String, Set{Tuple{Int, Int}}}()
    line_stations = Dict{String, Set{Int}}()

    for (id, neighbors) in station_neighbors
        lines = get(station_to_lines, id, String[])
        for neighbor_id in neighbors
            haskey(station_to_lines, neighbor_id) || continue
            neighbor_lines = get(station_to_lines, neighbor_id, String[])
            common_lines = intersect(lines, neighbor_lines)

            for line in common_lines
                if !haskey(line_edges, line)
                    line_edges[line] = Set{Tuple{Int, Int}}()
                    line_stations[line] = Set{Int}()
                end
                # Store edges with smaller ID first for consistency
                edge = id < neighbor_id ? (id, neighbor_id) : (neighbor_id, id)
                push!(line_edges[line], edge)
                push!(line_stations[line], id)
                push!(line_stations[line], neighbor_id)
            end
        end
    end

    # For each line, find triangles
    n_triangles = 0
    for (line, edges) in line_edges
        stations = collect(line_stations[line])
        n = length(stations)

        for i in 1:n
            for j in (i+1):n
                for k in (j+1):n
                    a, b, c = stations[i], stations[j], stations[k]

                    # Check if all three pairs are connected
                    ab = (min(a,b), max(a,b)) in edges
                    bc = (min(b,c), max(b,c)) in edges
                    ac = (min(a,c), max(a,c)) in edges

                    if ab && bc && ac
                        # Triangle found! Skip the longest edge
                        n_triangles += 1

                        # Calculate distances
                        lat_a, lon_a = station_coords[a]
                        lat_b, lon_b = station_coords[b]
                        lat_c, lon_c = station_coords[c]

                        dist_ab = haversine(lat_a, lon_a, lat_b, lon_b)
                        dist_bc = haversine(lat_b, lon_b, lat_c, lon_c)
                        dist_ac = haversine(lat_a, lon_a, lat_c, lon_c)

                        # Find longest edge
                        max_dist = max(dist_ab, dist_bc, dist_ac)
                        if dist_ab == max_dist
                            push!(shortcuts_to_skip, (a, b, line))
                            push!(shortcuts_to_skip, (b, a, line))
                            name_a = get(station_names, a, "?")
                            name_b = get(station_names, b, "?")
                            println("  Triangle on $line: skipping $name_a ↔ $name_b (longest edge: $(round(dist_ab, digits=2)) km)")
                        elseif dist_bc == max_dist
                            push!(shortcuts_to_skip, (b, c, line))
                            push!(shortcuts_to_skip, (c, b, line))
                            name_b = get(station_names, b, "?")
                            name_c = get(station_names, c, "?")
                            println("  Triangle on $line: skipping $name_b ↔ $name_c (longest edge: $(round(dist_bc, digits=2)) km)")
                        else
                            push!(shortcuts_to_skip, (a, c, line))
                            push!(shortcuts_to_skip, (c, a, line))
                            name_a = get(station_names, a, "?")
                            name_c = get(station_names, c, "?")
                            println("  Triangle on $line: skipping $name_a ↔ $name_c (longest edge: $(round(dist_ac, digits=2)) km)")
                        end
                    end
                end
            end
        end
    end

    println("  Found $n_triangles triangles, $(length(shortcuts_to_skip) ÷ 2) edges to skip")
    return shortcuts_to_skip
end

"""
Build expanded stations DataFrame with coordinates.
Transfer stations get central node + line-specific nodes with offset positions.

Returns: DataFrame with columns:
  - node_id: e.g., "Metro_Xujiahui" or "Metro_Xujiahui_L01"
  - station_id: numeric ID from stationInfo.csv (for debugging)
  - station_name: human-readable name
  - node_type: "regular", "central", or "line"
  - line: "" for central/regular, "Line X" for line nodes
  - lon, lat: coordinates (with offsets for line nodes)
"""
function build_expanded_stations(station_coords::Dict{Int, Tuple{Float64, Float64}},
                                  station_names::Dict{Int, String},
                                  station_to_lines::Dict{Int, Vector{String}})
    println("Building expanded station list...")

    expanded = DataFrame(
        node_id = String[],
        station_id = Int[],
        station_name = String[],
        node_type = String[],
        line = String[],
        lon = Float64[],
        lat = Float64[]
    )

    # Offset in degrees (~50m at Shanghai's latitude)
    # 1 degree ≈ 111km, so 0.0005 ≈ 55m
    offset = 0.0005

    n_regular = 0
    n_transfer = 0
    n_line_nodes = 0

    for (id, name) in station_names
        lat, lon = station_coords[id]  # Note: station_coords stores (lat, lon)
        lines = get(station_to_lines, id, String[])
        metro_name = format_station_name(name)

        if length(lines) <= 1
            # Regular station: single node (store the line if it has one)
            line_str = isempty(lines) ? "" : lines[1]
            push!(expanded, (metro_name, id, name, "regular", line_str, lon, lat))
            n_regular += 1
        else
            # Transfer station: central + line nodes
            push!(expanded, (metro_name, id, name, "central", "", lon, lat))
            n_transfer += 1

            # Create line nodes with offset positions distributed around center
            n_lines = length(lines)
            for (i, line) in enumerate(lines)
                angle = 2π * (i - 1) / n_lines  # Distribute evenly around center
                line_suffix = format_line_suffix(line)
                node_id = "$(metro_name)_$(line_suffix)"
                offset_lon = lon + offset * cos(angle)
                offset_lat = lat + offset * sin(angle)
                push!(expanded, (node_id, id, name, "line", line, offset_lon, offset_lat))
                n_line_nodes += 1
            end
        end
    end

    println("  Regular stations: $n_regular")
    println("  Transfer stations: $n_transfer")
    println("  Line nodes created: $n_line_nodes")
    println("  Total nodes: $(nrow(expanded))")

    return expanded
end

"""
Build expanded network arcs with proper parallel track handling.

For shared track segments (e.g., Lines 3 & 4), creates separate arcs for each line.
Adds pedestrian transfer arcs at transfer stations.
Skips edges identified as triangle shortcuts.

Returns: DataFrame with columns: origin, destination, capacity, traveltime, category
"""
function build_expanded_arcs(station_coords::Dict{Int, Tuple{Float64, Float64}},
                              station_names::Dict{Int, String},
                              station_neighbors::Dict{Int, Vector{Int}},
                              station_to_lines::Dict{Int, Vector{String}},
                              line_capacities::Dict{String, Int},
                              shortcuts_to_skip::Set{Tuple{Int, Int, String}}=Set{Tuple{Int, Int, String}}())
    println("Building expanded network arcs...")

    arcs = DataFrame(
        origin = String[],
        destination = String[],
        capacity = Int[],
        traveltime = Int[],
        category = String[]
    )

    n_track_arcs = 0
    n_entry_exit_arcs = 0
    n_transfer_arcs = 0

    for (id, neighbors) in station_neighbors
        name = station_names[id]
        metro_name = format_station_name(name)
        lines = get(station_to_lines, id, String[])
        is_transfer = length(lines) > 1

        # Process track arcs to neighbors
        for neighbor_id in neighbors
            if !haskey(station_names, neighbor_id)
                continue
            end

            neighbor_name = station_names[neighbor_id]
            neighbor_metro = format_station_name(neighbor_name)
            neighbor_lines = get(station_to_lines, neighbor_id, String[])
            neighbor_is_transfer = length(neighbor_lines) > 1

            # Find common lines between stations
            common_lines = intersect(lines, neighbor_lines)

            # Calculate travel time from distance
            lat1, lon1 = station_coords[id]
            lat2, lon2 = station_coords[neighbor_id]
            distance = haversine(lat1, lon1, lat2, lon2)
            travel_time = estimate_travel_time(distance)

            # Create arc for EACH common line (this handles parallel tracks!)
            for line in common_lines
                # Skip if this edge was identified as a triangle shortcut
                if (id, neighbor_id, line) in shortcuts_to_skip
                    continue
                end

                line_suffix = format_line_suffix(line)
                capacity = get(line_capacities, line, 500)

                # Determine origin/destination node names
                origin_node = is_transfer ? "$(metro_name)_$(line_suffix)" : metro_name
                dest_node = neighbor_is_transfer ? "$(neighbor_metro)_$(line_suffix)" : neighbor_metro

                push!(arcs, (origin_node, dest_node, capacity, travel_time, line))
                n_track_arcs += 1
            end
        end

        # Add pedestrian arcs for transfer stations
        if is_transfer
            # Central ↔ Line nodes (entry/exit)
            for line in lines
                line_suffix = format_line_suffix(line)
                line_node = "$(metro_name)_$(line_suffix)"

                # Bidirectional entry/exit
                push!(arcs, (metro_name, line_node, PEDESTRIAN_CAPACITY, TRANSFER_TIME, "Pedestrian"))
                push!(arcs, (line_node, metro_name, PEDESTRIAN_CAPACITY, TRANSFER_TIME, "Pedestrian"))
                n_entry_exit_arcs += 2
            end

            # Line ↔ Line nodes (cross-platform transfers)
            for (i, line1) in enumerate(lines)
                for line2 in lines[i+1:end]
                    ln1 = format_line_suffix(line1)
                    ln2 = format_line_suffix(line2)
                    node1 = "$(metro_name)_$(ln1)"
                    node2 = "$(metro_name)_$(ln2)"

                    # Bidirectional transfer
                    push!(arcs, (node1, node2, PEDESTRIAN_CAPACITY, TRANSFER_TIME, "Pedestrian"))
                    push!(arcs, (node2, node1, PEDESTRIAN_CAPACITY, TRANSFER_TIME, "Pedestrian"))
                    n_transfer_arcs += 2
                end
            end
        end
    end

    println("  Track arcs: $n_track_arcs")
    println("  Entry/exit arcs: $n_entry_exit_arcs")
    println("  Transfer arcs: $n_transfer_arcs")
    println("  Total arcs: $(nrow(arcs))")

    return arcs
end

"""
Load station data and build lookup structures.
"""
function load_stations(station_file::String)
    println("Loading station data from $station_file...")
    stations = CSV.read(station_file, DataFrame)

    station_coords = Dict{Int, Tuple{Float64, Float64}}()
    station_names = Dict{Int, String}()
    station_neighbors = Dict{Int, Vector{Int}}()

    for row in eachrow(stations)
        sid = row.stationID
        station_coords[sid] = (row.lat, row.lon)
        station_names[sid] = row.name
        station_neighbors[sid] = parse_neighbors(string(row.neighbour))
    end

    println("  Loaded $(length(station_names)) stations")
    return station_coords, station_names, station_neighbors
end

"""
Build directed arcs from neighbor relationships.
Returns: arcs dict with (origin_id, dest_id) => (distance_km, travel_time, line)
"""
function build_network_arcs(station_coords, station_names, station_neighbors, station_to_lines)
    println("Building network arcs...")

    arcs = Dict{Tuple{Int, Int}, Tuple{Float64, Int, String}}()

    for (origin_id, neighbors) in station_neighbors
        origin_lat, origin_lon = station_coords[origin_id]

        for dest_id in neighbors
            if !haskey(station_coords, dest_id)
                @warn "Neighbor station $dest_id not found for station $origin_id"
                continue
            end

            dest_lat, dest_lon = station_coords[dest_id]
            distance = haversine(origin_lat, origin_lon, dest_lat, dest_lon)
            travel_time = estimate_travel_time(distance)
            line = get_arc_line(origin_id, dest_id, station_to_lines)

            arcs[(origin_id, dest_id)] = (distance, travel_time, line)
        end
    end

    println("  Built $(length(arcs)) directed arcs")
    return arcs
end

"""
Build a Graphs.jl graph for shortest path computation.
Returns: graph, node_id_to_station_id mapping, station_id_to_node_id mapping
"""
function build_graph(station_names, arcs)
    println("Building graph for shortest path computation...")

    # Create mappings
    station_ids = sort(collect(keys(station_names)))
    station_to_node = Dict(sid => i for (i, sid) in enumerate(station_ids))
    node_to_station = Dict(i => sid for (i, sid) in enumerate(station_ids))

    # Build graph
    n = length(station_ids)
    g = SimpleDiGraph(n)
    weights = zeros(Int, n, n)

    for ((origin_id, dest_id), (_, travel_time, _)) in arcs
        if haskey(station_to_node, origin_id) && haskey(station_to_node, dest_id)
            o_node = station_to_node[origin_id]
            d_node = station_to_node[dest_id]
            add_edge!(g, o_node, d_node)
            weights[o_node, d_node] = travel_time
        end
    end

    println("  Graph: $(nv(g)) nodes, $(ne(g)) edges")
    return g, node_to_station, station_to_node, weights
end


"""
Precompute shortest paths between all station pairs.
Returns: paths dict with (origin_station_id, dest_station_id) => [arc1, arc2, ...]
         where each arc is (from_station_id, to_station_id)
"""
function precompute_paths(g, node_to_station, station_to_node, weights)
    println("Precomputing shortest paths...")

    n = nv(g)
    paths = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()

    # Compute all-pairs shortest paths using Dijkstra from each source
    @showprogress "Computing paths: " for source_node in 1:n
        source_station = node_to_station[source_node]

        # Dijkstra from this source
        state = dijkstra_shortest_paths(g, source_node, weights)

        for dest_node in 1:n
            if dest_node == source_node
                continue
            end
            dest_station = node_to_station[dest_node]

            # Reconstruct path
            if state.dists[dest_node] < typemax(Int)
                path_arcs = Tuple{Int, Int}[]
                current = dest_node
                while state.parents[current] != 0
                    parent = state.parents[current]
                    push!(path_arcs, (node_to_station[parent], node_to_station[current]))
                    current = parent
                end
                reverse!(path_arcs)
                paths[(source_station, dest_station)] = path_arcs
            end
        end
    end

    println("  Computed $(length(paths)) paths")
    return paths
end

"""
Parse date integer (YYYYMMDD) to get hour from startTime (HHMMSS).
"""
function get_hour(start_time_int::Int)::Int
    return div(start_time_int, 10000)
end

"""
Compute empirical arc capacities from OD flow data.
Aggregates flow to hourly level, finds peak, divides by 60 for per-minute capacity.
Ensures consistent capacity per line using line maximum.
"""
function compute_empirical_capacity(od_file::String, paths, arcs, station_names)
    println("\nComputing empirical capacity from OD data...")
    println("  This may take 20-30 minutes for the 11GB file...")

    # Arc hourly flows: arc => Dict(date_hour => flow)
    arc_hourly_flows = Dict{Tuple{Int, Int}, Dict{Tuple{Int, Int}, Float64}}()
    for arc in keys(arcs)
        arc_hourly_flows[arc] = Dict{Tuple{Int, Int}, Float64}()
    end

    # Build station ID lookup from names
    name_to_id = Dict(format_station_name(name) => id for (id, name) in station_names)

    # Process OD file
    println("  Loading OD data...")
    df = CSV.read(od_file, DataFrame; normalizenames=true,
        types=Dict(:date => Int, :timeslot => Int, :startTime => Int, :endTime => Int,
                   :originStation => Int, :destinationStation => Int, :Flow => Int))

    println("  Loaded $(nrow(df)) OD records")
    println("  Routing flows through network...")

    missing_paths = 0
    routed_flows = 0

    @showprogress "Processing OD flows: " for row in eachrow(df)
        origin_id = row.originStation
        dest_id = row.destinationStation
        flow = row.Flow
        date = row.date
        hour = get_hour(row.startTime)

        # Skip if same station
        if origin_id == dest_id
            continue
        end

        # Get path for this OD pair
        path_key = (origin_id, dest_id)
        if !haskey(paths, path_key)
            missing_paths += 1
            continue
        end

        path = paths[path_key]
        routed_flows += 1

        # Add flow to each arc in path
        date_hour = (date, hour)
        for arc in path
            if haskey(arc_hourly_flows, arc)
                if !haskey(arc_hourly_flows[arc], date_hour)
                    arc_hourly_flows[arc][date_hour] = 0.0
                end
                arc_hourly_flows[arc][date_hour] += flow
            end
        end
    end

    println("  Routed $routed_flows flows, $missing_paths had no path")

    # Find peak hourly flow for each arc
    println("  Computing peak flows per arc...")
    arc_peak_flows = Dict{Tuple{Int, Int}, Float64}()
    for (arc, hourly_flows) in arc_hourly_flows
        if isempty(hourly_flows)
            arc_peak_flows[arc] = 0.0
        else
            arc_peak_flows[arc] = maximum(values(hourly_flows))
        end
    end

    # Aggregate by line and find maximum
    println("  Computing per-line capacity (maximum across arcs)...")
    line_peak_flows = Dict{String, Float64}()
    for (arc, peak) in arc_peak_flows
        (origin_id, dest_id) = arc
        if haskey(arcs, arc)
            (_, _, line) = arcs[arc]
            if !haskey(line_peak_flows, line)
                line_peak_flows[line] = 0.0
            end
            line_peak_flows[line] = max(line_peak_flows[line], peak)
        end
    end

    # Convert to per-minute capacity (hourly peak / 60)
    line_capacities = Dict{String, Int}()
    for (line, peak) in line_peak_flows
        capacity = max(1, round(Int, peak / 60))
        line_capacities[line] = capacity
        println("    $line: peak=$(@sprintf("%.0f", peak))/hour → $capacity/min")
    end

    return line_capacities
end

function build_metroarcs(; output_file="metroarcs_shanghai.csv",
                          stations_file="stations_shanghai.csv",
                          od_file="metroData_ODFlow.csv",
                          station_lines_file="station_lines_2017.csv")
    println("=" ^ 60)
    println("Shanghai Metro Arcs Builder (Expanded Network Model)")
    println("=" ^ 60)
    println()
    println("Parameters:")
    println("  Train speed: $(TRAIN_SPEED_KMH) km/h")
    println("  Dwell time: $(DWELL_TIME_MIN * 60) seconds")
    println("  Transfer model: EXPANDED (central + line nodes)")
    println("  Transfer time: $(TRANSFER_TIME) minute")
    println("  Pedestrian capacity: $(PEDESTRIAN_CAPACITY) pax/min")
    println()

    # Step 1: Load station data
    station_coords, station_names, station_neighbors = load_stations("stationInfo.csv")

    # Step 2: Load line assignments from CSV and match to stations
    station_lines = load_station_lines(station_lines_file)
    station_to_lines = build_station_to_lines(station_names, station_lines)

    # Step 3: Build simple network arcs (for capacity calculation via shortest paths)
    # We still need this for routing OD flows through the network
    simple_arcs = build_network_arcs(station_coords, station_names, station_neighbors, station_to_lines)

    # Step 4: Build graph for shortest paths
    g, node_to_station, station_to_node, weights = build_graph(station_names, simple_arcs)

    # Step 5: Precompute shortest paths
    paths = precompute_paths(g, node_to_station, station_to_node, weights)

    # Step 6: Compute empirical capacity from OD data
    if isfile(od_file)
        line_capacities = compute_empirical_capacity(od_file, paths, simple_arcs, station_names)
    else
        println("\nWarning: OD file not found ($od_file)")
        println("Using default capacity of 500 passengers/min for all lines")
        line_capacities = Dict{String, Int}()
    end

    # Step 7: Build EXPANDED network
    println("\n" * "-" ^ 60)
    println("Building expanded network model...")
    println("-" ^ 60)

    # Detect triangle shortcuts to skip (prevents incorrect direct connections)
    println("\nDetecting triangle shortcuts...")
    shortcuts_to_skip = detect_triangle_shortcuts(station_neighbors, station_to_lines,
                                                   station_coords, station_names)

    # Build expanded stations with offset coordinates for line nodes
    expanded_stations = build_expanded_stations(station_coords, station_names, station_to_lines)

    # Build expanded arcs with parallel tracks and pedestrian transfers
    arc_data = build_expanded_arcs(station_coords, station_names, station_neighbors,
                                    station_to_lines, line_capacities, shortcuts_to_skip)

    # Sort and deduplicate arcs
    sort!(arc_data, [:category, :origin])
    unique!(arc_data, [:origin, :destination, :category])

    # Save stations file
    CSV.write(stations_file, expanded_stations)
    println("\nSaved: $stations_file")

    # Save arcs file
    CSV.write(output_file, arc_data)
    println("Saved: $output_file")

    # Summary statistics
    println("\n" * "=" ^ 60)
    println("SUMMARY")
    println("=" ^ 60)

    # Station counts
    n_regular = count(x -> x == "regular", expanded_stations.node_type)
    n_central = count(x -> x == "central", expanded_stations.node_type)
    n_line = count(x -> x == "line", expanded_stations.node_type)

    println("\nStations:")
    println("  Physical stations: $(n_regular + n_central)")
    println("    - Regular (1 line): $n_regular")
    println("    - Transfer (2+ lines): $n_central")
    println("  Line nodes: $n_line")
    println("  Total nodes: $(nrow(expanded_stations))")

    # Arc counts
    n_track = count(x -> x != "Pedestrian", arc_data.category)
    n_pedestrian = count(x -> x == "Pedestrian", arc_data.category)

    println("\nArcs:")
    println("  Track arcs: $n_track")
    println("  Pedestrian arcs: $n_pedestrian")
    println("  Total arcs: $(nrow(arc_data))")

    # Travel time statistics (track arcs only)
    track_arcs = filter(row -> row.category != "Pedestrian", arc_data)
    println("\nTravel time statistics (track arcs only):")
    println("  - Min: $(minimum(track_arcs.traveltime)) min")
    println("  - Max: $(maximum(track_arcs.traveltime)) min")
    println("  - Mean: $(round(mean(track_arcs.traveltime), digits=1)) min")

    # Capacity statistics (track arcs only, excluding pedestrian)
    println("\nCapacity statistics (track arcs only):")
    println("  - Min: $(minimum(track_arcs.capacity)) pax/min")
    println("  - Max: $(maximum(track_arcs.capacity)) pax/min")
    println("  - Mean: $(round(mean(track_arcs.capacity), digits=0)) pax/min")

    # Per-line capacities
    println("\nPer-line capacities:")
    for line in sort(collect(keys(line_capacities)))
        println("  - $line: $(line_capacities[line]) pax/min")
    end

    # Arc counts per category
    println("\nArcs per category:")
    category_counts = combine(groupby(arc_data, :category), nrow => :count)
    sort!(category_counts, :category)
    for row in eachrow(category_counts)
        println("  - $(row.category): $(row.count) arcs")
    end

    println("\nDone!")

    return arc_data, expanded_stations
end

if abspath(PROGRAM_FILE) == @__FILE__
    build_metroarcs()
end

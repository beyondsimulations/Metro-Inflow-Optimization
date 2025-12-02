using CSV
using DataFrames
using Statistics

const DATA_DIR = "data_public/Shanghai"

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
Estimate travel time in minutes from distance (returns integer).
Assumes average metro speed of 35 km/h (including acceleration, stops).
Adds 0.5 min base time for station dwell. Minimum 1 minute.
"""
function estimate_travel_time(distance_km)::Int
    avg_speed_kmh = 35.0
    base_dwell = 0.5  # minutes
    travel = (distance_km / avg_speed_kmh) * 60  # convert to minutes
    return max(1, round(Int, travel + base_dwell))
end

"""
Parse the neighbor string from stationInfo.csv.
Format: "[312, 314, 2012, 2038]"
"""
function parse_neighbors(neighbor_str::AbstractString)
    # Remove brackets and split
    cleaned = strip(neighbor_str, ['[', ']', ' '])
    if isempty(cleaned)
        return Int[]
    end
    parts = split(cleaned, ",")
    return [parse(Int, strip(p)) for p in parts if !isempty(strip(p))]
end

"""
Shanghai Metro line information based on station groupings.
Maps station IDs to their primary line.
"""
function get_line_assignments()
    # Shanghai Metro line assignments based on known station IDs
    # This is a comprehensive mapping of station IDs to lines
    # Lines 1-17 existed in 2017 Shanghai metro

    line_map = Dict{Int, String}()

    # Line 1 (Red): Xinzhuang - Fujin Road
    line1 = [112, 113, 114, 119, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,
             312, 314, 2011, 2012, 2025, 2027, 2029, 2038, 2048]
    for s in line1; line_map[s] = "Line 1"; end

    # Line 2 (Green): East Xujing - Pudong Airport / Guanglan Road
    line2 = [234, 237, 238, 239, 240, 247, 248, 250, 251, 253, 254, 255, 256, 257, 258, 259,
             260, 261, 262, 263, 2010, 2013, 2014, 2019, 2035, 2043, 2046]
    for s in line2; line_map[s] = "Line 2"; end

    # Line 3 (Yellow): Shanghai South Railway Station - North Jiangyang Road
    line3 = [325, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339,
             837, 838, 2005, 2021, 2027, 2030, 2034]
    for s in line3; line_map[s] = "Line 3"; end

    # Line 4 (Purple loop): Inner ring / Outer ring
    line4 = [413, 415, 416, 420, 421, 423, 426, 2002, 2003, 2005, 2006, 2007, 2010,
             2015, 2017, 2025, 2032, 2033, 2047]
    for s in line4; line_map[s] = "Line 4"; end

    # Line 5 (Magenta): Xinzhuang - Minhang Development Zone
    line5 = [502, 503, 505, 507, 508, 509, 510, 511, 512, 513, 2011]
    for s in line5; line_map[s] = "Line 5"; end

    # Line 6 (Pink): Gangcheng Road - Oriental Sports Center
    line6 = [622, 623, 624, 625, 626, 628, 629, 633, 634, 635, 636, 637, 638, 639, 640,
             642, 643, 644, 645, 646, 647, 648, 2004, 2010, 2015, 2023, 2040]
    for s in line6; line_map[s] = "Line 6"; end

    # Line 7 (Orange): Meilan Lake - Huamu Road
    line7 = [721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735,
             738, 744, 745, 747, 749, 750, 751, 753, 2016, 2020, 2023, 2024, 2035,
             2041, 2045, 2050, 2052, 2055]
    for s in line7; line_map[s] = "Line 7"; end

    # Line 8 (Blue): Shiguang Road - Shendu Highway
    line8 = [820, 821, 822, 823, 824, 825, 827, 828, 830, 834, 840, 842, 843, 844, 845,
             846, 847, 848, 849, 2002, 2005, 2008, 2016, 2018, 2035, 2040]
    for s in line8; line_map[s] = "Line 8"; end

    # Line 9 (Light Blue): Songjiang South Railway Station - Caolu
    line9 = [918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932,
             937, 940, 941, 943, 2006, 2025, 2031, 2033, 2036, 2044, 2051, 2010]
    for s in line9; line_map[s] = "Line 9"; end

    # Line 10 (Light Purple): Hongqiao Terminal 1 - Xinjiangwan City
    line10 = [1043, 1044, 1045, 1046, 1047, 1048, 1051, 1055, 1058, 1060, 1062, 1063,
              1064, 1065, 1066, 1067, 1068, 2008, 2013, 2019, 2026, 2037, 2044, 2047,
              2050, 2054, 2056]
    for s in line10; line_map[s] = "Line 10"; end

    # Line 11 (Dark Red): North Jiading/Disney - Luoshan Road
    line11 = [1114, 1115, 1116, 1117, 1118, 1119, 1120, 1131, 1132, 1133, 1134, 1135,
              1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1150, 1152, 1153, 1155,
              1156, 1157, 1159, 1161, 1162, 1163, 2009, 2012, 2025, 2043, 2054, 2055, 2057]
    for s in line11; line_map[s] = "Line 11"; end

    # Line 12 (Dark Green): Jinhai Road - Caoxi Road
    line12 = [1220, 1221, 1222, 1223, 1224, 1225, 1226, 1238, 1239, 1241, 1242, 1243,
              1244, 1245, 1246, 1248, 1249, 1250, 2004, 2022, 2026, 2032, 2038, 2048, 2050]
    for s in line12; line_map[s] = "Line 12"; end

    # Line 13 (Pink/Light Red): Jinyun Road - Shibo Avenue
    line13 = [1321, 1322, 1323, 1324, 1325, 1326, 1329, 1331, 1333, 1338, 1339,
              2009, 2020, 2028, 2031, 2039, 2041, 2055, 2056, 2057]
    for s in line13; line_map[s] = "Line 13"; end

    # Line 16 (Cyan): Longyang Road - Dishui Lake
    line16 = [1622, 1624, 1625, 1626, 1627, 1628, 1629, 1630, 1631, 1632, 1633, 2042, 2052]
    for s in line16; line_map[s] = "Line 16"; end

    # Line 17: Hongqiao Railway Station - Oriental Land (opened late 2017)
    line17 = [1018, 1019, 1020, 2013, 2014, 2043, 2044, 2051]
    for s in line17; line_map[s] = "Line 17"; end

    return line_map
end

"""
Get line name for an arc between two stations.
If both stations are on the same line, use that line.
Otherwise, try to find a common line or use the first station's line.
"""
function get_arc_line(station1_id::Int, station2_id::Int, line_map::Dict{Int, String})
    line1 = get(line_map, station1_id, nothing)
    line2 = get(line_map, station2_id, nothing)

    if isnothing(line1) && isnothing(line2)
        return "Metro"
    elseif isnothing(line1)
        return line2
    elseif isnothing(line2)
        return line1
    elseif line1 == line2
        return line1
    else
        # At transfer stations, prefer the first station's line
        return line1
    end
end

"""
Line capacity configuration for Shanghai Metro.
Based on train frequency and car capacity.
"""
function get_line_capacities()
    return Dict(
        "Line 1" => 300,    # High frequency trunk line
        "Line 2" => 300,    # High frequency trunk line
        "Line 3" => 250,    # Medium frequency
        "Line 4" => 250,    # Loop line
        "Line 5" => 200,    # Branch line
        "Line 6" => 200,    # Medium frequency
        "Line 7" => 280,    # High frequency
        "Line 8" => 280,    # High frequency
        "Line 9" => 280,    # High frequency
        "Line 10" => 250,   # Medium frequency
        "Line 11" => 250,   # Long line with branches
        "Line 12" => 250,   # Medium frequency
        "Line 13" => 250,   # Medium frequency
        "Line 16" => 200,   # Lower frequency suburban
        "Line 17" => 200,   # Newer suburban line
        "Metro" => 200      # Default
    )
end

"""
Format station name to match NYC convention: Metro_StationName
"""
function format_station_name(name::AbstractString)
    # Remove special characters and spaces
    formatted = replace(name, " " => "")
    formatted = replace(formatted, "'" => "")
    formatted = replace(formatted, "-" => "")
    formatted = replace(formatted, "/" => "_")
    formatted = replace(formatted, "(" => "")
    formatted = replace(formatted, ")" => "")
    return "Metro_" * formatted
end

"""
Build metro arcs from Shanghai station data.
"""
function build_metroarcs(; output_file="data_public/Shanghai/metroarcs_shanghai.csv")
    println("Loading Shanghai station data...")
    stations = CSV.read(joinpath(DATA_DIR, "stationInfo.csv"), DataFrame)

    # Get line assignments and capacities
    line_map = get_line_assignments()
    line_capacities = get_line_capacities()

    # Build lookup dictionaries
    station_coords = Dict{Int, Tuple{Float64, Float64}}()  # id -> (lat, lon)
    station_names = Dict{Int, String}()  # id -> name

    for row in eachrow(stations)
        station_coords[row.stationID] = (row.lat, row.lon)
        station_names[row.stationID] = row.name
    end

    println("Building arcs from neighbor relationships...")

    # Store arcs: (origin_id, dest_id) -> (distance_km, travel_time, line)
    arcs = Dict{Tuple{Int, Int}, Tuple{Float64, Int, String}}()

    for row in eachrow(stations)
        origin_id = row.stationID
        neighbors = parse_neighbors(string(row.neighbour))

        origin_lat, origin_lon = station_coords[origin_id]

        for dest_id in neighbors
            if !haskey(station_coords, dest_id)
                @warn "Neighbor station $dest_id not found for station $origin_id"
                continue
            end

            dest_lat, dest_lon = station_coords[dest_id]

            # Calculate distance and travel time
            distance = haversine(origin_lat, origin_lon, dest_lat, dest_lon)
            travel_time = estimate_travel_time(distance)

            # Determine line
            line = get_arc_line(origin_id, dest_id, line_map)

            # Store arc (both directions should be present due to neighbor list)
            arc_key = (origin_id, dest_id)
            if !haskey(arcs, arc_key)
                arcs[arc_key] = (distance, travel_time, line)
            end
        end
    end

    println("Found $(length(arcs)) directed arcs")

    # Build output DataFrame
    arc_data = DataFrame(
        origin = String[],
        destination = String[],
        capacity = Int[],
        traveltime = Int[],
        category = String[]
    )

    for ((origin_id, dest_id), (distance, travel_time, line)) in arcs
        origin_name = format_station_name(station_names[origin_id])
        dest_name = format_station_name(station_names[dest_id])
        capacity = get(line_capacities, line, 200)

        push!(arc_data, (
            origin = origin_name,
            destination = dest_name,
            capacity = capacity,
            traveltime = travel_time,
            category = line
        ))
    end

    # Add transfer connections at multi-line stations
    println("Adding transfer connections...")

    # Find stations with multiple lines
    station_lines = Dict{Int, Set{String}}()
    for ((origin_id, dest_id), (_, _, line)) in arcs
        if !haskey(station_lines, origin_id)
            station_lines[origin_id] = Set{String}()
        end
        push!(station_lines[origin_id], line)
    end

    transfer_time = 3  # minutes
    transfer_capacity = 500  # pedestrian capacity

    for (station_id, lines) in station_lines
        if length(lines) > 1
            station_name = format_station_name(station_names[station_id])
            lines_list = collect(lines)

            for i in 1:length(lines_list)
                for j in (i+1):length(lines_list)
                    origin_platform = station_name * "_" * replace(lines_list[i], " " => "")
                    dest_platform = station_name * "_" * replace(lines_list[j], " " => "")

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

    # Remove duplicates
    unique!(arc_data, [:origin, :destination, :category])

    # Save
    mkpath(dirname(output_file))
    CSV.write(output_file, arc_data)

    println("\nDone! Saved $(nrow(arc_data)) arcs to $output_file")
    println("  - $(count(r -> r.category != "Pedestrian Path", eachrow(arc_data))) train arcs")
    println("  - $(count(r -> r.category == "Pedestrian Path", eachrow(arc_data))) transfer arcs")

    # Print distance statistics
    train_arcs = filter(r -> r.category != "Pedestrian Path", arc_data)
    println("\nTravel time statistics:")
    println("  - Min: $(minimum(train_arcs.traveltime)) minutes")
    println("  - Max: $(maximum(train_arcs.traveltime)) minutes")
    println("  - Mean: $(round(mean(train_arcs.traveltime), digits=1)) minutes")

    return arc_data
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Shanghai Metro Arcs Builder")
    println("=" ^ 50)
    println()

    station_file = joinpath(DATA_DIR, "stationInfo.csv")
    if !isfile(station_file)
        error("Station info file not found: $station_file")
    end

    build_metroarcs()
end

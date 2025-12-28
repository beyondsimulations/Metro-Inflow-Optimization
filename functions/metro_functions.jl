"""
    load_demand(config::RegionConfig)

Purpose: Loads passenger demand data from CSV files for the specified date range and time window.
Details: Reads multiple daily CSV files from the configured region directory, combines them, and filters the data to only include entries between start_time and end_time.
Returns: A DataFrame containing the filtered demand data.
"""
function load_demand(config::RegionConfig)
    demand = DataFrame()

    # Iterate over each day in the date range
    @showprogress desc="Loading OD files..." for day in eachindex(daterange)
        filepath = get_od_filepath(config, daterange[day])

        # Check if file exists
        if !isfile(filepath)
            @warn "OD file not found, skipping: $filepath"
            continue
        end

        # Read and concatenate demand data
        day_demand = CSV.read(filepath, DataFrame)
        demand = isempty(demand) ? day_demand : vcat(demand, day_demand)
    end

    if isempty(demand)
        error("No demand data found for date range $(daterange[1]) to $(daterange[end]) in $(config.base_dir)")
    end

    # Filter the demand data to include only the data within the specified time range (start_time and end_time)
    filter!(row -> row.datetime >= start_time, demand)
    filter!(row -> row.datetime <= end_time, demand)

    return demand::DataFrame
end

"""
    aggregate_demand(config::RegionConfig)

Purpose: Processes raw demand data into a summarized format grouped by origin, destination, and time period.
Details: Calls load_demand() to get raw data, maps each minute to a specific period, and then groups and sums demand values.
Returns: A DataFrame with aggregated demand values by origin, destination, and time period.
"""
function aggregate_demand(config::RegionConfig)
    println("Preparing demand data...")

    # Load the demand data for the specified time range
    demand = load_demand(config)
    println("Loaded $(nrow(demand)) demand records.")

    # Create a mapping from periods to date ranges
    timesteps = start_time:Minute(1):end_time
    dict_periods = Dict(timesteps[x] => divrem(x, minutes_in_period)[1] + 1 for x in eachindex(timesteps))

    # Assign each row in the demand DataFrame to a particular period based on its datetime value
    println("Assigning periods...")
    demand.period .= 0
    @showprogress desc="Mapping to periods..." for row in eachrow(demand)
        row.period = dict_periods[row.datetime]
    end

    # Group the demand DataFrame by origin, destination, and period
    println("Aggregating demand...")
    demand = groupby(demand, [:origin, :destination, :period])

    # Sum the value column for each grouping
    demand = combine(demand, :value => sum => :value)
    println("Aggregated to $(nrow(demand)) OD-period combinations.")

    return demand
end

"""
    compute_shift(demand_od)

Purpose: Maps passenger journeys across the metro network over time.
Details: Builds a directed graph of the metro network, calculates shortest paths between all stations, and tracks which arcs (connections) passengers would use in each minute based on their origin, destination, and entry period.
         Only computes shift for OD pairs that have non-zero demand (optimization).
Returns: Four items: shift data (optimized mapping of passengers to arcs over time), original shift data (unoptimized mapping), shift_start_end (mapping of journeys by origin-destination), and the longest path length in the network.
"""
function compute_shift(demand_od::Array{Float64,3})
    println("Preparing shifting data...")

    # Precompute which OD pairs have any demand across all periods
    # Skip OD pairs with zero demand to save computation and memory
    has_demand = [sum(demand_od[o, d, :]) > 0 for o in 1:nr_nodes, d in 1:nr_nodes]
    active_pairs = sum(has_demand)
    total_pairs = nr_nodes * (nr_nodes - 1)
    println("Active OD pairs: $active_pairs / $total_pairs ($(round(100 * active_pairs / total_pairs, digits=1))%)")

    # Prepare the graph of the metro network
    G = DiGraph(nr_nodes)
    distance_matrix = zeros(Int64, nr_nodes, nr_nodes)
    @showprogress desc="Adding arcs..." for arc in metroarcs
        add_edge!(G, d_node_id[arc.origin], d_node_id[arc.destination])
        distance_matrix[d_node_id[arc.origin], d_node_id[arc.destination]] = ceil(Int64, arc.traveltime)
    end

    longest_path = 0
    shift_original = [Tuple{Int64,Int64,Int64}[] for _ in 1:nr_arcs, _ in 1:nr_minutes]
    shift_start_end = [Tuple{Int64,Int64}[] for _ in 1:nr_nodes, _ in 1:nr_nodes]

    @showprogress desc="Computing shortest paths..." for origin in 1:nr_nodes # origin at which people are allowed into the metro
        state = dijkstra_shortest_paths(G, origin, distance_matrix; trackvertices=true)
        all_distances = state.dists
        all_paths = enumerate_paths(state) # holds all nodes on the path from origin to destination
        for path in all_paths
            if path != []
                minute_at_arc_start = 1
                previous_stop = origin # initialize the first stop as origin of the path
                for next_node in path # we now iterate over all nodes on a path
                    if next_node != origin # we skip the first node in the path, as we need the arc
                        arc_id = d_arc_id[(previous_stop, next_node)]
                        push!(shift_start_end[origin, path[end]], (arc_id, minute_at_arc_start))
                        minute_at_arc_start += metroarcs[arc_id].traveltime # add the traveltime to the arc
                    end
                    previous_stop = next_node # set the previous node to the current node
                end
            end
        end
        longest_path = max(maximum(all_distances), longest_path)
        for period in 1:nr_periods # period in which people are allowed into the metro
            for destination in 1:nr_nodes # destination that the people allowed at the origin want to reach
                if origin != destination && has_demand[origin, destination] # skip zero-demand OD pairs
                    previous_stop = origin # initialize the first stop as origin of the path
                    minute_at_arc_start = (period - 1) * minutes_in_period + 1 # first minute of the current period
                    start_at_origin = minute_at_arc_start # used for the assertation of the path
                    for next_node in all_paths[destination] # we now iterate over all nodes on a path
                        if next_node != origin # we skip the first node in the path, as we need the arc
                            arc_id = d_arc_id[(previous_stop, next_node)]
                            traveltime = metroarcs[arc_id].traveltime
                            upper_minute = min(minute_at_arc_start + minutes_in_period - 1, nr_minutes)
                            for minute in minute_at_arc_start:upper_minute
                                push!(shift_original[arc_id, minute], (origin, destination, period))
                            end
                            minute_at_arc_start += traveltime
                            previous_stop = next_node # set the previous node to the current node
                            if next_node == destination
                                expected_distance = start_at_origin + all_distances[destination]
                                @assert expected_distance == minute_at_arc_start "Shortest Paths: Duration computed $minute_at_arc_start, expected $expected_distance"
                            end
                        end
                    end
                end
            end
        end
    end

    # Build shift directly — keep only LAST occurrence of each unique pattern per arc
    shift = [Tuple{Int64,Int64,Int64}[] for _ in 1:nr_arcs, _ in 1:nr_minutes]
    @showprogress desc="Deduplicating shifts..." for a in 1:nr_arcs
        seen = Set{Vector{Tuple{Int64,Int64,Int64}}}()
        # Scan backwards: keep only LAST occurrence of each pattern
        for t in nr_minutes:-1:1
            content = shift_original[a, t]
            if !isempty(content) && content ∉ seen
                push!(seen, content)
                shift[a, t] = content  # Share reference, not copy
            end
        end
    end

    return shift, shift_original, shift_start_end, longest_path
end

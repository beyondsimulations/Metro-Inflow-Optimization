"""
    load_demand()

Purpose: Loads passenger demand data from CSV files for the specified date range and time window.
Details: Reads multiple daily CSV files, combines them, and filters the data to only include entries between start_time and end_time.
Returns: A DataFrame containing the filtered demand data.
"""
function load_demand()
    demand = []

    # Iterate over each day in the date range
    for day in eachindex(daterange)
        # If this is the first day, read the demand data from the corresponding CSV file
        if day == 1
            demand = CSV.read("data_demand/OD_$(daterange[day]).csv", DataFrame)

            # Otherwise, concatenate the demand data from the previous days with the current day's data
        else
            demand = vcat(demand, CSV.read("data_demand/OD_$(daterange[day]).csv", DataFrame))
        end
    end

    # Filter the demand data to include only the data within the specified time range (start_time and end_time)
    filter!(row -> row.datetime >= start_time, demand)
    filter!(row -> row.datetime <= end_time, demand)

    return demand::DataFrame
end

"""
    aggregate_demand()

Purpose: Processes raw demand data into a summarized format grouped by origin, destination, and time period.
Details: Calls load_demand() to get raw data, maps each minute to a specific period, and then groups and sums demand values.
Returns: A DataFrame with aggregated demand values by origin, destination, and time period.
"""
function aggregate_demand()

    # Load the demand data for the specified time range
    demand = load_demand()

    # Create a mapping from periods to date ranges
    timesteps = start_time:Minute(1):end_time
    dict_periods = Dict(timesteps[x] => divrem(x, minutes_in_period)[1] + 1 for x in eachindex(timesteps))

    # Assign each row in the demand DataFrame to a particular period based on its datetime value
    demand.period .= 0
    for row in eachrow(demand)
        row.period = dict_periods[row.datetime]
    end

    # Group the demand DataFrame by origin, destination, and period
    demand = groupby(demand, [:origin, :destination, :period])

    # Sum the value column for each grouping
    demand = combine(demand, :value => sum => :value)

    return demand
end

"""
    compute_shift()

Purpose: Maps passenger journeys across the metro network over time.
Details: Builds a directed graph of the metro network, calculates shortest paths between all stations, and tracks which arcs (connections) passengers would use in each minute based on their origin, destination, and entry period.
Returns: Four items: shift data (optimized mapping of passengers to arcs over time), original shift data (unoptimized mapping), shift_start_end (mapping of journeys by origin-destination), and the longest path length in the network.
"""
function compute_shift()
    # Prepare the graph of the metro network
    G = DiGraph(nr_nodes)
    distance_matrix = zeros(Int64, nr_nodes, nr_nodes)
    for arc in metroarcs
        add_edge!(G, d_node_id[arc.origin], d_node_id[arc.destination])
        distance_matrix[d_node_id[arc.origin], d_node_id[arc.destination]] = ceil(Int64, arc.traveltime)
    end

    longest_path = 0
    shift_original = [Tuple{Int64,Int64,Int64}[] for _ in 1:nr_arcs, _ in 1:nr_minutes]
    shift_start_end = [Tuple{Int64,Int64}[] for _ in 1:nr_nodes, _ in 1:nr_nodes]
    for origin in 1:nr_nodes # origin at which people are allowed into the metro
        state = dijkstra_shortest_paths(G, origin, distance_matrix; trackvertices=true)
        all_distances = state.dists
        all_paths = enumerate_paths(state) # holds all nodes on the path from origin to destination
        for path in all_paths
            if path != []
                minute_at_arc_start = 1
                previous_stop = origin # initialize the first stop as origin of the path
                for next_node in path # we now iterate over all nodes on a path
                    if next_node != origin # we skip the first node in the path, as we need the arc
                        push!(shift_start_end[origin, path[end]], (d_arc_id[(previous_stop, next_node)], minute_at_arc_start))
                        minute_at_arc_start += metroarcs[d_arc_id[(previous_stop, next_node)]].traveltime # add the traveltime to the arc
                    end
                    previous_stop = next_node # set the previous node to the current node
                end
            end
        end
        longest_path = max(maximum(all_distances), longest_path)
        for period in 1:nr_periods # period in which people are allowed into the metro
            for destination in 1:nr_nodes # destination that the people allowed at the origin want to reach
                if origin != destination # skips cases where no travel is necessary
                    previous_stop = origin # initialize the first stop as origin of the path
                    minute_at_arc_start = (period - 1) * minutes_in_period + 1 # first minute of the current period
                    start_at_origin = minute_at_arc_start # used for the assertation of the path
                    for next_node in all_paths[destination] # we now iterate over all nodes on a path
                        if next_node != origin # we skip the first node in the path, as we need the arc
                            for minute in minute_at_arc_start:(minute_at_arc_start+minutes_in_period-1) # access all minutes of that period
                                if minute <= nr_minutes # make sure that we don't exeed the overall timeframe
                                    push!(shift_original[d_arc_id[(previous_stop, next_node)], minute], (origin, destination, period))
                                end
                            end
                            minute_at_arc_start += metroarcs[d_arc_id[(previous_stop, next_node)]].traveltime # add the traveltime to the arc
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

    shift = copy(shift_original)
    for a in 1:nr_arcs
        for t in 1:nr_minutes
            for tt in t+1:nr_minutes
                if shift[a, tt] == shift[a, t]
                    shift[a, t] = []
                end
            end
        end
    end

    return shift, shift_original, shift_start_end, longest_path
end

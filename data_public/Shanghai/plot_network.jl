#!/usr/bin/env julia
"""
Interactive Shanghai Metro Network Visualization using PlotlyJS.

Usage:
    julia --project=../../metroflow plot_network.jl

Opens an interactive plot in your browser with:
- Hover over stations to see name, lines, and coordinates
- Hover over arcs to see line, travel time, and capacity
- Zoom and pan with mouse
- Click legend items to toggle lines on/off
"""

using CSV
using DataFrames
using PlotlyJS

# Shanghai Metro official line colors
const LINE_COLORS = Dict(
    "Line 1" => "#E4002B",   # Red
    "Line 2" => "#97D700",   # Green
    "Line 3" => "#FFD100",   # Yellow
    "Line 4" => "#652D90",   # Purple
    "Line 5" => "#A6217C",   # Magenta
    "Line 6" => "#E50082",   # Pink
    "Line 7" => "#F57E00",   # Orange
    "Line 8" => "#00A0E9",   # Blue
    "Line 9" => "#87CEEB",   # Light Blue
    "Line 10" => "#C6A5CD", # Light Purple
    "Line 11" => "#8B4513", # Brown
    "Line 12" => "#007A60", # Teal
    "Line 13" => "#E898A9", # Light Pink
    "Line 16" => "#32CD32", # Lime Green
    "Pedestrian" => "#AAAAAA",  # Grey for pedestrian/transfer arcs
    "Metro" => "#808080"    # Gray (fallback)
)

function load_station_data()
    println("Loading station data...")

    # Try to load expanded stations file first
    stations_file = "stations_shanghai.csv"
    if isfile(stations_file)
        println("  Using expanded network: $stations_file")
        station_df = CSV.read(stations_file, DataFrame)

        # Build station data structure from expanded format
        stations = Dict{String, NamedTuple}()
        for row in eachrow(station_df)
            node_id = row.node_id
            node_type = row.node_type
            line = ismissing(row.line) || row.line == "" ? "" : row.line
            station_id = hasproperty(row, :station_id) ? row.station_id : 0

            # For line nodes, extract the line as an array; for others, use station_lines lookup
            if node_type == "line"
                lines = [line]
            elseif node_type == "central"
                # Transfer station - will have multiple lines from line nodes
                lines = String[]  # Will be shown as "Transfer"
            else
                # Regular station - use line from CSV if available
                lines = line == "" ? String[] : [line]
            end

            stations[node_id] = (
                name = row.station_name,
                station_id = station_id,
                lon = row.lon,
                lat = row.lat,
                node_type = node_type,
                line = line,
                lines = lines
            )
        end

        n_regular = count(row -> row.node_type == "regular", eachrow(station_df))
        n_central = count(row -> row.node_type == "central", eachrow(station_df))
        n_line = count(row -> row.node_type == "line", eachrow(station_df))
        println("  Loaded $(length(stations)) nodes (regular: $n_regular, central: $n_central, line: $n_line)")
        return stations
    end

    # Fallback: load from original files
    println("  Using original format: stationInfo.csv")
    station_info = CSV.read("stationInfo.csv", DataFrame)
    station_lines = CSV.read("station_lines_2017.csv", DataFrame)

    # Build lookup: station name -> lines
    name_to_lines = Dict{String, Vector{String}}()
    for row in eachrow(station_lines)
        name = row[Symbol("Station Name")]
        lines_str = row[Symbol("Metro Line(s)")]
        lines = [strip(String(l)) for l in split(lines_str, ",")]
        name_to_lines[name] = lines
    end

    # Build station data structure
    stations = Dict{String, NamedTuple}()
    for row in eachrow(station_info)
        ismissing(row.name) && continue
        name = row.name
        metro_name = "Metro_" * replace(name, r"['\s\-\.]" => "")
        lines = get(name_to_lines, name, String[])
        stations[metro_name] = (
            name = name,
            lon = row.lon,
            lat = row.lat,
            node_type = "regular",
            line = "",
            lines = lines
        )
    end

    println("  Loaded $(length(stations)) stations")
    return stations
end

function load_arc_data()
    println("Loading arc data...")
    arcs = CSV.read("metroarcs_shanghai.csv", DataFrame)
    println("  Loaded $(nrow(arcs)) arcs")
    return arcs
end

function create_detailed_arc_plot(stations, arcs)
    """Create a plot where arcs also have hover information."""
    println("Creating detailed interactive plot with arc hover info...")

    traces = GenericTrace[]

    # Group arcs by line, with Pedestrian at the end
    lines_in_data = unique(arcs.category)
    sort!(lines_in_data, by = l -> begin
        if l == "Pedestrian"
            return 1000  # Put pedestrian last
        end
        m = match(r"Line (\d+)", l)
        isnothing(m) ? 999 : parse(Int, m.captures[1])
    end)

    drawn_arcs = Set{String}()

    for line in lines_in_data
        line_arcs = filter(row -> row.category == line, arcs)
        color = get(LINE_COLORS, line, LINE_COLORS["Metro"])

        # Different styling for pedestrian arcs
        is_pedestrian = line == "Pedestrian"
        line_width = is_pedestrian ? 2 : 4
        line_dash = is_pedestrian ? "dot" : "solid"

        # Draw each arc as a separate trace for proper hover
        first_arc = true
        for row in eachrow(line_arcs)
            orig = get(stations, row.origin, nothing)
            dest = get(stations, row.destination, nothing)
            (isnothing(orig) || isnothing(dest)) && continue

            # For pedestrian arcs, don't deduplicate (show all transfer paths)
            if !is_pedestrian
                arc_key = join(sort([row.origin, row.destination]), "-") * "-" * line
                arc_key in drawn_arcs && continue
                push!(drawn_arcs, arc_key)
            end

            if is_pedestrian
                hover_text = "<b>Pedestrian Transfer</b><br>$(orig.name) ↔ $(dest.name)<br>Walk: $(row.traveltime) min"
            else
                hover_text = "<b>$(line)</b><br>$(orig.name) → $(dest.name)<br>Travel: $(row.traveltime) min<br>Capacity: $(row.capacity) pax/min"
            end

            push!(traces, scatter(
                x = [orig.lon, dest.lon],
                y = [orig.lat, dest.lat],
                mode = "lines",
                name = line,
                line = attr(color = color, width = line_width, dash = line_dash),
                text = [hover_text, hover_text],
                hoverinfo = "text",
                legendgroup = line,
                showlegend = first_arc
            ))
            first_arc = false
        end
    end

    # Draw stations with different styles based on node type
    drawn_stations = Set{String}()

    # Collect nodes by type for different marker styles
    regular_lons, regular_lats, regular_texts, regular_colors = Float64[], Float64[], String[], String[]
    central_lons, central_lats, central_texts = Float64[], Float64[], String[]
    line_lons, line_lats, line_texts, line_colors = Float64[], Float64[], String[], String[]

    for row in eachrow(arcs)
        for station_name in [row.origin, row.destination]
            station_name in drawn_stations && continue

            s = get(stations, station_name, nothing)
            isnothing(s) && continue
            push!(drawn_stations, station_name)

            node_type = hasproperty(s, :node_type) ? s.node_type : "regular"

            station_id = hasproperty(s, :station_id) ? s.station_id : 0

            if node_type == "central"
                push!(central_lons, s.lon)
                push!(central_lats, s.lat)
                hover = "<b>$(s.name)</b><br>Type: Transfer Hub<br>ID: $(station_id)<br>Node: $(station_name)"
                push!(central_texts, hover)
            elseif node_type == "line"
                push!(line_lons, s.lon)
                push!(line_lats, s.lat)
                line_name = hasproperty(s, :line) ? s.line : ""
                color = get(LINE_COLORS, line_name, "#808080")
                push!(line_colors, color)
                hover = "<b>$(s.name)</b><br>Type: Line Node<br>Line: $(line_name)<br>ID: $(station_id)<br>Node: $(station_name)"
                push!(line_texts, hover)
            else
                push!(regular_lons, s.lon)
                push!(regular_lats, s.lat)
                lines_str = isempty(s.lines) ? "Unknown" : join(s.lines, ", ")
                primary_line = isempty(s.lines) ? "Metro" : s.lines[1]
                push!(regular_colors, get(LINE_COLORS, primary_line, "#808080"))
                hover = "<b>$(s.name)</b><br>Lines: $(lines_str)<br>ID: $(station_id)<br>Node: $(station_name)"
                push!(regular_texts, hover)
            end
        end
    end

    # Add regular stations
    if !isempty(regular_lons)
        push!(traces, scatter(
            x = regular_lons,
            y = regular_lats,
            mode = "markers",
            name = "Stations",
            marker = attr(size = 8, color = "white", line = attr(color = regular_colors, width = 2)),
            text = regular_texts,
            hoverinfo = "text"
        ))
    end

    # Add central (transfer hub) nodes
    if !isempty(central_lons)
        push!(traces, scatter(
            x = central_lons,
            y = central_lats,
            mode = "markers",
            name = "Transfer Hubs",
            marker = attr(size = 12, color = "white", line = attr(color = "#333333", width = 3)),
            text = central_texts,
            hoverinfo = "text"
        ))
    end

    # Add line nodes (platform nodes at transfer stations)
    if !isempty(line_lons)
        push!(traces, scatter(
            x = line_lons,
            y = line_lats,
            mode = "markers",
            name = "Platform Nodes",
            marker = attr(size = 6, color = line_colors, line = attr(color = "white", width = 1)),
            text = line_texts,
            hoverinfo = "text"
        ))
    end

    layout = Layout(
        title = attr(
            text = "Shanghai Metro Network (Expanded Model) - Hover for Details",
            font = attr(size = 20)
        ),
        xaxis = attr(
            title = "Longitude",
            scaleanchor = "y",
            scaleratio = 1,
            showgrid = true,
            gridcolor = "#f0f0f0",
            zeroline = false
        ),
        yaxis = attr(
            title = "Latitude",
            showgrid = true,
            gridcolor = "#f0f0f0",
            zeroline = false
        ),
        hovermode = "closest",
        legend = attr(
            title = attr(text = "Metro Lines (click to toggle)"),
            itemclick = "toggle",
            itemdoubleclick = "toggleothers",
            x = 1.02,
            y = 1
        ),
        paper_bgcolor = "white",
        plot_bgcolor = "#fafafa",
        width = 1400,
        height = 1000,
        margin = attr(r = 150)
    )

    return plot(traces, layout)
end

function main()
    stations = load_station_data()
    arcs = load_arc_data()

    println("\nSummary:")
    println("  Stations: $(length(stations))")
    println("  Arcs: $(nrow(arcs))")

    # Count arcs per line
    line_counts = combine(groupby(arcs, :category), nrow => :count)
    sort!(line_counts, :category)
    println("\nArcs per line:")
    for row in eachrow(line_counts)
        println("  $(row.category): $(row.count) arcs")
    end

    # Check for "Metro" fallback (unassigned lines)
    metro_arcs = filter(row -> row.category == "Metro", arcs)
    if nrow(metro_arcs) > 0
        println("\n⚠️  Warning: $(nrow(metro_arcs)) arcs with 'Metro' fallback (no line assigned)")
        for row in eachrow(metro_arcs)
            orig = get(stations, row.origin, nothing)
            dest = get(stations, row.destination, nothing)
            orig_name = isnothing(orig) ? row.origin : orig.name
            dest_name = isnothing(dest) ? row.destination : dest.name
            println("    $orig_name → $dest_name")
        end
    else
        println("\n✓ All arcs have proper line assignments")
    end

    println("\nCreating interactive visualization...")
    p = create_detailed_arc_plot(stations, arcs)

    # Save to HTML and open
    savefig(p, "network_plot.html")
    println("\n✓ Saved to network_plot.html")
    println("  Opening in browser...")

    # Open in default browser
    if Sys.isapple()
        run(`open network_plot.html`)
    elseif Sys.islinux()
        run(`xdg-open network_plot.html`)
    elseif Sys.iswindows()
        run(`start network_plot.html`)
    end

    return p
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

using HTTP
using URIs
using CSV
using DataFrames

"""
Download MTA Subway Origin-Destination Ridership data from NY Open Data.
Data source: https://data.ny.gov/Transportation/MTA-Subway-Origin-Destination-Ridership-Estimate-2/jsu2-fbtj
"""

function download_mta_data(output_file="data/mta_od_ridership.csv"; chunk_size=50000)
    base_url = "https://data.ny.gov/resource/jsu2-fbtj.csv"
    mkpath(dirname(output_file))

    offset = 0
    first_chunk = true
    total_rows = 0

    while true
        url = "$(base_url)?\$limit=$(chunk_size)&\$offset=$(offset)&\$order=:id"
        println("Downloading offset $(format_number(offset))...")

        response = HTTP.get(url)
        chunk = String(response.body)

        if isempty(strip(chunk))
            break
        end

        lines = split(strip(chunk), '\n')
        rows_in_chunk = length(lines) - 1

        open(output_file, first_chunk ? "w" : "a") do f
            if first_chunk
                write(f, chunk)
                first_chunk = false
            else
                if length(lines) > 1
                    write(f, "\n" * join(lines[2:end], "\n"))
                end
            end
        end

        total_rows += rows_in_chunk
        println("  Downloaded $(format_number(total_rows)) rows so far...")

        offset += chunk_size

        if rows_in_chunk < chunk_size
            break
        end
    end

    println("\nDone! Saved $(format_number(total_rows)) rows to $output_file")

    size_mb = filesize(output_file) / (1024 * 1024)
    println("File size: $(round(size_mb, digits=2)) MB")

    return output_file
end

function format_number(n::Integer)
    str = string(n)
    parts = String[]
    while length(str) > 3
        pushfirst!(parts, str[end-2:end])
        str = str[1:end-3]
    end
    pushfirst!(parts, str)
    return join(parts, ",")
end

"""
Download filtered MTA data in chunks to avoid memory issues.

Arguments:
- year, month: Filter by year/month fields
- start_date, end_date: Filter by date range (format: "YYYY-MM-DD")
- chunk_size: Rows per batch (default 50000)
- output_file: Output path (auto-generated if not specified)
"""
function download_mta_filtered(;
    year=nothing,
    month=nothing,
    start_date=nothing,
    end_date=nothing,
    chunk_size=50000,
    output_file=nothing
)
    base_url = "https://data.ny.gov/resource/jsu2-fbtj.csv"

    conditions = String[]
    if !isnothing(year)
        push!(conditions, "year=$year")
    end
    if !isnothing(month)
        push!(conditions, "month=$month")
    end
    if !isnothing(start_date)
        push!(conditions, "timestamp>='$(start_date)T00:00:00'")
    end
    if !isnothing(end_date)
        push!(conditions, "timestamp<='$(end_date)T23:59:59'")
    end

    where_clause = isempty(conditions) ? "" : join(conditions, " AND ")

    if isnothing(output_file)
        parts = filter(!isnothing, [year, month, start_date, end_date])
        suffix = isempty(parts) ? "" : "_" * join(parts, "_")
        output_file = "data/mta_od_ridership$(suffix).csv"
    end

    mkpath(dirname(output_file))

    offset = 0
    first_chunk = true
    total_rows = 0

    while true
        params = Dict(
            "\$limit" => string(chunk_size),
            "\$offset" => string(offset),
            "\$order" => ":id"
        )
        if !isempty(where_clause)
            params["\$where"] = where_clause
        end

        query_string = join(["$(URIs.escapeuri(k))=$(URIs.escapeuri(v))" for (k, v) in params], "&")
        url = "$(base_url)?$(query_string)"

        println("Downloading offset $(format_number(offset))...")

        response = HTTP.get(url)
        chunk = String(response.body)
        response = nothing  # Free response memory

        if isempty(strip(chunk))
            break
        end

        lines = split(strip(chunk), '\n')
        rows_in_chunk = length(lines) - 1

        if rows_in_chunk == 0
            break
        end

        open(output_file, first_chunk ? "w" : "a") do f
            if first_chunk
                write(f, chunk)
                first_chunk = false
            else
                if length(lines) > 1
                    write(f, "\n" * join(lines[2:end], "\n"))
                end
            end
        end

        total_rows += rows_in_chunk
        println("  Downloaded $(format_number(total_rows)) rows so far...")

        chunk = nothing
        lines = nothing
        GC.gc()

        offset += chunk_size

        if rows_in_chunk < chunk_size
            break
        end
    end

    println("\nDone! Saved $(format_number(total_rows)) rows to $output_file")

    size_mb = filesize(output_file) / (1024 * 1024)
    println("File size: $(round(size_mb, digits=2)) MB")

    return output_file
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("MTA Subway Origin-Destination Ridership Data Downloader")
    println("=" ^ 50)
    println()

    download_mta_filtered(start_date="2024-10-07", end_date="2024-10-13", output_file="data/mta_week_oct_2024.csv")
end

"""
Region Configuration Module

Provides configuration loading and path utilities for multi-region metro optimization.
"""

using TOML
using Dates

"""
    RegionConfig

Holds all region-specific configuration settings loaded from a TOML file.

Fields include region identification, file paths, time settings (including analysis_dates
for specifying which dates to process), optimization parameters, and plotting filters.
"""
struct RegionConfig
    # Region identification
    name::String
    display_name::String

    # Paths
    base_dir::String
    metroarcs_file::String
    stations_file::String
    od_file_pattern::String

    # Time settings
    interval_minutes::Int
    closed_hours::Vector{Int}
    date_format::String
    analysis_dates::Vector{Date}

    # Optimization parameters
    safety_factors::Vector{Float64}
    min_enter::Vector{Int}
    max_enter::Vector{Int}
    scaling_factors::Vector{Float64}
    past_periods::Vector{Int}
    minutes_in_period::Vector{Int}
    kind_opt::Vector{String}
    kind_queue::Vector{String}

    # Plotting parameters (optional - if set, filter plots to these specific values)
    plot_safety::Union{Nothing, Float64}
    plot_period_length::Union{Nothing, Int}
    plot_past_periods::Union{Nothing, Int}
    plot_max_enter::Union{Nothing, Int}
    plot_min_enter::Union{Nothing, Int}
end

"""
    load_config(config_path::String) -> RegionConfig

Load and parse a TOML configuration file, returning a RegionConfig struct.

# Arguments
- `config_path`: Path to the TOML configuration file

# Returns
- `RegionConfig`: Parsed configuration struct

# Example
```julia
config = load_config("config/doha.toml")
println(config.display_name)  # "Doha Metro"
```
"""
function load_config(config_path::String)::RegionConfig
    if !isfile(config_path)
        error("Configuration file not found: $config_path")
    end

    toml = TOML.parsefile(config_path)

    # Extract required sections
    region = toml["region"]
    paths = toml["paths"]
    time = toml["time"]
    opt = get(toml, "optimization", Dict())
    plotting = get(toml, "plotting", Dict())

    # Default optimization parameters (for backward compatibility)
    base_interval = time["interval_minutes"]
    default_minutes = [base_interval * i for i in 1:4 if base_interval * i <= 60]

    # Parse analysis dates
    date_fmt = DateFormat(time["date_format"])
    analysis_dates = haskey(time, "analysis_dates") ?
        [Date(d, date_fmt) for d in time["analysis_dates"]] :
        Date[]

    # Helper to get optional plotting parameter
    get_optional(d, key, T) = haskey(d, key) ? convert(T, d[key]) : nothing

    return RegionConfig(
        # Region
        region["name"],
        region["display_name"],
        # Paths
        paths["base_dir"],
        paths["metroarcs_file"],
        paths["stations_file"],
        paths["od_file_pattern"],
        # Time
        time["interval_minutes"],
        convert(Vector{Int}, time["closed_hours"]),
        time["date_format"],
        analysis_dates,
        # Optimization (with defaults)
        convert(Vector{Float64}, get(opt, "safety_factors", [0.9])),
        convert(Vector{Int}, get(opt, "min_enter", [0])),
        convert(Vector{Int}, get(opt, "max_enter", [120])),
        convert(Vector{Float64}, get(opt, "scaling_factors", [1.0])),
        convert(Vector{Int}, get(opt, "past_periods", [4])),
        convert(Vector{Int}, get(opt, "minutes_in_period", default_minutes)),
        convert(Vector{String}, get(opt, "kind_opt", ["linweight"])),
        convert(Vector{String}, get(opt, "kind_queue", ["shift_periods"])),
        # Plotting (optional filters)
        get_optional(plotting, "safety", Float64),
        get_optional(plotting, "period_length", Int),
        get_optional(plotting, "past_periods", Int),
        get_optional(plotting, "max_enter", Int),
        get_optional(plotting, "min_enter", Int)
    )
end

"""
    get_metroarcs_path(config::RegionConfig) -> String

Get the full path to the metro arcs CSV file.

# Example
```julia
path = get_metroarcs_path(config)  # "data_public/Doha/metroarcs_doha.csv"
```
"""
function get_metroarcs_path(config::RegionConfig)::String
    return joinpath(config.base_dir, config.metroarcs_file)
end

"""
    get_stations_path(config::RegionConfig) -> String

Get the full path to the stations CSV file.
"""
function get_stations_path(config::RegionConfig)::String
    return joinpath(config.base_dir, config.stations_file)
end

"""
    get_od_filepath(config::RegionConfig, date::Date) -> String

Get the full path to the OD demand file for a specific date.

# Arguments
- `config`: Region configuration
- `date`: The date for which to get the OD file path

# Example
```julia
path = get_od_filepath(config, Date(2022, 11, 29))
# Returns: "data_public/Doha/OD_2022-11-29.csv"
```
"""
function get_od_filepath(config::RegionConfig, date::Date)::String
    date_str = Dates.format(date, config.date_format)
    filename = replace(config.od_file_pattern, "{date}" => date_str)
    return joinpath(config.base_dir, filename)
end

"""
    is_closed_hour(config::RegionConfig, hour::Int) -> Bool

Check if the metro is closed during the given hour.
"""
function is_closed_hour(config::RegionConfig, hour::Int)::Bool
    return hour in config.closed_hours
end

"""
Region Configuration Module

Provides configuration loading and path utilities for multi-region metro optimization.
"""

using TOML
using Dates

"""
    RegionConfig

Holds all region-specific configuration settings loaded from a TOML file.
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
        time["date_format"]
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

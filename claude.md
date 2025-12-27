# Metro-Inflow-Optimization Project Context

## Overview
This project optimizes passenger inflow into metro systems to prevent overcrowding. It supports multiple regions (Doha, Shanghai) through a TOML-based configuration system.

## Multi-Region Configuration

### Config Files
- `config/doha.toml` - Doha Metro (15-minute intervals, dates: 2022-11-27 to 2022-11-30)
- `config/shanghai.toml` - Shanghai Metro (10-minute intervals, dates: 2017-05-01 to 2017-08-31)

### Key Config Fields
```toml
[region]
name = "doha"                    # Used in output filenames
display_name = "Doha Metro"

[paths]
base_dir = "data_public/Doha"    # Data directory
metroarcs_file = "metroarcs_doha.csv"
stations_file = "stations_doha.csv"
od_file_pattern = "OD_{date}.csv"

[time]
interval_minutes = 15            # Base time interval (15 for Doha, 10 for Shanghai)
closed_hours = [3, 4, 5]         # Metro closure hours
date_format = "yyyy-mm-dd"
```

## Running the Framework

### Parallel Framework (recommended)
```bash
# Doha (default)
julia metro_framework_parallel.jl 15 '2022-11-29T05:00:00' '2022-11-30T04:59:00'

# Shanghai
julia metro_framework_parallel.jl --config config/shanghai.toml 10 '2017-05-08T05:00:00' '2017-05-09T04:59:00'
```

### Simple Framework (quick testing)
Edit `CONFIG_FILE` and dates in `metro_framework.jl`, then:
```bash
julia metro_framework.jl
```

### Shell Script
Edit `CONFIG` variable in `run_parallel_metro.sh`, then:
```bash
./run_parallel_metro.sh
```

## Shanghai Data Preprocessing

See `data_public/Shanghai/README.md` for detailed documentation.

### Quick Start

```bash
cd data_public/Shanghai

# 1. Build network (generates stations_shanghai.csv + metroarcs_shanghai.csv)
julia --project=../../metroflow build_metroarcs.jl

# 2. Verify network visually (opens interactive HTML plot)
julia --project=../../metroflow plot_network.jl

# 3. Transform OD data to daily files
julia --project=../../metroflow transform_od.jl 2017-05-15 2017-05-21

# 4. Verify transformation
julia --project=../../metroflow verify_od.jl 2017-05-15 2017-05-21
```

### Expanded Network Model

The Shanghai network uses an expanded model for transfer stations:
- **Central node**: Entry/exit point (e.g., `Metro_Xujiahui`)
- **Line nodes**: Platform per line (e.g., `Metro_Xujiahui_L01`, `Metro_Xujiahui_L09`)
- **Pedestrian arcs**: 1-minute walking connections between platforms

This properly handles parallel lines sharing tracks (separate capacity) and transfer times.

### Input Files

**Included** (with corrections to original data):
- `stationInfo.csv`: Station coordinates and neighbor relationships
- `station_lines_2017.csv`: Which lines serve each station

**Download required** from [original dataset](https://www.nature.com/articles/s41597-025-05416-8):
- `metroData_ODFlow.csv` (11 GB)
- `metroData_InOutFlow.csv` (217 MB)

**Note**: May 8-9, 2017 have incomplete OD data. Use May 10+ for complete data.

## Data Directory Structure
```
data_public/
├── Doha/
│   ├── OD_2022-11-27.csv ... OD_2022-11-30.csv
│   ├── metroarcs_doha.csv
│   └── stations_doha.csv
└── Shanghai/
    ├── README.md                    # Detailed preprocessing documentation
    ├── metroData_ODFlow.csv         # Raw OD data (11GB, 267M trips)
    ├── metroData_InOutFlow.csv      # Station in/out flows (217MB)
    ├── stationInfo.csv              # Station coords + neighbors (input)
    ├── station_lines_2017.csv       # Station-to-line mapping (input)
    ├── stations_shanghai.csv        # Expanded station nodes (generated)
    ├── metroarcs_shanghai.csv       # Network arcs (generated)
    ├── network_plot.html            # Interactive visualization (generated)
    ├── build_metroarcs.jl           # Network builder script
    ├── plot_network.jl              # Visualization script
    ├── transform_od.jl              # OD transformer script
    ├── verify_od.jl                 # Verification script
    └── OD_*.csv                     # Daily OD files (generated)
```

## Output Files

Results are saved with region prefix to prevent overwrites:
- `results/logfile_{region}_{start_time}_{minutes}.csv` - Aggregated summary
- `results/queues/sim_queues_{region}_*.csv` - Queue simulation data
- `results/arcs/sim_arcs_{region}_*.csv` - Arc utilization data
- `results/log_{region}_mip_*.txt` - Run journal

## Key Source Files

| File | Purpose |
|------|---------|
| `functions/config.jl` | Configuration loading (RegionConfig struct) |
| `functions/metro_functions.jl` | Data loading (load_demand, aggregate_demand) |
| `functions/metro_simulation.jl` | Simulation (uses config.interval_minutes) |
| `functions/metro_heuristic.jl` | Optimization heuristic |
| `functions/metro_model.jl` | JuMP optimization models |
| `metro_data_summary.jl` | Analysis and visualization |

## Important Notes

1. **minutes_in_period must be a multiple of config.interval_minutes**
   - Doha: 15, 30, 45, 60
   - Shanghai: 10, 20, 30, 40, 50, 60

2. **Old data directories removed**: `data_demand/` and `data_metro/` no longer exist. All data is in `data_public/{region}/`.

3. **Renaming old result files**: See `results/README_file_renaming.md` for scripts to add region prefix to legacy files.

4. **Station naming**: Both regions use `Metro_StationName` format (e.g., `Metro_Msheireb`, `Metro_LongcaoRoad`).

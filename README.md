# Metro Inflow Optimization

This repository contains the code for optimizing passenger inflow into metro systems to prevent overcrowding. It supports multiple regions through a TOML-based configuration system.

## Supported Regions

| Region | Interval | Data Period | Data Availability |
|--------|----------|-------------|-------------------|
| Doha Metro | 15 min | Nov 27-30, 2022 | Not included (confidential) |
| Shanghai Metro | 10 min | May-Aug 2017 | Publicly available |

**Note**: Doha Metro data cannot be provided due to confidentiality agreements. The framework supports Doha but users must supply their own data.

## Quick Start

### Prerequisites

- Julia 1.9+
- Required packages (see `metroflow/Project.toml`)

### Running the Framework

```bash
# Doha (default)
julia metro_framework_parallel.jl 15 '2022-11-29T05:00:00' '2022-11-30T04:59:00'

# Shanghai
julia metro_framework_parallel.jl --config config/shanghai.toml 10 '2017-05-15T05:00:00' '2017-05-16T04:59:00'
```

## Directory Structure

```
├── config/                      # Region configuration files
│   ├── doha.toml
│   └── shanghai.toml
├── data_public/
│   ├── Doha/                    # Doha Metro data
│   │   ├── metroarcs_doha.csv
│   │   ├── stations_doha.csv
│   │   └── OD_*.csv
│   └── Shanghai/                # Shanghai Metro data
│       ├── README.md            # Detailed preprocessing docs
│       ├── stationInfo.csv      # Station data (included, with fixes)
│       ├── station_lines_2017.csv
│       ├── stations_shanghai.csv    # Generated
│       ├── metroarcs_shanghai.csv   # Generated
│       └── OD_*.csv                 # Generated
├── functions/                   # Core Julia modules
│   ├── config.jl                # Configuration loading
│   ├── metro_functions.jl       # Data loading utilities
│   ├── metro_model.jl           # Optimization models
│   ├── metro_heuristic.jl       # Heuristic algorithms
│   └── metro_simulation.jl      # Simulation logic
├── results/                     # Output directory
├── metroflow/                   # Julia project environment
├── metro_framework.jl           # Simple framework (single-threaded)
└── metro_framework_parallel.jl  # Parallel framework (recommended)
```

## Shanghai Data Setup

The Shanghai data requires preprocessing. Station data files are included with corrections; OD flow files must be downloaded separately.

### 1. Download OD Data

Download from the [original dataset](https://www.nature.com/articles/s41597-025-05416-8):
- `metroData_ODFlow.csv` (11 GB)
- `metroData_InOutFlow.csv` (217 MB)

Place in `data_public/Shanghai/`.

### 2. Build Network

```bash
cd data_public/Shanghai
julia --project=../../metroflow build_metroarcs.jl
```

This generates:
- `stations_shanghai.csv` - Expanded network nodes
- `metroarcs_shanghai.csv` - Network arcs with capacities

### 3. Visualize Network (optional)

```bash
julia --project=../../metroflow plot_network.jl
```

Opens an interactive HTML plot for data verification.

### 4. Transform OD Data

```bash
julia --project=../../metroflow transform_od.jl 2017-05-15 2017-05-21
```

See `data_public/Shanghai/README.md` for detailed documentation.

## Configuration Parameters

Edit config files in `config/` or pass via command line:

| Parameter | Description |
|-----------|-------------|
| `set_safety` | Safety factor limiting arc capacity |
| `set_max_enter` | Maximum passengers allowed to enter |
| `set_min_enter` | Minimum passengers allowed to enter |
| `set_scaling` | Demand scaling factor |
| `minutes_in_period` | Optimization period length (must be multiple of interval) |

## Output

Results are saved with region prefix:
- `results/logfile_{region}_*.csv` - Aggregated summary
- `results/queues/sim_queues_{region}_*.csv` - Queue data
- `results/arcs/sim_arcs_{region}_*.csv` - Arc utilization

## License

This project is licensed under the MIT License.

## Citation

If you use this code, please cite our upcoming research paper on metro inflow optimization.

For the Shanghai OD data, please cite:
> https://www.nature.com/articles/s41597-025-05416-8

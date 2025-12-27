# Shanghai Metro Data

This directory contains Shanghai Metro network data and preprocessing scripts for the Metro Inflow Optimization framework.

## Data Sources

### Included in This Repository (with corrections)

| File | Description |
|------|-------------|
| `stationInfo.csv` | Station coordinates and neighbor relationships (302 stations) |
| `station_lines_2017.csv` | Station-to-line mapping (which lines serve each station) |

These files contain corrections to the original data (fixed station-line assignments, neighbor relationships).

### Download Required

The large OD flow files must be downloaded from the original data repository:

**Source**: [Shanghai Metro OD Flow Dataset](https://www.nature.com/articles/s41597-025-05416-8)

| File | Size | Description |
|------|------|-------------|
| `metroData_ODFlow.csv` | 11 GB | Origin-destination flow records (267M trips, May-Aug 2017) |
| `metroData_InOutFlow.csv` | 217 MB | Station entry/exit counts per 10-min interval |

Place these files in this directory (`data_public/Shanghai/`) after downloading.

### Generated Data (Output)

| File | Description |
|------|-------------|
| `stations_shanghai.csv` | Expanded station nodes with coordinates |
| `metroarcs_shanghai.csv` | Network arcs with travel times and capacities |
| `OD_YYYY-MM-DD.csv` | Daily OD demand files for the framework |
| `network_plot.html` | Interactive network visualization |

## Preprocessing Pipeline

### Step 1: Edit Station Data (Manual)

Before running the build scripts, verify and edit the input files as needed:

- **`stationInfo.csv`**: Station coordinates and neighbor relationships
- **`station_lines_2017.csv`**: Which metro lines serve each station

### Step 2: Build Network (`build_metroarcs.jl`)

Generates the expanded network model with proper transfer station handling.

```bash
cd data_public/Shanghai
julia --project=../../metroflow build_metroarcs.jl
```

**Runtime**: ~7-10 minutes (processes 11GB OD file for capacity estimation)

**What it does**:
1. Loads station data and line assignments
2. Computes empirical arc capacities from OD flow data
3. Detects and skips triangle shortcuts (see below)
4. Creates expanded network with transfer stations
5. Generates `stations_shanghai.csv` and `metroarcs_shanghai.csv`

**Expanded Network Model**:

Transfer stations (serving 2+ lines) are expanded into multiple nodes:
- **Central node**: Entry/exit point (e.g., `Metro_Xujiahui`)
- **Line nodes**: Platform for each line (e.g., `Metro_Xujiahui_L01`, `Metro_Xujiahui_L09`)

Arc types:
- **Track arcs**: Train connections between stations (line-specific capacity)
- **Pedestrian arcs**: Walking connections within transfer stations (1 min, high capacity)

This model properly handles:
- Parallel lines sharing track segments (separate capacity per line)
- Transfer times between lines (1 minute walking)
- Correct shortest-path routing

**Triangle Detection**:

The script automatically detects "triangle shortcuts" where three stations on the same line form a triangle. The longest edge (likely an incorrect shortcut) is skipped. Example: Line 13 goes Hanzhong Road → Natural History Museum → West Nanjing Road, not directly.

### Step 3: Verify Network (`plot_network.jl`)

Generates an interactive HTML visualization for data verification.

```bash
julia --project=../../metroflow plot_network.jl
```

Opens `network_plot.html` in browser with:
- Hover over stations: name, lines, station ID, node type
- Hover over arcs: line, travel time, capacity
- Click legend items to toggle lines on/off
- Zoom and pan to inspect details
- Grey dotted lines = pedestrian/transfer arcs

### Step 4: Transform OD Data (`transform_od.jl`)

Converts raw OD data to daily files for the framework.

```bash
julia --project=../../metroflow transform_od.jl 2017-05-15 2017-05-21
```

**Note**: May 8-9, 2017 have incomplete data (<2% of usual ridership). Use May 10+ for complete data.

### Step 5: Verify Transformation (`verify_od.jl`)

Compares generated OD files against InOutFlow data.

```bash
julia --project=../../metroflow verify_od.jl 2017-05-15 2017-05-21
```

Expected: ~2% difference is normal (OD contains complete trips only).

## Network Statistics

After preprocessing:

| Metric | Count |
|--------|-------|
| Physical stations | 302 |
| Total nodes | ~440 (including line nodes at transfers) |
| Transfer stations | ~60 |
| Track arcs | ~720 |
| Pedestrian arcs | ~450 |
| Metro lines | 14 (Lines 1-13, 16) |

## Capacity Calculation

Arc capacities are derived from observed OD flows:

1. Route each OD trip through the network via shortest path
2. Aggregate flows per arc per hour
3. Take peak hour flow / 60 = passengers per minute
4. Apply consistent capacity across all arcs on the same line

## File Formats

### stations_shanghai.csv
```csv
node_id,station_id,station_name,node_type,line,lon,lat
Metro_Xujiahui,1234,Xujiahui,central,,121.43,31.19
Metro_Xujiahui_L01,1234,Xujiahui,line,Line 1,121.431,31.19
```

- `node_type`: "regular" (single-line), "central" (transfer hub), "line" (platform node)
- `station_id`: Original ID from stationInfo.csv (for debugging)

### metroarcs_shanghai.csv
```csv
origin,destination,capacity,traveltime,category
Metro_Xujiahui_L01,Metro_HengshanRoad,1129,2,Line 1
Metro_Xujiahui,Metro_Xujiahui_L01,9999,1,Pedestrian
```

- `capacity`: Passengers per minute
- `traveltime`: Minutes
- `category`: Line name or "Pedestrian"

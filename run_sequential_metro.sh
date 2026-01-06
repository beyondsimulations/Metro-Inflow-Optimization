#!/bin/bash

# run_sequential_metro.sh - Script to run metro framework SEQUENTIALLY in a tmux session
# Runs one experiment at a time to minimize memory usage
# Each minutes_in_period value runs after the previous one completes
# Uses threading for faster execution

set -e  # Exit on any error

# =============================================================================
# CONFIGURATION - Edit these values for your overnight run
# =============================================================================

# Region configuration file (change this for different regions)
# Options: config/doha.toml, config/shanghai.toml
CONFIG="config/shanghai.toml"

# Number of threads for Julia (use "auto" for automatic detection)
JULIA_THREADS="auto"

# Read analysis dates from TOML config and convert to time periods
# Each date becomes a 24-hour period: YYYY-MM-DDT05:00:00 to next day T04:59:00
echo "Reading dates from $CONFIG..."
DATES=($(grep 'analysis_dates' "$CONFIG" | sed 's/.*\[\(.*\)\].*/\1/' | tr -d '"' | tr ',' ' '))

if [ ${#DATES[@]} -eq 0 ]; then
    echo "Error: No analysis_dates found in $CONFIG"
    exit 1
fi

# Convert dates to time periods (05:00 to 04:59 next day)
TIME_PERIODS=()
for date in "${DATES[@]}"; do
    # Calculate next day (macOS and Linux compatible)
    next_day=$(date -j -v+1d -f "%Y-%m-%d" "$date" "+%Y-%m-%d" 2>/dev/null || date -d "$date + 1 day" "+%Y-%m-%d")
    TIME_PERIODS+=("${date}T05:00:00.00,${next_day}T04:59:00.00")
done

# Read minutes_in_period values from TOML config
INTERVAL=$(grep '^interval_minutes' "$CONFIG" | sed 's/.*= *\([0-9]*\).*/\1/')
MINUTES_VALUES=($(grep '^minutes_in_period' "$CONFIG" | sed 's/.*\[\(.*\)\].*/\1/' | tr -d ' ' | tr ',' ' '))

if [ ${#MINUTES_VALUES[@]} -eq 0 ]; then
    echo "Error: No minutes_in_period found in $CONFIG"
    exit 1
fi

# Verify all values are divisible by interval_minutes
for minutes in "${MINUTES_VALUES[@]}"; do
    if [ $((minutes % INTERVAL)) -ne 0 ]; then
        echo "Error: minutes_in_period value $minutes is not divisible by interval_minutes ($INTERVAL)"
        exit 1
    fi
done

# Run mode: "bound" for optimized runs, "unbound" for baseline (no optimization)
# Set to "unbound" to generate baseline data for comparison
RUN_MODE="bound"  # Options: "bound", "unbound"

# =============================================================================

# Session name
SESSION_NAME="metro_sequential"

echo "=========================================="
echo "Metro Framework Sequential Execution (Threaded)"
echo "=========================================="
echo "Configuration: $CONFIG"
echo "Julia threads: $JULIA_THREADS"
echo "Run mode: $RUN_MODE"
echo "Time periods: ${#TIME_PERIODS[@]}"
echo "Minutes in period values: ${MINUTES_VALUES[*]}"
echo ""

if [ "$RUN_MODE" = "unbound" ]; then
    # Unbound mode: only 1 run per time period (scaling varies inside)
    total_runs=${#TIME_PERIODS[@]}
    echo "UNBOUND MODE: $total_runs time periods (scaling varies per config)"
else
    total_runs=$((${#TIME_PERIODS[@]} * ${#MINUTES_VALUES[@]}))
    echo "BOUND MODE: $total_runs experiments (sequentially, threaded)"
fi
echo ""

# Function to check if tmux session exists
session_exists() {
    tmux has-session -t "$1" 2>/dev/null
}

# Kill existing session if it exists
if session_exists "$SESSION_NAME"; then
    echo "Killing existing session: $SESSION_NAME"
    tmux kill-session -t "$SESSION_NAME"
fi

# Create the sequential run script
SCRIPT_FILE="/tmp/metro_sequential_runner.sh"
cat > "$SCRIPT_FILE" << 'SCRIPT_HEADER'
#!/bin/bash
set -e

CONFIG="$1"
shift
TIME_PERIODS=("$@")

SCRIPT_HEADER

# Add the minutes values, threads, and run mode to the script
echo "MINUTES_VALUES=(${MINUTES_VALUES[*]})" >> "$SCRIPT_FILE"
echo "JULIA_THREADS=\"${JULIA_THREADS}\"" >> "$SCRIPT_FILE"
echo "RUN_MODE=\"${RUN_MODE}\"" >> "$SCRIPT_FILE"

cat >> "$SCRIPT_FILE" << 'SCRIPT_BODY'

echo "=========================================="
echo "Metro Framework Sequential Runner (Threaded)"
echo "=========================================="
echo "Configuration: $CONFIG"
echo "Julia threads: $JULIA_THREADS"
echo "Run mode: $RUN_MODE"
echo "Started at: $(date)"
echo ""

run_count=0

if [ "$RUN_MODE" = "unbound" ]; then
    # UNBOUND MODE: Only 1 run per time period, uses base interval from config
    total_runs=${#TIME_PERIODS[@]}
    # Use base interval: 10 for Shanghai, 15 for Doha
    if [[ "$CONFIG" == *"shanghai"* ]]; then
        minutes=10
    else
        minutes=15
    fi

    for time_period in "${TIME_PERIODS[@]}"; do
        IFS=',' read -r start_time end_time <<< "$time_period"
        ((run_count++))

        echo "=========================================="
        echo "UNBOUND Experiment $run_count / $total_runs"
        echo "Time period: $start_time to $end_time"
        echo "Started at: $(date)"
        echo "=========================================="
        echo ""

        julia -t "$JULIA_THREADS" metro_framework_parallel.jl --config "$CONFIG" --unbound "$minutes" "$start_time" "$end_time"

        echo ""
        echo "Completed at: $(date)"
        echo ""
        echo "Memory cleanup..."
        sleep 10
    done
else
    # BOUND MODE: Full optimization with all minutes values
    total_runs=$((${#TIME_PERIODS[@]} * ${#MINUTES_VALUES[@]}))

    for time_period in "${TIME_PERIODS[@]}"; do
        IFS=',' read -r start_time end_time <<< "$time_period"

        for minutes in "${MINUTES_VALUES[@]}"; do
            ((run_count++))

            echo "=========================================="
            echo "BOUND Experiment $run_count / $total_runs"
            echo "Minutes in period: $minutes"
            echo "Time period: $start_time to $end_time"
            echo "Started at: $(date)"
            echo "=========================================="
            echo ""

            julia -t "$JULIA_THREADS" metro_framework_parallel.jl --config "$CONFIG" "$minutes" "$start_time" "$end_time"

            echo ""
            echo "Completed at: $(date)"
            echo ""
            echo "Memory cleanup..."
            sleep 10  # Brief pause between runs for memory cleanup
        done
    done
fi

echo "=========================================="
echo "ALL EXPERIMENTS COMPLETED"
echo "Finished at: $(date)"
echo "=========================================="
SCRIPT_BODY

chmod +x "$SCRIPT_FILE"

# Create new detached tmux session
echo "Creating tmux session: $SESSION_NAME"
tmux new-session -d -s "$SESSION_NAME"

# Set up the session
tmux send-keys -t "$SESSION_NAME" "cd $(pwd)" Enter
tmux send-keys -t "$SESSION_NAME" "echo 'Starting sequential metro framework runs...'" Enter
tmux send-keys -t "$SESSION_NAME" "echo ''" Enter

# Build the command with time periods as arguments
CMD="$SCRIPT_FILE '$CONFIG'"
for period in "${TIME_PERIODS[@]}"; do
    CMD="$CMD '$period'"
done

tmux send-keys -t "$SESSION_NAME" "$CMD" Enter

echo ""
echo "Sequential runner started in tmux session: $SESSION_NAME"
echo ""
echo "Commands:"
echo "  Attach to session:  tmux attach-session -t $SESSION_NAME"
echo "  Detach from session: Ctrl+B, then D"
echo "  Kill session:       tmux kill-session -t $SESSION_NAME"
echo "  List sessions:      tmux list-sessions"
echo ""
echo "Results will be saved to: results/"
echo ""

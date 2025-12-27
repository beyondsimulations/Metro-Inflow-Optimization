#!/bin/bash

# run_sequential_metro.sh - Script to run metro framework SEQUENTIALLY in a tmux session
# Runs one experiment at a time to minimize memory usage
# Each minutes_in_period value runs after the previous one completes

set -e  # Exit on any error

# Region configuration file (change this for different regions)
# Options: config/doha.toml, config/shanghai.toml
CONFIG="config/shanghai.toml"

# Define multiple time periods to analyze
# Each element should be in format "start_time,end_time"
# Doha dates:
# TIME_PERIODS=(
#     "2022-11-27T05:00:00.00,2022-11-28T04:59:00.00"
#     "2022-11-28T05:00:00.00,2022-11-29T04:59:00.00"
#     "2022-11-29T05:00:00.00,2022-11-30T04:59:00.00"
# )
# Shanghai dates (after running transform_od.jl):
TIME_PERIODS=(
    "2017-05-15T05:00:00.00,2017-05-16T04:59:00.00"
)

# Array of minutes_in_period values to test (you can modify this list)
# For Doha (15-min intervals): use multiples of 15
# For Shanghai (10-min intervals): use multiples of 10
MINUTES_VALUES=(10 20 30 40 50 60)

# Session name
SESSION_NAME="metro_sequential"

echo "=========================================="
echo "Metro Framework Sequential Execution"
echo "=========================================="
echo "Configuration: $CONFIG"
echo "Time periods: ${#TIME_PERIODS[@]}"
echo "Minutes in period values: ${MINUTES_VALUES[*]}"
echo ""

total_runs=$((${#TIME_PERIODS[@]} * ${#MINUTES_VALUES[@]}))
echo "Total experiments to run: $total_runs (sequentially)"
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

# Add the minutes values to the script
echo "MINUTES_VALUES=(${MINUTES_VALUES[*]})" >> "$SCRIPT_FILE"

cat >> "$SCRIPT_FILE" << 'SCRIPT_BODY'

echo "=========================================="
echo "Metro Framework Sequential Runner"
echo "=========================================="
echo "Configuration: $CONFIG"
echo "Started at: $(date)"
echo ""

run_count=0
total_runs=$((${#TIME_PERIODS[@]} * ${#MINUTES_VALUES[@]}))

for time_period in "${TIME_PERIODS[@]}"; do
    IFS=',' read -r start_time end_time <<< "$time_period"

    for minutes in "${MINUTES_VALUES[@]}"; do
        ((run_count++))

        echo "=========================================="
        echo "Experiment $run_count / $total_runs"
        echo "Minutes in period: $minutes"
        echo "Time period: $start_time to $end_time"
        echo "Started at: $(date)"
        echo "=========================================="
        echo ""

        julia metro_framework_parallel.jl --config "$CONFIG" "$minutes" "$start_time" "$end_time"

        echo ""
        echo "Completed at: $(date)"
        echo ""
        echo "Memory cleanup..."
        sleep 5  # Brief pause between runs for memory cleanup
    done
done

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

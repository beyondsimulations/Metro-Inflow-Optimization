#!/bin/bash

# run_parallel_metro.sh - Script to run metro framework in parallel tmux sessions with windows
# Each time period gets its own tmux session with multiple windows for different minutes_in_period values

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

# Base session name
SESSION_BASE="metro_framework"

echo "Starting Metro Framework parallel execution..."
echo "Configuration: $CONFIG"
echo "Time periods to analyze: ${#TIME_PERIODS[@]}"
for i in "${!TIME_PERIODS[@]}"; do
    IFS=',' read -r start_time end_time <<< "${TIME_PERIODS[$i]}"
    echo "  $((i+1)). $start_time to $end_time"
done
echo "Minutes in period values: ${MINUTES_VALUES[*]}"
echo ""

total_windows=$((${#TIME_PERIODS[@]} * ${#MINUTES_VALUES[@]}))
echo "Total sessions: ${#TIME_PERIODS[@]}"
echo "Windows per session: ${#MINUTES_VALUES[@]}"
echo "Total windows: $total_windows"
echo ""

# Function to check if tmux session exists
session_exists() {
    tmux has-session -t "$1" 2>/dev/null
}

# Function to create a tmux session with multiple windows
create_session_with_windows() {
    local start_time=$1
    local end_time=$2
    local time_index=$3

    # Create session name
    local session_name="${SESSION_BASE}_period${time_index}"

    echo "Creating tmux session: $session_name"
    echo "  Time period: $start_time to $end_time"

    # Kill existing session if it exists
    if session_exists "$session_name"; then
        echo "  Killing existing session: $session_name"
        tmux kill-session -t "$session_name"
    fi

    # Create new detached session with first window
    tmux new-session -d -s "$session_name"

    # Set up each window for different minutes_in_period values
    for i in "${!MINUTES_VALUES[@]}"; do
        local minutes="${MINUTES_VALUES[$i]}"
        local window_name="${minutes}min"

        if [ $i -eq 0 ]; then
            # Rename the first window (already exists)
            tmux rename-window -t "$session_name:0" "$window_name"
            window_target="$session_name:0"
        else
            # Create new window
            tmux new-window -t "$session_name" -n "$window_name"
            window_target="$session_name:$window_name"
        fi

        echo "    Window $((i+1))/${#MINUTES_VALUES[@]}: $window_name"

        # Set up the window
        tmux send-keys -t "$window_target" "cd $(pwd)" Enter
        tmux send-keys -t "$window_target" "echo '=========================================='" Enter
        tmux send-keys -t "$window_target" "echo 'Metro Framework - Time Period $time_index'" Enter
        tmux send-keys -t "$window_target" "echo 'Configuration: $CONFIG'" Enter
        tmux send-keys -t "$window_target" "echo 'Minutes in period: $minutes'" Enter
        tmux send-keys -t "$window_target" "echo 'Start time: $start_time'" Enter
        tmux send-keys -t "$window_target" "echo 'End time: $end_time'" Enter
        tmux send-keys -t "$window_target" "echo '=========================================='" Enter
        tmux send-keys -t "$window_target" "echo 'Starting computation...'" Enter

        # Run the Julia script with config
        tmux send-keys -t "$window_target" "julia metro_framework_parallel.jl --config $CONFIG $minutes '$start_time' '$end_time'" Enter

        # Small delay between window setups
        sleep 0.5
    done

    echo "  Session $session_name created with ${#MINUTES_VALUES[@]} windows"
}

# Create all sessions with windows
echo "Creating sessions..."
echo ""

for time_index in "${!TIME_PERIODS[@]}"; do
    IFS=',' read -r start_time end_time <<< "${TIME_PERIODS[$time_index]}"
    create_session_with_windows "$start_time" "$end_time" "$((time_index + 1))"
    echo ""
    sleep 1  # Small delay between session creation
done

echo "All ${#TIME_PERIODS[@]} sessions created successfully!"
echo ""

# Show all created sessions and windows
echo "Created sessions and windows:"
for time_index in "${!TIME_PERIODS[@]}"; do
    IFS=',' read -r start_time end_time <<< "${TIME_PERIODS[$time_index]}"
    session_name="${SESSION_BASE}_period$((time_index + 1))"
    echo "📁 $session_name ($start_time to $end_time)"
    for minutes in "${MINUTES_VALUES[@]}"; do
        echo "  └── ${minutes}min"
    done
    echo ""
done

echo "Navigation commands:"
echo "List all sessions:"
echo "  tmux list-sessions"
echo ""
echo "Attach to a session:"
for time_index in "${!TIME_PERIODS[@]}"; do
    session_name="${SESSION_BASE}_period$((time_index + 1))"
    echo "  tmux attach-session -t $session_name"
done
echo ""
echo "Kill specific session:"
for time_index in "${!TIME_PERIODS[@]}"; do
    session_name="${SESSION_BASE}_period$((time_index + 1))"
    echo "  tmux kill-session -t $session_name"
done
echo ""
echo "Kill all sessions:"
echo "  tmux kill-server"
echo ""

# Function to show session status
show_status() {
    echo ""
    echo "Current session status:"

    running_sessions=0
    total_sessions=${#TIME_PERIODS[@]}

    for time_index in "${!TIME_PERIODS[@]}"; do
        session_name="${SESSION_BASE}_period$((time_index + 1))"
        if session_exists "$session_name"; then
            echo "✓ $session_name - RUNNING (${#MINUTES_VALUES[@]} windows)"
            ((running_sessions++))
        else
            echo "✗ $session_name - NOT RUNNING"
        fi
    done

    echo ""
    echo "Running: $running_sessions/$total_sessions sessions"
}

# Show initial status
show_status

echo ""
echo "Press Ctrl+C to exit this script (sessions will continue running)"
echo "Monitoring session status every 30 seconds..."

# Monitor loop (optional - user can Ctrl+C out)
while true; do
    sleep 30
    show_status

    # Check if all sessions are still running
    all_running=true
    for time_index in "${!TIME_PERIODS[@]}"; do
        session_name="${SESSION_BASE}_period$((time_index + 1))"
        if ! session_exists "$session_name"; then
            all_running=false
            break
        fi
    done

    if [ "$all_running" = false ]; then
        echo ""
        echo "Some sessions have completed or crashed."
        echo "Check individual session logs for details."
        break
    fi
done

echo ""
echo "Monitoring stopped. Sessions may still be running."
echo "Use 'tmux list-sessions' to check current status."

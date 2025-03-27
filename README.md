# Metro Inflow Problem

This repository contains the code associated with an upcoming research paper on controlling the inflow into metro systems while preventing overcrowding.

## Project Overview

The primary objective of this project is to optimize metro inflow by balancing passenger flows in urban rail transport systems. The repository is structured as follows:

### Directory Structure

- **`metro_framework.jl`**: Main framework file that sets up the environment and includes other necessary files.
- **`metro_functions.jl`**: Utility functions for data loading and processing.
- **`metro_model.jl`**: Defines the optimization model for the metro system.
- **`metro_heuristic.jl`**: Implements algorithms for solving the optimization problem.
- **`metro_simulation.jl`**: Simulates the metro system based on the optimization results.
- **`metro_visuals.jl`**: Generates visualizations of the results.
- **`metro_data_summary.jl`**: Summarizes and analyzes the input data and results for the research paper.
- **`data_demand/`**: Directory containing demand data files.
- **`data_metro/`**: Directory containing metro system data files.
- **`results`**: Directory for storing optimization results with queues.
- **`results_paper/`**: Directory for storing optimization results we aim to publish in our paper.
- **`visuals/`**: Directory for storing generated visualizations.

## Configuration

You can customize various parameters in the `metro_framework.jl` file to tailor the optimization and simulation to your needs:

- **`set_safety`**: Safety factor limiting the arc capacity.
- **`set_max_enter`**: Maximum number of people allowed to enter from outside.
- **`set_min_enter`**: Minimum number of people allowed to enter.
- **`set_scaling`**: Scaling of the metro queue to test lower or higher demand.
- **`set_past_periods`**: Timeframe to consider from the past during the optimization.
- **`set_kind_opt`**: Kind of optimization ("regularSqr", "linweight").
- **`set_kind_queue`**: Kind of queue ("shift_periods", "lag_periods").
- **`kind_sim`**: Kind of simulation ("bound", "inflow", "unbound").
- **`minutes_in_period`**: Minutes in each period (in 15-minute intervals).
- **`start_time`**: Start time of the observed time horizon.
- **`end_time`**: End time of the observed time horizon.

## Results

The results of the optimization and simulation will be saved in the `results/` directory. Visualizations will be generated and stored in the `visuals/` directory.

## License

This project is licensed under the MIT License.

function plot_optimization!(results_arcs)
    ```
    plot_optimization!(results_arcs)
    Purpose: Creates animated visualizations of arc utilization from optimization results.
    Details: Generates two animated GIFs showing:
    1. Aggregated utilization: How the overall utilization of metro connections evolves over time
    2. Period utilization: How the utilization within specific time periods changes

    The function:
    - Groups arc data by timestamp
    - For each timestamp, creates a bar chart showing utilization levels
    - Color-codes bars by metro line
    - Includes a horizontal line showing capacity/safety threshold
    - Saves the animations as GIF files with appropriate naming conventions
    ```
    Plots.scalefontsizes()
    Plots.scalefontsizes(2.2)
    timesteps = unique(results_arcs.datetime)
    gp_results_arcs = groupby(results_arcs, :datetime)
    get_group(gdf, keys...) = gdf[(keys...,)]

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_arcs, dt).utilization_aggregated,
                group=get_group(gp_results_arcs, dt).line,
                size=(1920, 1080),
                title="Real Metroarc Usage at timestep $dt",
                ylims=(0, 2),
                ylabel="Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=25mm,
                xformatter=_ -> "",
            )
            hline!([safety], label="Capacity", linewidth=2, colour=:grey60,)
        end
    end
    gif(anim, "visuals/opt_utilization_fps30_aggregated_$(kind_opt)_$(kind_queue)_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).gif", fps=30)

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_arcs, dt).utilization_period,
                group=get_group(gp_results_arcs, dt).line,
                size=(1920, 1080),
                title="Real Metroarc Usage at timestep $dt",
                ylims=(0, 2),
                ylabel="Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=25mm,
                xformatter=_ -> "",
            )
            hline!([safety], label="Capacity", linewidth=2, colour=:grey60,)
        end
    end
    gif(anim, "visuals/opt_utilization_fps30_period_$(kind_opt)_$(kind_queue)_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).gif", fps=30)
end

function plot_simulation!(sim_queues, sim_arcs, kind_sim)
    ```
    plot_simulation!(sim_queues, sim_arcs, kind_sim)
    Purpose: Creates animated visualizations of simulation results.
    Details: Generates three different animated GIFs:
    1. Queue lengths: Shows how passenger queues at each station evolve over time
    2. Entry allowances: Shows how many passengers are allowed to enter each station per minute
    3. Arc utilization: Shows how the utilization of metro connections changes over time

    The function:
    - Groups data by timestamp
    - For each timestamp, creates appropriate bar charts with color coding
    - Adds titles, legends, and formatting to make visualizations clear
    - Saves the animations as GIF files with names that include configuration parameters

    These visualizations help analyze how passenger flows, queues, and network utilization evolve throughout the simulation period under different optimization strategies.
    ```
    timesteps = unique(sim_arcs.datetime)
    gp_results_queues = groupby(sim_queues, :datetime)
    gp_results_arcs = groupby(sim_arcs, :datetime)
    get_group(gdf, keys...) = gdf[(keys...,)]

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_queues, dt).queued,
                size=(1920, 1080),
                title="Queue at timestep $dt",
                ylims=(0, 50000),
                label="People in Queue",
                legend=:topright,
                colour=:steelblue2,
                margin=25mm,
                xticks=(1:length(nodes), get_group(gp_results_queues, dt).station),
                xrotation=30,
                xtickfontsize=10,
            )
        end
    end
    gif(anim, "visuals/sim_queues_fp30_$(kind_opt)_$(kind_sim)_$(kind_queue)_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).gif", fps=30)

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_queues, dt).allowed,
                label="People allowed to Enter per Minute",
                size=(1920, 1080),
                title="Entry per Minute at timestep $dt",
                ylims=(0, max_enter * 1.1),
                legend=:topright,
                margin=25mm,
                xticks=(1:length(nodes), get_group(gp_results_queues, dt).station),
                xrotation=30,
                xtickfontsize=10,
                colour=:orange1,
            )
        end
    end
    gif(anim, "visuals/sim_entry_fps30_$(kind_opt)_$(kind_sim)_$(kind_queue)_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).gif", fps=30)

    anim = @animate for dt in timesteps
        begin
            new_plot = bar(
                get_group(gp_results_arcs, dt).utilization,
                group=get_group(gp_results_arcs, dt).line,
                size=(1920, 1080),
                title="Real Metroarc Usage at timestep $dt",
                ylims=(0, 2),
                ylabel="Arc Utilization",
                colour=[:goldenrod1 :aquamarine3 :seashell2 :crimson],
                legend=:bottomright,
                margin=25mm,
                xformatter=_ -> "",
            )
            hline!([safety], label="Capacity", linewidth=2, colour=:grey60,)
        end
    end
    gif(anim, "visuals/sim_utilization_fps30_$(kind_opt)_$(kind_sim)_$(kind_queue)_mip-$(minutes_in_period)_pam-$(past_minutes)-$(start_time).gif", fps=30)
end

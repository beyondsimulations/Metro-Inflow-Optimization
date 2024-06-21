# Build a MetroInstance object for solving the demand-based routing problem
function build_model_instance(cum_demand_od,demand_od)
    # Preparing model instance.
    println("Preparing model instance.")

    # Create a MetroInstance object with the following parameters:
    # * nr_nodes: number of nodes in the network (i.e., metro areas)
    # * nr_arcs: number of arcs between nodes in the network
    # * nr_periods: number of time periods to consider in the problem
    # * nr_minutes: number of minutes per period
    # * getproperty.(metroarcs, :capacity): list of capacities for each arc 
    # * safety: ratio of max. arc utilization
    # * max_enter: maximum number of people that can enter a node at any minute
    # * cum_demand_od: cummulated demand from origin to destination per period
    # * demand_od: demand from origin to destination per period
    # * shift: shift data that maps minutes and periods and arcs
    
    modelInstance = MetroInstance(
        nr_nodes,
        nr_arcs,
        nr_periods,
        nr_minutes,
        getproperty.(metroarcs, :capacity),
        safety,
        max_enter,
        cum_demand_od,
        demand_od,
        shift
    )

    return modelInstance
end

# Function to build the optimization model
function build_optimization_model(modelInstance)
    # Initialize optimization model object
    im = Model(HiGHS.Optimizer)

    # Set attributes for presolving, time limit, and mip gap
    set_attribute(im, "presolve", "on")
    set_attribute(im, "time_limit", 120.0)
    set_attribute(im, "mip_rel_gap", 0.0)

    println("Preparing optimization model.")
    @variable(im, 
        0 .<= X[o=1:modelInstance.nr_nodes,p=1:modelInstance.nr_periods] .<= modelInstance.max_entry_origin * safety
    )

    println("Preparing objective function.")
    @objective(im, Min, 
        sum((sum(modelInstance.cum_demand_od_in_period[o,d,p] for d in 1:modelInstance.nr_nodes) - X[o,p] *  minutes_in_period)^2 for o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods)
    )

    println("Preparing capacity constraints.")
    @constraint(im, capacity[a in 1:modelInstance.nr_arcs,t in 1:modelInstance.nr_minutes, p_shifts in 0:ceil(Int,past_minutes/minutes_in_period); shift[a,t] != []],
        sum(X[o,p] * modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
    )

    return im,X
end

# Function to build the restricted optimization model that only consideres a slice of the overall periods
function build_restricted_optimization_model(modelInstance,past_minutes,minutes_in_period,current_period)
    # Initialize optimization model object
    im = Model(HiGHS.Optimizer)

    # Set attributes for presolving, time limit, and mip gap
    set_attribute(im, "presolve", "on")
    set_attribute(im, "time_limit", 120.0)
    set_attribute(im, "mip_rel_gap", 0.0)

    println("Preparing optimization model.")
    @variable(im, 
        0 .<= X[o=1:modelInstance.nr_nodes,p=1:modelInstance.nr_periods] .<= modelInstance.max_entry_origin * safety
    )

    # Prepare the upper and lower bound of the observed time slice
    lower_period::Int64 = max(1,current_period - floor(past_minutes/minutes_in_period))
    upper_period::Int64 = min(current_period + ceil(past_minutes/minutes_in_period),modelInstance.nr_periods)
    current_minute::Int64 = ceil(minutes_in_period * current_period)

    println("Preparing objective function.")
    @objective(im, Min, 
        sum((sum(modelInstance.cum_demand_od_in_period[o,d,p] for d in 1:modelInstance.nr_nodes) - X[o,p] *  minutes_in_period)^2 for o in 1:modelInstance.nr_nodes, p in lower_period:upper_period)
    )
    println("Preparing capacity constraints for periodical demand.")
    @constraint(im, capacity_period[a in 1:modelInstance.nr_arcs,t in max(1,current_minute-past_minutes):min(modelInstance.nr_minutes,modelInstance.nr_minutes+minutes_in_period), p_shifts in 0:ceil(Int,past_minutes/minutes_in_period); shift[a,t] != []],
        sum(X[o,p] * modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
    )
    return im,X
end
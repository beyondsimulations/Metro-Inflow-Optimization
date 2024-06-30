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
        0 .<= X[o=1:modelInstance.nr_nodes,p=1:modelInstance.nr_periods] .<= modelInstance.max_entry_origin * modelInstance.safety_factor
    )

    println("Preparing objective function.")
    @objective(im, Min, 
        sum((sum(modelInstance.cum_demand_od_in_period[o,d,p] for d in 1:modelInstance.nr_nodes) - (X[o,p] *  modelInstance.minutes_in_period))^2 for o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods)
    )

    println("Preparing capacity constraints.")
    @constraint(im, capacity[a in 1:modelInstance.nr_arcs,t in 1:modelInstance.nr_minutes, p_shifts in 0:ceil(Int,modelInstance.past_minutes/modelInstance.minutes_in_period); shift[a,t] != []],
        sum(X[o,p] * modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
    )

    return im,X
end

# Function to build the restricted optimization model that only consideres a slice of the overall periods
function build_restricted_optimization_model(modelInstance,current_period,queue_period_age,inflow_raw)
    # Initialize optimization model object
    im = Model(HiGHS.Optimizer)

    # Set attributes for presolving, time limit, and mip gap
    set_attribute(im, "presolve", "on")
    set_attribute(im, "time_limit", 120.0)
    set_attribute(im, "mip_rel_gap", 0.0)

    # Prepare the upper and lower bound of the observed time slice
    lower_period::Int64 = max(1,current_period - floor((60+modelInstance.minutes_in_period)/modelInstance.minutes_in_period))
    upper_period::Int64 = min(current_period + ceil((60+modelInstance.minutes_in_period)/modelInstance.minutes_in_period)+1,modelInstance.nr_periods)
    current_minute::Int64 = ceil(modelInstance.minutes_in_period * current_period)

    println("Preparing optimization model.")
    @variable(im, 
        0 .<= X[o=1:modelInstance.nr_nodes,p=lower_period:upper_period] .<= modelInstance.max_entry_origin * modelInstance.safety_factor
    )

    if modelInstance.kind_opt == "regular"
        println("Preparing regular objective function.")
        @objective(im, Min, 
            sum((sum(modelInstance.cum_demand_od_in_period[o,:,p]) - (X[o,p] * modelInstance.minutes_in_period))^2 for o in 1:modelInstance.nr_nodes, p in lower_period:upper_period)
        )
    elseif modelInstance.kind_opt == "weight"
        println("Preparing weighted objective function.")
        @objective(im, Min, 
            sum((sum(modelInstance.cum_demand_od_in_period[o,d,p] * (queue_period_age[o,p]+1) for d in 1:modelInstance.nr_nodes) - X[o,p]  * (queue_period_age[o,p]+1) *  modelInstance.minutes_in_period)^2 for o in 1:modelInstance.nr_nodes, p in lower_period:upper_period)
        )
    elseif modelInstance.kind_opt == "linear"
        println("Preparing linear objective function.")
        @objective(im, Min, 
            sum((sum(modelInstance.cum_demand_od_in_period[o,:,p]) - (X[o,p] * minutes_in_period)) for o in 1:modelInstance.nr_nodes, p in lower_period:upper_period)
        )

        println("Prepare constraint to prevent a negative dispatch.")
        @constraint(
            im, queue[o in 1:modelInstance.nr_nodes, p in 1:modelInstance.nr_periods],
            sum(sum(modelInstance.cum_demand_od_in_period[o,d,p] for d in 1:modelInstance.nr_nodes) - X[o,p] *  modelInstance.minutes_in_period) >= 0
        )
    elseif modelInstance.kind_opt == "linwei"
        println("Preparing linear objective function.")
        @objective(im, Min, 
            sum((sum(modelInstance.cum_demand_od_in_period[o,:,p]) * (queue_period_age[o,p]+1) - (X[o,p] * modelInstance.minutes_in_period  * (queue_period_age[o,p]+1))) for o in 1:modelInstance.nr_nodes, p in lower_period:upper_period)
        )

    end

    if modelInstance.kind_queue == "shift_cum"
        println("Preparing capacity constraints for periodical demand.")
        @constraint(im, capacity_period[a in 1:modelInstance.nr_arcs,t in max(1,current_minute-(modelInstance.minutes_in_period)):min(current_minute+(60+modelInstance.minutes_in_period),modelInstance.nr_minutes), p_shifts in 0:ceil(Int,modelInstance.past_minutes/modelInstance.minutes_in_period); shift[a,t] != []],
            sum(X[o,p] * modelInstance.cum_demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.cum_demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.cum_demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
        )
    elseif modelInstance.kind_queue == "shift_per"
        println("Preparing capacity constraints for periodical demand.")
        @constraint(im, capacity_period[a in 1:modelInstance.nr_arcs,t in max(1,current_minute-(minutes_in_period)):min(current_minute+(longest_path+minutes_in_period),nr_minutes), p_shifts in 0:ceil(Int,modelInstance.past_minutes/minutes_in_period); shift[a,t] != []],
            sum(X[o,p] * modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
        )
    elseif modelInstance.kind_queue == "shift_dyn"
        println("Preparing capacity constraints for periodical demand.")
        @constraint(im, capacity_period[a in 1:modelInstance.nr_arcs,t in max(1,current_minute-(minutes_in_period)):min(current_minute+(longest_path+minutes_in_period),nr_minutes), p_shifts in 0:min(maximum(queue_period_age), ceil(Int,modelInstance.past_minutes/minutes_in_period)); shift[a,t] != []],
            sum(X[o,p] * modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)]/sum(modelInstance.demand_od_in_period[o,:,max(1,p-p_shifts)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,max(1,p-p_shifts)] > 0) <= modelInstance.capacity_arcs[a] * modelInstance.safety_factor
        )
    elseif modelInstance.kind_queue == "lag_static"
        smart_intervall = ceil(Int,modelInstance.past_minutes/modelInstance.minutes_in_period)
        adjusted_capacity_arcs = zeros(Float64,modelInstance.nr_minutes,modelInstance.nr_arcs)
        minute_range =  max(1,current_minute-(modelInstance.minutes_in_period)):min(current_minute+(60+modelInstance.minutes_in_period),modelInstance.nr_minutes)
        for t in minute_range
            adjusted_capacity_arcs[t,:] = modelInstance.capacity_arcs * modelInstance.safety_factor
        end
        for a in 1:modelInstance.nr_arcs
            for t in minute_range
                if shift[a,t] != []
                    for p_shifts in 0:smart_intervall     
                        arcweight = 0.0
                        for (o,d,p) in shift[a,t]
                            if modelInstance.demand_od_in_period[o,d,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)] > 0
                                if p < current_period
                                    arcweight += ceil(inflow_raw[o,p] * modelInstance.demand_od_in_period[o,d,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)]/sum(modelInstance.demand_od_in_period[o,:,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)]),digits=5)
                                else
                                    if sum(modelInstance.demand_od_in_period[o,:,p]) > 0
                                        arcweight += ceil(modelInstance.min_entry_origin * modelInstance.demand_od_in_period[o,d,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)]/sum(modelInstance.demand_od_in_period[o,:,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)]),digits=5)
                                    end
                                end
                            end
                        end
                        if arcweight > adjusted_capacity_arcs[t,a]
                            adjusted_capacity_arcs[t,a] = arcweight
                        end
                    end
                end
            end
        end

        println("Preparing capacity constraints for periodical demand.")
        @constraint(im, capacity_period[a in 1:modelInstance.nr_arcs,t in minute_range, p_shifts in 0:smart_intervall; shift[a,t] != []],
            sum(X[o,p] * modelInstance.demand_od_in_period[o,d,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)]/sum(modelInstance.demand_od_in_period[o,:,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)]) for (o,d,p) in modelInstance.shift[a,t] if modelInstance.demand_od_in_period[o,d,min(current_period,p-queue_period_age[o,p]+p_shifts,modelInstance.nr_periods)] > 0) <= adjusted_capacity_arcs[t,a]
        )
    end

    println("Restrict dispatch increase.")
    @constraint(im, min_dispatch[p in current_period:upper_period, o in 1:modelInstance.nr_nodes; sum(modelInstance.demand_od_in_period[o,:,p]) > 0],
        X[o,p] >= modelInstance.min_entry_origin
        )
    
    println("Restrict dispatch increase.")
    @constraint(im, max_increase[p in max(2,current_period):upper_period, o in 1:modelInstance.nr_nodes; sum(modelInstance.demand_od_in_period[o,:,p]) > 0],
        X[o,p] - X[o,p-1] <= modelInstance.max_entry_origin/2
        )

    return im,X
end
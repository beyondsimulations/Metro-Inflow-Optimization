function create_model()
    inflow_model = Model(HiGHS.Optimizer)
    set_attribute(inflow_model, "presolve", "on")
    set_attribute(inflow_model, "time_limit", 60.0)
    return inflow_model
end

function build_mainmodel(inflow_model,nodes,on_path,max_enter,queue,minutes_in_step,horizon)
    @variable(
        inflow_model, 
        0 .<= X[eachindex(nodes),1:size(on_path,4)].<= max_enter
        )
    for p in 0:(horizon-1)
        @constraint(
            inflow_model,
            [
                n = eachindex(nodes), 
                t1 in (p*minutes_in_step+1):((p+1)*minutes_in_step+1), 
                t2 in (p*minutes_in_step+1):((p+1)*minutes_in_step+1), 
                
                ],
            X[n,t1] == X[n,t2]
            )
        if p > 1
            @constraint(
                inflow_model,
                [
                    n = eachindex(nodes), 
                    t1 in p*minutes_in_step+1, 
                    t2 in (p+1)*minutes_in_step+1, 
                    ],
                X[n,t2] <= X[n,t1] + max_change
            )
        end
    end
    return X
end

function update_metro_model(
    inflow_model,
    X,
    nodes,
    queue,
    on_path,
    safety,
    arc_capacity,
    grapharcs,
    max_change,
    previous_entry,
    )
    @objective(
        inflow_model, 
        Min,
        sum((sum(queue[o,d] for d in eachindex(nodes)) - X[o,t] for o in eachindex(nodes),t in 1:(minutes_in_step*horizon)).^2)
        )
    @constraint(
        inflow_model,
        [
            t in axes(arc_capacity,1),
            a in axes(arc_capacity,2)
            ],
        sum(X[o,m]*(queue[o,d]/(sum(queue[o,e] for e in eachindex(nodes)))) for o in eachindex(nodes), d in eachindex(nodes), m in axes(on_path,5) if on_path[o,d,a,t,m] == true && queue[o,d] > 0) <= arc_capacity[t,a] - (1-safety) * grapharcs.capacity[a]
        )
    @constraint(
        inflow_model,
        [n = eachindex(nodes)],
        X[n,1] <= previous_entry[n] + max_change
        )
end


function new_model()
    inflow_model = Model(HiGHS.Optimizer)
    set_attribute(inflow_model, "presolve", "on")
    set_attribute(inflow_model, "time_limit", 60.0)
    return inflow_model
end







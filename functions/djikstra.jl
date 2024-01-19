function djisktra(
    nodes,
    nodeid,
    grapharcs
    )
    println("Computing distances with Djikstras algorithm.")
    nr_nodes = length(nodes)
    distance = fill(Inf,nr_nodes,nr_nodes)
    previous_nodes = zeros(Int64,nr_nodes,nr_nodes)
    for node in nodes
        distance[nodeid[node],nodeid[node]] = 0.0
        queue = PriorityQueue{Int64, Float64}()
        enqueue!(queue, nodeid[node], 0.01)
        while !(isempty(queue))
            origin = dequeue!(queue)
            for arc in axes(grapharcs,1)
                if nodeid[grapharcs.origin[arc]] == origin
                    destination = nodeid[grapharcs.destination[arc]]
                    new_distance = distance[nodeid[node],origin] + grapharcs.traveltime[arc]
                    if new_distance < distance[nodeid[node],destination]
                        distance[nodeid[node],destination] = new_distance
                        previous_nodes[nodeid[node],destination] = origin
                        if !(haskey(queue, destination))
                            enqueue!(queue,destination,new_distance)
                        else
                            delete!(queue, destination)
                            enqueue!(queue,destination,new_distance)
                        end
                    end
                end
            end
        end
    end
    return distance, 
    previous_nodes
end
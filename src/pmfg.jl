export planar_maximally_filtered_graph
function planar_maximally_filtered_graph(g, distmx)
    #Return the PMFG, the graph with the least possible weight 
    #assumes that the graph is:
    # - complete 
    # - connected 
    #assumes distmx is 
    # - real 
    # - symmetric
    # - non-negative

    N = nv(g)
    E = ne(g)
    #check that g is complete
    if E != N*(N-1)÷2
        correct_E = N*(N-1)÷2
        error("the graph g is not complete. It should have $(correct_E) edges and it only has $E.")
    elseif !is_connected(g) #check that g is connected
        error("the graph g is not connected.")
    end

    #distmx checks
    R = eltype(distmx)
    if !(R <: Real)
        error("the eltype of distmx is not real.")
    elseif !issymmetric(distmx)
        error("the distmx is not symmetric.")
    elseif sum(distmx .< 0) > 0
        error("distmx has negative elements")
    end

    if nv(g) <= 4
        #if g has fewer than four nodes, it can always be embedded on the plane 
        return g 
    end 

    #otherwise, allocate the empty graph we are going to use 
    out_g = SimpleWeightedGraph{T, R}(N)

    #gather edge data
    T = eltype(g)
    edge_data = Matrix{T}(undef, E, 2)
    weight_data = Vector{R}(undef)
    i = 1
    for edge in edges(g)
        u, v = (src(edge), dst(edge))
        w = distmx[u, v]
        edge_data[i, :] = [u v]
        weight_data[i] = w
        i += 1
    end
    
    #get perm - we want to add smallest weights possible 
    perm = sortperm(weight_data)
    permute!(weight_data, perm)
    permute!(edge_data, perm)

    for i ∈ axes(edge_data, 1)
        #try adding the edge to the pmfg 
        u, v = edge_data[i, :]
        w = weight_data[i]
        add_edge!(out_g, u, v, w)
        if !is_planar(out_g)
            rem_edge!(out_g, u, v)
        end
        (ne(out_g) == 3*(N-2)) && break 
    end

    return out_g
    
end

planar_maximally_filtered_graph(g::SimpleWeightedGraph) = planar_maximally_filtered_graph(g, weights(g))
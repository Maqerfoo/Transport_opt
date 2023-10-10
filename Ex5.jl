struct Instance
    n_tasks::Int64
    avg_speed::Float64
    workhours::Float64
    volume_capacity::Float64
    weight_capacity::Float64
    coord_x::Array{Float64,1}
    coord_y::Array{Float64,1}
    service_time::Array{Float64,1}
    parts_weight::Array{Float64,1}
    parts_volume::Array{Float64,1}
    distance::Array{Float64,2}
end

struct Edge
    from::Int64
    to::Int64
    cost::Any
    res_usage::Array{Any,1}
end

struct Graph
    V::Array{Int64,1}
    A::Array{Array{Edge,1},1}
end

struct Label
    vertex::Int64
    s::Int64 # visited count
    visited::Array{Bool,1} # visited set
    cost::Int64 # cost
    prev# previous label, or the value "missing" if none exist
    res_usage::Array{Float64,1}
end

function instance_reader(filename)
    file = open(filename)

    vals = split(readline(file))
    n_tasks = parse(Int64,vals[1])
    avg_speed = parse(Float64,vals[2])
    workhours = parse(Float64,vals[3])
    volume_cap = parse(Float64,vals[4])
    weight_cap = parse(Float64,vals[5])

    coord_x = zeros(n_tasks)
    coord_y = zeros(n_tasks)
    service_time = zeros(n_tasks)
    parts_weight = zeros(n_tasks)
    parts_volume = zeros(n_tasks)
    dist = zeros(n_tasks,n_tasks)

    for t in 1:n_tasks
        vals = split(readline(file))
        coord_x[t] = parse(Float64,vals[1])
        coord_y[t] = parse(Float64,vals[2])
        parts_volume[t] = parse(Float64,vals[3])
        parts_weight[t] = parse(Float64,vals[4])
        service_time[t] = parse(Float64,vals[5])
    end

    for t in 1:n_tasks
        vals = split(readline(file))
        dist[t,:] = parse.(Float64,vals)
    end

    close(file)

    return Instance(n_tasks,avg_speed,workhours,
                    volume_cap, weight_cap,
                    coord_x, coord_y,
                    service_time,parts_weight,
                    parts_volume,dist)
end

function graph_constructor(inst)
    n_vertices = inst.n_tasks
    res_high = [inst.volume_capacity, inst.weight_capacity, inst.workhours]
    res_no = length(res_high)
    prio_no = floor(Int, inst.n_tasks/2)
    prio_mult = vcat(2*ones(prio_no), ones(inst.n_tasks - prio_no))
    println(prio_mult)
    #n_arcs = n_vertices^2
    A = [Edge[] for v in 1:n_vertices+1]
    for i in 1:n_vertices, j in 1:n_vertices 
        if i!= j
            push!(A[i],Edge(i,j,-1*prio_mult[j],[inst.parts_volume[j], inst.parts_weight[j], inst.distance[i,j]/inst.avg_speed + inst.service_time[i]]))
        end
    end
    #Add dummy edge
    for i in 2:n_vertices
        push!(A[i],Edge(i,n_vertices+1,0, [0,0, inst.distance[i,1]/inst.avg_speed + inst.service_time[i]] ))
    end
    return Graph(collect(1:n_vertices+1),A), res_high
end

function dominates(a::Label, b::Label)
    return a.cost <= b.cost && a.res_usage <= b.res_usage && a.visited ⊆ b.visited
end

# Label constructor used to generate the first label
function Label(source,V, res_no)
    v = zeros(Int64,length(V))
    v[source] = 1
    return Label(source,1,v,0,missing, zeros(Float64,res_no))
end

# Extends a label from an edge toward the vertex v
function extend(e::Edge,λ::Label,v::Int64,res_max)
    if λ.visited[v] || any(λ.res_usage + e.res_usage .> res_max)
        return Nothing
    end
    visited = deepcopy(λ.visited)
    visited[v] = 1
    return Label(v,λ.s+1,visited,λ.cost+e.cost,λ,λ.res_usage+e.res_usage)
end

# Applies the dominance rule. Can be improved by keeping the labels sorted
function dominance(Λ::Array{Label,1},l::Label)
    changed = false
    A = []
    if isempty(Λ)
        changed = true
        push!(A,l)
    else
        for λ in Λ
            if dominates(λ,l)
                return Λ, false
            elseif dominates(l,λ)
                if !changed
                    changed = true
                    push!(A,l)
                end
            else
                push!(A,λ)
            end
        end
    end
    if !changed
        changed = true
        push!(A,l)
    end
    return A, changed
end

function labeling_ERCSP(G,source, res_max)

    # INITIALIZE-SINGLE-SOURCE
    # Generate the set of empty labels for each vertex in the graph
    res_no = length(res_max)
    Λ = [Label[] for v in G.V]
    # Add the initialal label for the source
    push!(Λ[source],Label(source,G.V,res_no))

    # MAIN-BODY
    # The set of vertices to extend
    Q = Set{Int64}()
    # Adding the source to the vertices to extend
    push!(Q,source)

    # while we still have vertices to extend
    while !isempty(Q)
        # get the next vertex to extend
        u = pop!(Q)

        # for each edge going out of the vertex
        for e in G.A[u]
            # get the next vertex
            v = e.to
            # flag to check if the labels changed
            changed = false
            # for each label of the vertex to extend
            for λ in Λ[u]
                # extend the label
                l=extend(e,λ,v,res_max)
                if l!=Nothing # we could extend
                    # apply the dominance rule
                    Λ[v], ch = dominance(Λ[v],l)
                    # update the changed flag
                    changed = changed || ch
                end
            end
            # if the labels have been updated
            if changed
                # add vertex v to the list of nodes to extend
                push!(Q,v)
            end
        end
    end
    return Λ
end


inst = instance_reader("Data/instance1.txt")
G, res_max = graph_constructor(inst)
labels = labeling_ERCSP(G,1,res_max)


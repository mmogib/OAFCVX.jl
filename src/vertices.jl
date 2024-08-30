abstract type RandomVertex <: NMVertex end
abstract type IdealPointVertex <: NMVertex end
abstract type InnerPointVertex <: NMVertex end
abstract type ClusterVertex <: NMVertex end

function symbolToVertex(s::Symbol)
    return if s == :RandomVertex
        RandomVertex
    elseif s == :IdealPointVertex
        IdealPointVertex
    elseif s == :InnerPointVertex
        InnerPointVertex
    else
        ClusterVertex
    end
end

struct VertexOptions
    problem::MOAProblem
    Vk::Vector
    V_used::Set
    Iter::Int
    C::Union{Nothing,Set}
    t::Int
end
VertexOptions(vo::VertexOptions, Vk::Vector, V_used::Set, iter::Int) = begin
    if iter == 3
        s = Set()
        VkS = VkToSet(Vk)
        for v in VkS
            push!(s, v)
        end
        return VertexOptions(vo.problem, Vk, V_used, iter, s, vo.t)
    else
        return VertexOptions(vo.problem, Vk, V_used, iter, vo.C, vo.t)
    end
end
VertexOptions(vo::VertexOptions, t::Int) =
    VertexOptions(vo.problem, vo.Vk, vo.V_used, vo.Iter, vo.C, t)



function getVertex(::Type{RandomVertex}, options::VertexOptions)
    Vk = options.Vk
    vt = filter(x -> !x[3], Vk)
    if length(vt) == 0
        return nothing, nothing, options.t
    end
    rv = rand(vt)
    v = Float64.(rv[1])
    vindx = rv[2]
    return v, vindx, options.t
end

function getVertex(::Type{IdealPointVertex}, options::VertexOptions)
    Vk = options.Vk
    vt = filter(x -> !x[3], Vk)
    if length(vt) == 0
        return nothing, nothing, options.t
    end
    yI = options.problem.IdealPoint
    VkwithYI =
        map(Vk) do vk
            return (vk[1], vk[2], vk[3], norm(vk[1] - yI))
        end |> Y -> sort(Y, by = x -> x[4])
    rv = VkwithYI[1]
    v = rv[1]

    vindx = rv[2]
    return v, vindx, options.t
end


function getVertex(::Type{InnerPointVertex}, options::VertexOptions)
    Vk = options.Vk
    vt = filter(x -> !x[3], Vk)
    if length(vt) == 0
        return nothing, nothing, options.t
    end
    p = options.problem.InnerPoint
    VkwithYI =
        map(Vk) do vk
            return (vk[1], vk[2], vk[3], norm(vk[1] - p))
        end |> Y -> sort(Y, by = x -> x[4], rev = true)
    rv = VkwithYI[1]
    v = rv[1]
    vindx = rv[2]
    return v, vindx, options.t
end
function getVertex(::Type{ClusterVertex}, options::VertexOptions)
    # function select_vertex(t::Int, Vk::Vector{Vector{Float64}}, V_used::Set{Vector{Float64}}, 
    #  C::Vector{Vector{Float64}}) :: Tuple{Vector{Float64}, Int}
    # Initialize clusters
    C = options.C
    if isnothing(C)
        return getVertex(RandomVertex, options)
    end
    clusters = [Set{Vector}() for _ = 1:length(C)]

    # Assign vertices to the nearest cluster
    Vk = options.Vk
    Vks = VkToSet(Vk)
    V_used = options.V_used
    t = options.t
    for v in setdiff(Vks, V_used)
        distances = [norm(Float64.(v) - Float64.(ci)) for ci in C]
        i = argmin(distances)
        push!(clusters[i], v)
    end

    # Calculate the current cluster index
    current = (t + 1) % length(C)
    current = current == 0 ? length(C) : current # Adjust for Julia's 1-based indexing

    # Find the next non-empty cluster
    while isempty(clusters[current])
        t += 1
        current = (t + 1) % length(C)
        current = current == 0 ? length(C) : current # Adjust for Julia's 1-based indexing
    end

    # Pick an arbitrary vertex from the current cluster
    v = first(clusters[current])
    t += 1
    vinx = getVkVertexIndex(Vk, v)
    return v, vinx, t
end

export RandomVertex, IdealPointVertex, InnerPointVertex, ClusterVertex, symbolToVertex
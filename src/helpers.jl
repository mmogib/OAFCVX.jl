function getPurityIndex(
    results::Vector{MOASuccessResult};
    vertex::Union{Nothing,Symbol} = :RandomVertex,
    direction::Union{Nothing,Symbol} = nothing,
)
    condition(x::MOASuccessResult) = begin
        cond1 = isnothing(vertex) ? true : x.name.vertex == String(vertex)
        cond2 = isnothing(direction) ? true : x.name.direction == String(direction)
        return cond1 && cond2
    end
    Fp = solutions2Matrix(results)
    random_vertex_solutions = filter(condition, results)
    random_vertex_solutions_total = sum(map(x -> size(x.X, 2), random_vertex_solutions))
    ddim = results[1].problem.ddim
    FpR = Matrix{Number}(undef, ddim, random_vertex_solutions_total)
    i = 0
    for sol in random_vertex_solutions
        X = sol.X
        l = size(X, 2)
        FpR[:, i+1:i+l] = X
        i = i + l
    end
    FpRu = unique(FpR, dims = 2)
    purity_index = size(Fp, 2) / length(FpRu)
    purity_index
end
function solutions2Matrix(results::Vector{MOASuccessResult}; is_unique::Bool = true)
    total_solutsion = sum(map(x -> size(x.X, 2), results))
    ddim = results[1].problem.ddim
    all_solustions = Matrix{Number}(undef, ddim, total_solutsion)
    i = 0
    for sol in results
        X = sol.X
        l = size(X, 2)
        all_solustions[:, i+1:i+l] = X
        i = i + l
    end
    Fp = is_unique ? unique(all_solustions, dims = 2) : all_solustions
    Fp
end

function hausdorffDistance(Vk::Matrix, Opt::Function)
    z = Vector{Number}(undef, size(Vk, 1))
    for (i, v) in enumerate(eachrow(Vk))
        z[i] = Opt(v[:])
    end
    maximum(z)
end

function hausdorffDistance(Vk::Vector, Opt::Function)
    z = Vector{Number}(undef, length(Vk))
    for (i, v) in enumerate(Vk)
        z[i] = Opt(v)
    end
    maximum(z)
end
function Vk2vec(Vk::Vector)
    map(x -> x[1], Vk)
end
function set2vec(S::Set)
    V = Vector{Vector{Number}}(undef, length(S))
    for (i, v) in enumerate(S)
        V[i] = v
    end
    V
end
function VkToSet(Vk::Vector)
    S = Set()
    foreach(Vk) do v
        push!(S, v[1])
    end
    S
end
function getPk(::Type{NMScalar}, W::Matrix{<:Number}, X::Matrix{<:Number}, f::Function)
    m = size(W, 2)
    model = Model()
    @variable(model, y[1:m])
    @constraint(model, pl[i in 1:m], dot(W[:, i], y - f(X[1, i])) >= 0)
    poly = polyhedron(model, CDDLib.Library(:exact))
    removehredundancy!(poly)

    poly
end
function getPk(::Type{NMVector}, W::Matrix{<:Number}, X::Matrix{<:Number}, f::Function)
    m = size(W, 2)
    model = Model()
    @variable(model, y[1:m])
    @constraint(model, pl[i in 1:m], dot(W[:, i], y - f(X[:, i])) >= 0)
    poly = polyhedron(model, CDDLib.Library(:exact))
    removehredundancy!(poly)

    poly
end
function getPk(P::Polyhedron, H::HalfSpace)
    nP = P ∩ H
    # removehredundancy!(nP)
    nP
end



function getVertices(P::Polyhedron)
    hv = vrep(P) |> MixedMatVRep
    hv.V, hv.R
end


function updateVkUsed(Vk::Vector, v::Vector{<:Number}, ϵ::Float64 = 1e-6)
    map(x -> norm(x[1] - v) <= ϵ ? (x[1], x[2], true) : x, Vk)
end
function getVkVertexIndex(Vk::Vector, v::Vector{<:Number}, ϵ::Float64 = 1e-6)
    fv = findfirst(x -> norm(x[1] - v) <= ϵ, Vk)
    Vk[fv][2]
end

function getnewVk(Vk, V::Matrix{<:Number}, ϵ::Float64 = 1e-6)
    newvs = map(enumerate(eachrow(V))) do (i, r)
        row = r[:]
        found = findfirst(x -> norm(x[1] - row) <= ϵ, Vk)
        if isnothing(found)
            (row, i, false)
        else
            (row, i, true)
        end
    end
    newvs

end

# function getnewVk(Vk, V::Matrix{<:Number}, ϵ::Float64 = 1e-2)
#     current_count = length(Vk)
#     newvs = map(eachrow(V)) do r
#         row = Vector(r)
#         found = findfirst(x -> norm(x[1] - row) <= ϵ, Vk)
#         if isnothing(found)
#             current_count += 1
#             (row, current_count, false)
#         else
#             nothing
#         end
#     end
#     newvs = filter(x -> !isnothing(x), newvs)
#     if length(newvs) == 0
#         Vk
#     else
#         # uones = unique(i -> newvs[i][1], eachindex(newvs))
#         Vk = vcat(Vk, newvs)
#         # uniqueVertices(Vk)
#         Vk
#     end
# end
function uniqueVertices(Vk, digs = 6)
    unique(x -> norm(round.(x[1], digits = digs)), Vk)
end



function getVecrtexIndex(v::Vector, poly::Polyhedron{<:Number})
    allinedices = eachindex(points(poly))
    i = findfirst(map(indx -> get(poly, indx) == v, allinedices))
    for (j, indx) in enumerate(allinedices)
        if i == j
            return indx
        end
    end
end

function getAdjacentVecrtices(v::Vector{<:Number}, poly::Polyhedron{<:Number})
    # vs = Matrix{Float64}(undef, npoints(poly), length(v))

    # w = [incidentpoints(poly, hidx) for hidx in eachindex(halfspaces(poly))]

    for hs in halfspaces(poly)
        hp = Polyhedra.hyperplane(hs)
        facet = intersect(poly, hp)
        for indx in eachindex(points(facet))
            w = get(facet, indx)
            if w == v
                println(facet)
            end
        end
        # println("--", facet)
    end
end




export getAdjacentVecrtices,
    getVecrtexIndex, hausdorffDistance, getVertices, getPurityIndex, solutions2Matrix
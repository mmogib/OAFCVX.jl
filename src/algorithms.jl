function checkNewx(x::Vector{<:Number}, X::Matrix{<:Number})
    isnothing(findfirst(y -> y ≈ x, eachcol(X)))
end
function paa(
    direction::Type{T} where {T<:NMDirection},
    vertex::Type{S} where {S<:NMVertex},
    problem::MOAProblem,
    options::MOAOptions,
)
    n = problem.pdim
    m = problem.ddim
    W = problem.Cp.A
    f = problem.f
    # WS = problem.WS
    PS = problem.PS
    DPS = problem.DPS
    ϵ = options.ϵ
    # maxiters = options.maxiters
    X = problem.X

    SC = m
    Pk = n == 1 ? getPk(NMScalar, W, X, f) : getPk(NMVector, W, X, f)
    V, R = getVertices(Pk)
    Vraws = collect(eachrow(V))
    Vk = Vector{Tuple{Vector{Rational{BigInt}},Int,Bool}}(undef, length(Vraws))
    foreach(enumerate(eachrow(V))) do (i, r)
        Vk[i] = (Rational{BigInt}.(r), i, false)
    end
    counter = 1
    V_used = Set()
    V_all = VkToSet(Vk)
    vertexOptions = VertexOptions(problem, Vk, V_used, 1, nothing, 0)
    while true
        card = size(X, 2)
        if isempty(setdiff(V_all, V_used))
            return MOASuccessResult(
                MOAResultName(String(Symbol(direction)), String(Symbol(vertex))),
                Message("Solved due to no new vertices."),
                X,
                Vk2vec(Vk),
                Pk,
                SC,
                card,
                problem,
            )
        end
        if options.stopping_criteria(counter, card)
            return MOASuccessResult(
                MOAResultName(String(Symbol(direction)), String(Symbol(vertex))),
                Message("Solved due to stopping criteria."),
                X,
                Vk2vec(Vk),
                Pk,
                SC,
                card,
                problem,
            )
        end
        counter += 1
        try
            vertexOptions = VertexOptions(vertexOptions, Vk, V_used, counter)
            v, vindx, t = getVertex(vertex, vertexOptions)
            vertexOptions = VertexOptions(vertexOptions, t)
            if isnothing(v)
                break
            end
            push!(V_used, v)
            directionOptoins = DirectionOptions(problem, Vk, R, v, vindx, ϵ)
            d = getDirection(direction, directionOptoins)
            xv, zv = PS(v, d)
            wv = DPS(v, d)
            SC += 2
            Vk = updateVkUsed(Vk, v)
            X = checkNewx(xv, X) ? hcat(X, xv) : X

            V, R, Vk = if zv > ϵ
                H = HalfSpace(-wv, -dot(wv, v) - zv)
                Pk = getPk(Pk, H)
                V, R = getVertices(Pk)
                Vk = getnewVk(Vk, V)
                V_all = VkToSet(Vk)
                V, R, Vk
            else
                V, R, Vk
            end
        catch e
            # println(e)
            msg = if hasfield(typeof(e), :msg)
                e.msg
            elseif isa(e, MethodError)
                "Error calling $(e.f) with $(e.args) ..."
            else
                "Unknown error"
            end
            return MOAFailureResult(msg)
        end
    end
    return MOASuccessResult(
        MOAResultName(String(Symbol(direction)), String(Symbol(vertex))),
        Message("Solved"),
        X,
        Vk2vec(Vk),
        Pk,
        SC,
        size(X, 2),
        problem,
    )
end

# function paaold(
#     direction::Type{T} where {T<:NMDirection},
#     vertex::Type{S} where {S<:NMVertex},
#     problem::MOAProblem,
#     options::MOAOptions,
# )
#     n = problem.pdim
#     m = problem.ddim
#     W = problem.Cp.A
#     f = problem.f
#     # WS = problem.WS
#     PS = problem.PS
#     DPS = problem.DPS
#     ϵ = options.ϵ
#     # maxiters = options.maxiters
#     X = problem.X

#     SC = m
#     Pk = n == 1 ? getPk(NMScalar, W, X, f) : getPk(NMVector, W, X, f)
#     V, R = getVertices(Pk)
#     Vraws = collect(eachrow(V))
#     Vk = Vector{Tuple{Vector{Rational{BigInt}},Int,Bool}}(undef, length(Vraws))
#     foreach(enumerate(eachrow(V))) do (i, r)
#         Vk[i] = (Rational{BigInt}.(r), i, false)
#     end
#     counter = 1
#     while true
#         card = size(X, 2)
#         if options.stopping_criteria(counter, card)
#             return MOASuccessResult(
#                 MOAResultName(String(Symbol(direction)), String(Symbol(vertex))),
#                 Message("Solved due to stopping criteria."),
#                 X,
#                 Vk,
#                 Pk,
#                 SC,
#                 card,
#                 problem,
#             )
#         end
#         counter += 1
#         try
#             vertexOptions = VertexOptions(problem, Vk)
#             v, vindx = getVertex(vertex, vertexOptions)
#             if isnothing(v)
#                 break
#             end
#             directionOptoins = DirectionOptions(problem, Vk, R, v, vindx, ϵ)
#             d = getDirection(direction, directionOptoins)
#             xv, zv = PS(v, d)
#             wv = DPS(v, d)
#             SC += 2
#             Vk = updateVkUsed(Vk, v)
#             X = hcat(X, xv)
#             V, R, Vk = if zv > ϵ
#                 H = HalfSpace(-wv, -dot(wv, v) - zv)
#                 Pk = getPk(Pk, H)
#                 V, R = getVertices(Pk)
#                 Vk = getnewVk(Vk, V)
#                 # println(counter, "::::", Vk, "\n\n")
#                 V, R, Vk
#             else
#                 V, R, Vk
#             end
#         catch e
#             # println(e.msg)
#             msg = if hasfield(typeof(e), :msg)
#                 e.msg
#             elseif isa(e, MethodError)
#                 "Error calling $(e.f) with $(e.args) ..."
#             else
#                 "Unknown error"
#             end
#             return MOAFailureResult(msg)
#         end
#     end
#     return MOASuccessResult(
#         MOAResultName(String(Symbol(direction)), String(Symbol(vertex))),
#         Message("Solved"),
#         X,
#         Vk,
#         Pk,
#         SC,
#         size(X, 2),
#         problem,
#     )
# end

export paa
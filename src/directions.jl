abstract type FixedDirection <: NMDirection end
abstract type FixedPointDirection <: NMDirection end
abstract type IdealPointDirection <: NMDirection end
abstract type ApproximateAdjDirection <: NMDirection end

function symbolToDirection(s::Symbol)
    return if s == :FixedDirection
        FixedDirection
    elseif s == :FixedPointDirection
        FixedPointDirection
    elseif s == :IdealPointDirection
        IdealPointDirection
    else
        ApproximateAdjDirection
    end
end

struct DirectionOptions
    problem::MOAProblem
    Vk::Vector
    R::Array
    v::Vector{<:Number}
    vindx::Int
    ϵ::Float64
end

function getDirection(::Type{FixedDirection}, options::DirectionOptions)
    return options.problem.InteriorDirection
end
function getDirection(::Type{FixedPointDirection}, options::DirectionOptions)
    d = options.problem.InnerPoint - options.v
    if any(d .< 0)
        d = options.problem.InteriorDirection
    end
    return d / norm(d)
end

function getDirection(::Type{IdealPointDirection}, options::DirectionOptions)
    yI = options.problem.IdealPoint
    v = options.v
    ϵ = options.ϵ
    d = 1 ./ (v - yI .+ ϵ)
    if any(d .< 0)
        d = options.problem.InteriorDirection
    end
    return d / norm(d)
end
function getDirection(::Type{ApproximateAdjDirection}, options::DirectionOptions)
    Vk = options.Vk
    R = options.R
    v = options.v
    vindx = options.vindx
    num_verteces = length(Vk)
    ddim = options.problem.ddim
    Vk = map(Vk) do (vr, i, selected)
        (vr, i, selected, norm(v - vr))
    end |> filter(x -> x[3] != vindx)
    Rk = map(enumerate(eachrow(R))) do (i, r)
        nv = v + r / norm(r)
        (nv, i + num_verteces, false, norm(nv))
    end

    Vd = vcat(Vk, Rk) |> (y -> sort(y, by = x -> x[4]; rev = false))
    Vdlen = length(Vd)
    if Vdlen < ddim
        return options.problem.InteriorDirection
    end
    A = Matrix{Number}(undef, ddim, ddim)
    temp_counter = 1
    temp_set = Set()
    while temp_counter <= Vdlen
        if length(temp_set) == ddim
            break
        end
        nd = Vd[temp_counter][1]
        push!(temp_set, nd)

        temp_counter += 1
    end
    if length(temp_set) < ddim
        return options.problem.InteriorDirection
    end
    for (i, nd) in enumerate(temp_set)
        A[i, :] = nd
    end

    d = A \ ones(ddim)
    thed = if all(d .>= 0)
        d / norm(d)
    else
        options.problem.InteriorDirection
    end
    return thed
end


# struct FixedDirection <: NMDirectoinFunction
#     d::Vector{<:Number}
#     func::Function
# end

# FixedDirection(d::Vector{<:Number}) = FixedDirection(d, () -> d)
# FixedDirection(FD::FixedDirection) =
#     FixedDirection(FD.d, () -> FD.d)


# struct FixedPointDirection <: NMDirectoinFunction
#     p::Vector{<:Number}
#     v::Vector{<:Number}
#     func::Function
# end
# FixedPointDirection(p::Vector{<:Number}, v::Vector{<:Number}) =
#     FixedPointDirection(p, v, () -> begin
#         d = p - v
#         d / norm(d)
#     end)
# FixedPointDirection(FP::FixedPointDirection, v::Vector{<:Number}) =
#     FixedPointDirection(FP.p, v, () -> begin
#         d = FP.p - v
#         d / norm(d)
#     end)


# struct IdealPointDirection <: NMDirectoinFunction
#     yI::Vector{<:Number}
#     v::Vector{<:Number}
#     ϵ::Float64
#     func::Function
# end

# IdealPointDirection(yI::Vector{<:Number}, v::Vector{<:Number}) =
#     IdealPointDirection(yI, v, 1e-5, () -> yI - v)
# IdealPointDirection(IP::IdealPointDirection, v::Vector{<:Number}) =
#     IdealPointDirection(IP.yI, v, IP.ϵ, () -> begin
#         d = 1 ./ (v - IP.yI .+ IP.ϵ)
#         d / norm(d)
#     end)

# struct ApproximateAdjDirection <: NMDirectoinFunction
#     d::Vector{<:Number}
#     ϵ::Float64
#     func::Function
# end
# ApproximateAdjDirection(d::Vector{<:Number}) = ApproximateAdjDirection(d, 1e-5, () -> d)

# ApproximateAdjDirection(AAdj::ApproximateAdjDirection, d::Vector{<:Number}) =
#     ApproximateAdjDirection(d, AAdj.ϵ, () -> d)


# function getDirection(dir_func::T where {T<:NMDirectoinFunction})
#     D = if isa(dir_func, FixedDirection)
#         FixedDirection(dir_func)
#     elseif isa(dir_func, FixedPointDirection)
#         FixedPointDirection(dir_func, v)
#     elseif isa(dir_func, IdealPointDirection)
#         IdealPointDirection(dir_func, v)
#     else
#         num_verteces = length(Vk)
#         ddim = length(v)
#         Vk = map(Vk) do (vr, i, selected)
#             (vr, i, selected, norm(v - vr))
#         end |> filter(x -> x[3] != vindx)
#         Rk = map(enumerate(eachrow(R))) do (i, r)
#             nv = v + r / norm(r)
#             (nv, i + num_verteces, false, norm(nv - v))
#         end

#         Vd = vcat(Vk, Rk) |> (y -> sort(y, by = x -> x[4]))
#         A = Matrix{Number}(undef, ddim, ddim)
#         foreach(1:ddim) do i
#             A[i, :] = Vd[i][1]
#         end

#         d = A \ ones(ddim)
#         thed = if all(d .>= 0)
#             d / norm(d)
#         else

#             dd = ones(ddim)
#             dd / norm(dd)
#         end
#         ApproximateAdjDirection(dir_func, thed)
#     end
#     return D
# end
# function getDirection(
#     dir_func::T where {T<:NMDirectoinFunction},
#     X::Union{Vector{<:Number},Matrix{<:Number}},
#     f::Function,
# )
#     D = if isa(dir_func, FixedPointDirection)
#         n = length(dir_func.v)
#         p = Vector{Float64}(undef, n)
#         ismatrix = isa(X, Matrix)
#         p = map(1:n) do i
#             ismatrix ? maximum(f(X[:, i])) : maximum(f(X[i]))
#         end
#         FixedPointDirection(2p, dir_func.v)
#     else
#         dir_func
#     end
#     return D
# end
export FixedDirection,
    FixedPointDirection, IdealPointDirection, ApproximateAdjDirection, symbolToDirection
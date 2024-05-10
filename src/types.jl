abstract type NMScalar end
abstract type NMVector end
abstract type NMDirection end
abstract type NMVertex end


struct NMCone{T}
    A::Matrix{T}
end

struct MOAProblem
    pdim::Int
    ddim::Int
    C::NMCone
    Cp::NMCone
    f::Function
    WS::Function
    PS::Function
    DPS::Function
    X::Matrix{<:Number}
    InteriorDirection::Vector{<:Number}
    IdealPoint::Union{Nothing,Vector{<:Number}}
    InnerPoint::Vector{<:Number}
end
MOAProblem(
    pdim::Int,
    ddim::Int,
    C::NMCone,
    Cp::NMCone,
    f::Function,
    WS::Function,
    PS::Function,
    DPS::Function,
    InteriorDirection::Vector{<:Number},
    IdealPoint::Union{Nothing,Vector{<:Number}},
) = begin
    W = Cp.A
    n = pdim
    m = ddim
    X = Matrix{Float64}(undef, n, m)
    if n == 1
        for i = 1:m
            X[1, i] = WS(W[:, i])
        end
    else
        for i = 1:m
            X[:, i] = WS(W[:, i])
        end
    end
    p = Vector{Float64}(undef, m)
    isonedim = n == 1
    p = map(1:m) do i
        isonedim ? maximum(f(X[1, i])) : maximum(f(X[:, i]))
    end
    MOAProblem(pdim, ddim, C, Cp, f, WS, PS, DPS, X, InteriorDirection, IdealPoint, 2p)
end

struct MOAOptions
    ϵ::Float64
    maxiters::Int
    stopping_criteria::Function
end
MOAOptions(ϵ::Float64, maxiters::Int) = MOAOptions(ϵ, maxiters, (c) -> c >= maxiters)

abstract type MOAResult end
struct Message
    msg::String
end

struct MOAResultName
    direction::String
    vertex::String
end
MOAResultName() = MOAResultName("", "")
struct MOASuccessResult <: MOAResult
    name::MOAResultName
    message::String
    X::Matrix{<:Number}
    Vk::Vector
    Pk::Polyhedron
    SC::Number
    Card::Number
    T::Number
    HD::Union{Nothing,Number}
    problem::MOAProblem
end
Base.show(io::IO, R::MOASuccessResult) = begin
    header = [
        "Direction"
        "Vertex"
        "n"
        "m"
        "Number of solutions"
        "Cardinality"
        "Number of vertices"
        "Number of Scalarization"
        "CPU Time"
        "Notes"
    ]
    Tbl =
        [R.name.direction R.name.vertex R.problem.pdim R.problem.ddim size(R.X, 2) R.Card length(
            R.Vk,
        ) R.SC R.T R.message]

    pretty_table(io, Tbl; header = header)
end


MOASuccessResult(
    name::MOAResultName,
    msg::Message,
    X::Matrix{<:Number},
    Vk::Vector,
    Pk::Polyhedron,
    SC::Number,
    Card::Number,
    p::MOAProblem,
) = MOASuccessResult(name, msg.msg, X, Vk, Pk, SC, Card, 0.0, nothing, p)
MOASuccessResult(R::MOASuccessResult, T::Float64) = MOASuccessResult(
    R.name,
    R.message,
    R.X,
    R.Vk,
    R.Pk,
    R.SC,
    R.Card,
    T,
    nothing,
    R.problem,
)
MOASuccessResult(R::MOASuccessResult, T::Float64, HD::Number) =
    MOASuccessResult(R.name, R.message, R.X, R.Vk, R.Pk, R.SC, R.Card, T, HD, R.problem)
MOASuccessResult(R::MOASuccessResult, direction::String, vertex::String, T::Float64) =
    MOASuccessResult(
        MOAResultName(direction, vertex),
        R.message,
        R.X,
        R.Vk,
        R.Pk,
        R.SC,
        R.Card,
        T,
        nothing,
        R.problem,
    )

struct MOAFailureResult <: MOAResult
    message::String
    T::Number
end
MOAFailureResult(M::String) = MOAFailureResult(M, 0.0)
MOAFailureResult(R::MOAFailureResult, T::Number) = MOAFailureResult(R.message, T)
export NMDirection,
    NMVertex,
    NMCone,
    MOAProblem,
    MOAOptions,
    MOASuccessResult,
    MOAFailureResult,
    MOAResult,
    MOAResultName
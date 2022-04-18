
struct OffsetLinearConstraint{S,P,W,T} <: StageConstraint
    n::Int
    m::Int
    A::SizedMatrix{P,W,T,2,Matrix{T}}
    b::SVector{P,T}
    sense::S
    l::Float64
    inds::SVector{3,Int}
    function OffsetLinearConstraint(n::Int, m::Int, A::StaticMatrix{P,W,T},
        b::StaticVector{P,T},
        sense::ConstraintSense,
        l::Float64, inds=1:3) where {P,W,T}
        @assert size(A, 1) == length(b) "Length of A, b must be equal"
        @assert W == 2 "Only support two dimension now"
        inds = SVector{3}(inds)
        new{typeof(sense),P,W,T}(n, m, A, b, sense, l, inds)
    end
end

function OffsetLinearConstraint(n::Int, m::Int, A::AbstractMatrix, b::AbstractVector,
    sense::S, l::Float64, inds=1:3) where {S<:ConstraintSense}
    @assert size(A, 1) == length(b) "Length of A, b must be equal"
    @assert size(A, 2) == 2 "Currently only support 2 dimensions"
    p, q = size(A)
    T = promote_type(eltype(A), eltype(b))
    A = SizedMatrix{p,q,T}(A)
    b = SVector{p,T}(b)
    OffsetLinearConstraint(n, m, A, b, sense, l, inds)
end

Base.copy(con::OffsetLinearConstraint) where S =
    OffsetLinearConstraint(con.n, con.m, copy(con.A), copy(con.b), S(), con.l, con.inds)

@inline TO.sense(con::OffsetLinearConstraint) = con.sense
@inline RD.output_dim(::OffsetLinearConstraint{<:Any,P}) where {P} = P
@inline RD.state_dim(con::OffsetLinearConstraint) = con.n
@inline RD.control_dim(con::OffsetLinearConstraint) = con.m
RD.functioninputs(::OffsetLinearConstraint) = RD.StateOnly()

function RD.evaluate(con::OffsetLinearConstraint{S,P}, X::RD.DataVector) where {S,P}
    A = con.A
    b = con.b
    l = con.l
    x, y, θ = X[con.inds]
    x0 = x + l * cos(θ)
    y0 = y + l * sin(θ)
    c = zeros(P)
    for i in 1:P
        c[i] = dot(A[i, 1:2], (x0, y0))
    end
    c .-= b
end

function RD.evaluate!(con::OffsetLinearConstraint{<:Any,P}, c, X::RD.DataVector) where {P}
    A = con.A
    b = con.b
    l = con.l
    x, y, θ = X[con.inds]
    x0 = x + l * cos(θ)
    y0 = y + l * sin(θ)
    for i in 1:P
        c[i] = dot(A[i, 1:2], (x0, y0))
    end
    c .-= b
    return nothing
end

function RD.jacobian!(con::OffsetLinearConstraint{<:Any,P}, ∇c, c, z::RD.AbstractKnotPoint) where {P}
    X = RD.state(z)
    A = con.A
    l = con.l
    xi, yi, θi = con.inds
    θ = X[θi]
    for i = 1:P
        ∇c[i, xi] = A[i, 1]
        ∇c[i, yi] = A[i, 2]
        ∇c[i, θi] = l * dot(A[i, 1:2], (-sin(θ), cos(θ)))
    end
    return nothing
end

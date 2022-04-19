struct OffsetEllipseConstraint{P,T} <: StageConstraint
  n::Int
  m::Int
  x::SVector{P,T}
  y::SVector{P,T}
  a::SVector{P,T}
  b::SVector{P,T}
  ψ::SVector{P,T}
  l::Float64
  inds::SVector{3,Int}
  function OffsetEllipseConstraint{P,T}(n::Int, m::Int, xc::AbstractVector, yc::AbstractVector,
    a::AbstractVector, b::AbstractVector, ψ::AbstractVector,
    l::Float64, inds=1:3) where {P,T}
    @assert length(xc) == length(yc) == length(a) == length(b) == length(ψ)
    "Length of xc, yc, a, b, ψ should be equal"
    new{P,T}(n, m, xc, yc, a, b, ψ, l, inds)
  end
end

function OffsetEllipseConstraint(n::Int, m::Int, xc::AbstractVector, yc::AbstractVector,
  a::AbstractVector, b::AbstractVector, ψ::AbstractVector, l::Float64, inds=1:3)
  T = promote_type(eltype(xc), eltype(yc), eltype(ψ))
  P = length(xc)
  OffsetEllipseConstraint{P,T}(n, m, xc, yc, a, b, ψ, l, inds)
end

@inline TO.sense(::OffsetEllipseConstraint) = Inequality()
@inline RD.state_dim(con::OffsetEllipseConstraint) = con.n
@inline RD.control_dim(con::OffsetEllipseConstraint) = con.m
@inline RD.output_dim(::OffsetEllipseConstraint{P}) where {P} = P
@inline RD.functioninputs(::OffsetEllipseConstraint) = RD.StateOnly()

function RD.evaluate(con::OffsetEllipseConstraint{P}, X::RD.DataVector) where {P}
  xc = con.x
  yc = con.y
  a = con.a
  b = con.b
  ψ = con.ψ
  c = zeros(P)
  xi, yi, θi = con.inds
  x = X[xi]
  y = X[yi]
  θ = X[θi]
  Δx = con.l * cos(θ)
  Δy = con.l * sin(θ)
  for i in 1:P
    sn, cn = sincos(ψ[i])
    h = xc[i]
    k = yc[i]
    x′ = cn * x + sn * y - cn * h - sn * k + cn * Δx + sn * Δy
    y′ = -sn * x + cn * y + sn * h - cn * k - sn * Δx + cn * Δy
    c[i] = -(x′ / a[i])^2 - (y′ / b[i])^2 + 1
  end
  c
end

function RD.evaluate!(con::OffsetEllipseConstraint{P}, c, X::RD.DataVector) where {P}
  xc = con.x
  yc = con.y
  a = con.a
  b = con.b
  ψ = con.ψ
  xi, yi, θi = con.inds
  x = X[xi]
  y = X[yi]
  θ = X[θi]
  Δx = con.l * cos(θ)
  Δy = con.l * sin(θ)
  for i in 1:P
    sn, cn = sincos(ψ[i])
    h = xc[i]
    k = yc[i]
    x′ = cn * x + sn * y - cn * h - sn * k + cn * Δx + sn * Δy
    y′ = -sn * x + cn * y + sn * h - cn * k - sn * Δx + cn * Δy
    c[i] = -(x′ / a[i])^2 - (y′ / b[i])^2 + 1
  end
  return nothing
end

function RD.jacobian!(con::OffsetEllipseConstraint{P}, ∇c, c, z::RD.AbstractKnotPoint) where {P}
  X = RD.state(z)
  xc = con.x
  yc = con.y
  a = con.a
  b = con.b
  ψ = con.ψ
  xi, yi, θi = con.inds
  x = X[xi]
  y = X[yi]
  θ = X[θi]
  Δx = con.l * cos(θ)
  Δy = con.l * sin(θ)
  for i = 1:P
    sn, cn = sincos(ψ[i])
    h = xc[i]
    k = yc[i]
    x′ = cn * x + sn * y - cn * h - sn * k + cn * Δx + sn * Δy
    y′ = -sn * x + cn * y + sn * h - cn * k - sn * Δx + cn * Δy
    dxx = cn
    dxy = sn
    dyx = -sn
    dyy = cn
    ∇c[i, xi] = -2 * (x′ / a[i]) * dxx - 2 * (y′ / b[i]) * dyx
    ∇c[i, yi] = -2 * (x′ / a[i]) * dxy - 2 * (y′ / b[i]) * dyy
    ∇c[i, θi] = -2 * (x′ / a[i]) * (cn * Δy + sn * Δx) -
                2 * (y′ / b[i]) * (sn * Δy + cn * Δx)
  end
  return nothing
end

################################################
struct EllipseConstraint{P,T} <: StageConstraint
  n::Int
  m::Int
  x::SVector{P,T}
  y::SVector{P,T}
  a::SVector{P,T}
  b::SVector{P,T}
  ψ::SVector{P,T}
  inds::SVector{2,Int}
  function EllipseConstraint{P,T}(n::Int, m::Int, xc::AbstractVector, yc::AbstractVector,
    a::AbstractVector, b::AbstractVector, ψ::AbstractVector,
    inds=1:2) where {P,T}
    @assert length(xc) == length(yc) == length(a) == length(b) == length(ψ)
    "Length of xc, yc, a, b, ψ should be equal"
    new{P,T}(n, m, xc, yc, a, b, ψ, inds)
  end
end

function EllipseConstraint(n::Int, m::Int, xc::AbstractVector, yc::AbstractVector,
  a::AbstractVector, b::AbstractVector, ψ::AbstractVector, inds=1:2)
  T = promote_type(eltype(xc), eltype(yc), eltype(ψ))
  P = length(xc)
  EllipseConstraint{P,T}(n, m, xc, yc, a, b, ψ, inds)
end

@inline TO.sense(::EllipseConstraint) = Inequality()
@inline RD.state_dim(con::EllipseConstraint) = con.n
@inline RD.control_dim(con::EllipseConstraint) = con.m
@inline RD.output_dim(::EllipseConstraint{P}) where {P} = P
@inline RD.functioninputs(::EllipseConstraint) = RD.StateOnly()

function RD.evaluate(con::EllipseConstraint{P}, X::RD.DataVector) where {P}
  xc = con.x
  yc = con.y
  a = con.a
  b = con.b
  ψ = con.ψ
  c = zeros(P)
  xi, yi = con.inds
  x = X[xi]
  y = X[yi]
  for i in 1:P
    sn, cn = sincos(ψ[i])
    h = xc[i]
    k = yc[i]
    x′ = cn * x + sn * y - cn * h - sn * k
    y′ = -sn * x + cn * y + sn * h - cn * k
    c[i] = -(x′ / a[i])^2 - (y′ / b[i])^2 + 1
  end
  c
end

function RD.evaluate!(con::EllipseConstraint{P}, c, X::RD.DataVector) where {P}
  xc = con.x
  yc = con.y
  a = con.a
  b = con.b
  ψ = con.ψ
  xi, yi = con.inds
  x = X[xi]
  y = X[yi]
  for i in 1:P
    sn, cn = sincos(ψ[i])
    h = xc[i]
    k = yc[i]
    x′ = cn * x + sn * y - cn * h - sn * k
    y′ = -sn * x + cn * y + sn * h - cn * k
    c[i] = -(x′ / a[i])^2 - (y′ / b[i])^2 + 1
  end
  return nothing
end

function RD.jacobian!(con::EllipseConstraint{P}, ∇c, c, z::RD.AbstractKnotPoint) where {P}
  X = RD.state(z)
  xc = con.x
  yc = con.y
  a = con.a
  b = con.b
  ψ = con.ψ
  xi, yi = con.inds
  x = X[xi]
  y = X[yi]
  for i = 1:P
    sn, cn = sincos(ψ[i])
    h = xc[i]
    k = yc[i]
    x′ = cn * x + sn * y - cn * h - sn * k
    y′ = -sn * x + cn * y + sn * h - cn * k
    dxx = cn
    dxy = sn
    dyx = -sn
    dyy = cn
    ∇c[i, xi] = -2 * (x′ / a[i]) * dxx - 2 * (y′ / b[i]) * dyx
    ∇c[i, yi] = -2 * (x′ / a[i]) * dxy - 2 * (y′ / b[i]) * dyy
  end
  return nothing
end

#####################################################
struct OffsetCircleConstraint{P,T} <: StageConstraint
  n::Int
  m::Int
  x::SVector{P,T}
  y::SVector{P,T}
  radius::SVector{P,T}
  l::Float64
  inds::SVector{3,Int}
  function OffsetCircleConstraint{P,T}(n::Int, m::Int, xc::AbstractVector, yc::AbstractVector,
    radius::AbstractVector, l::Float64, inds=1:3) where {P,T}
    @assert length(xc) == length(yc) == length(radius) == P "Length must be equal"
    new{P,T}(n, m, xc, yc, radius, l, inds)
  end
end

function OffsetCircleConstraint(n::Int, m::Int, xc::AbstractVector, yc::AbstractVector, radius::AbstractVector,
  l::Float64, inds=1:3)
  T = promote_type(eltype(xc), eltype(yc), eltype(radius))
  P = length(xc)
  OffsetCircleConstraint{P,T}(n, m, xc, yc, radius, l, inds)
end

@inline TO.sense(::OffsetCircleConstraint) = Inequality()
@inline RD.state_dim(con::OffsetCircleConstraint) = con.n
@inline RD.control_dim(con::OffsetCircleConstraint) = con.m
@inline RD.output_dim(::OffsetCircleConstraint{P}) where {P} = P
@inline RD.functioninputs(::OffsetCircleConstraint) = RD.StateOnly()

function RD.evaluate(con::OffsetCircleConstraint, X::RD.DataVector)
  xc = con.x
  yc = con.y
  r = con.radius
  xi, yi, θi = con.inds
  x = X[xi]
  y = X[yi]
  θ = X[θi]
  l = con.l
  -(x + l * cos(θ) .- xc) .^ 2 - (y + l * sin(θ) .- yc) .^ 2 + r .^ 2
end

function RD.evaluate!(con::OffsetCircleConstraint{P}, c, X::RD.DataVector) where {P}
  xc = con.x
  yc = con.y
  r = con.radius
  xi, yi, θi = con.inds
  x = X[xi]
  y = X[yi]
  θ = X[θi]
  l = con.l
  for i = 1:P
    c[i] = -(x + l * cos(θ) - xc[i])^2 - (y + l * sin(θ) - yc[i])^2 + r[i]^2
  end
  return nothing
end

function RD.jacobian!(con::OffsetCircleConstraint{P}, ∇c, c, z::RD.AbstractKnotPoint) where {P}
  X = RD.state(z)
  xc = con.x
  yc = con.y
  xi, yi, θi = con.inds
  x = X[xi]
  y = X[yi]
  θ = X[θi]
  l = con.l
  s, c = sincos(θ)
  for i = 1:P
    ∇c[i, xi] = -2 * (x + l * c - xc[i])
    ∇c[i, yi] = -2 * (y + l * s - yc[i])
    ∇c[i, θi] = l * dot((∇c[i, xi], ∇c[i, yi]), (-s, c))
  end
  return nothing
end

#########################################################
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

Base.copy(con::OffsetLinearConstraint) where {S} =
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
  s, c = sincos(θ)
  for i = 1:P
    ∇c[i, xi] = A[i, 1]
    ∇c[i, yi] = A[i, 2]
    ∇c[i, θi] = l * dot(A[i, 1:2], (-s, c))
  end
  return nothing
end

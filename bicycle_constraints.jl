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
  function OffsetEllipseConstraint{P,T}(
    n::Int,
    m::Int,
    xc::AbstractVector,
    yc::AbstractVector,
    a::AbstractVector,
    b::AbstractVector,
    ψ::AbstractVector,
    l::Float64,
    inds = 1:3,
  ) where {P,T}
    @assert length(xc) == length(yc) == length(a) == length(b) == length(ψ)
    "Length of xc, yc, a, b, ψ should be equal"
    new{P,T}(n, m, xc, yc, a, b, ψ, l, inds)
  end
end

function OffsetEllipseConstraint(
  n::Int,
  m::Int,
  xc::AbstractVector,
  yc::AbstractVector,
  a::AbstractVector,
  b::AbstractVector,
  ψ::AbstractVector,
  l::Float64,
  inds = 1:3,
)
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
  for i = 1:P
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
  for i = 1:P
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
    ∇c[i, θi] = -2 * (x′ / a[i]) * (cn * Δy + sn * Δx) - 2 * (y′ / b[i]) * (sn * Δy + cn * Δx)
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
  function EllipseConstraint{P,T}(
    n::Int,
    m::Int,
    xc::AbstractVector,
    yc::AbstractVector,
    a::AbstractVector,
    b::AbstractVector,
    ψ::AbstractVector,
    inds = 1:2,
  ) where {P,T}
    @assert length(xc) == length(yc) == length(a) == length(b) == length(ψ)
    "Length of xc, yc, a, b, ψ should be equal"
    new{P,T}(n, m, xc, yc, a, b, ψ, inds)
  end
end

function EllipseConstraint(
  n::Int,
  m::Int,
  xc::AbstractVector,
  yc::AbstractVector,
  a::AbstractVector,
  b::AbstractVector,
  ψ::AbstractVector,
  inds = 1:2,
)
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
  for i = 1:P
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
  for i = 1:P
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
  function OffsetCircleConstraint{P,T}(
    n::Int,
    m::Int,
    xc::AbstractVector,
    yc::AbstractVector,
    radius::AbstractVector,
    l::Float64,
    inds = 1:3,
  ) where {P,T}
    @assert length(xc) == length(yc) == length(radius) == P "Length must be equal"
    new{P,T}(n, m, xc, yc, radius, l, inds)
  end
end

function OffsetCircleConstraint(
  n::Int,
  m::Int,
  xc::AbstractVector,
  yc::AbstractVector,
  radius::AbstractVector,
  l::Float64,
  inds = 1:3,
)
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
struct OffsetLineSegmentConstraint{P,W,T} <: StageConstraint
  n::Int
  m::Int
  points::SizedMatrix{P,W,T,2,Matrix{T}}
  l::Float64
  inds::SVector{3,Int}
  dist::Float64
  side::Symbol
  function OffsetLineSegmentConstraint(
    n::Int,
    m::Int,
    points::StaticMatrix{P,W,T},
    l::Float64,
    inds = 1:3;
    dist::Float64 = 1.1,
    side::Symbol = :lower,
  ) where {P,W,T}
    @assert size(points, 2) == W == 5 "Length of segment points coord error"
    inds = SVector{3}(inds)
    new{P,W,T}(n, m, points, l, inds, dist, side)
  end
end

function OffsetLineSegmentConstraint(
  n::Int,
  m::Int,
  points::AbstractMatrix,
  l::Float64,
  inds = 1:3;
  dist = 1.0,
  side = :lower,
)
  @assert size(points, 2) == 5 "Points coord size error"
  @assert side == :lower || side == :upper
  p, q = size(points)
  T = promote_type(eltype(points))
  points = SizedMatrix{p,q,T}(points)
  OffsetLineSegmentConstraint(n, m, points, l, inds, dist = dist, side = side)
end

@inline TO.sense(::OffsetLineSegmentConstraint) = Inequality()
@inline RD.output_dim(::OffsetLineSegmentConstraint{P}) where {P} = P
@inline RD.state_dim(con::OffsetLineSegmentConstraint) = con.n
@inline RD.control_dim(con::OffsetLineSegmentConstraint) = con.m
RD.functioninputs(::OffsetLineSegmentConstraint) = RD.StateOnly()

function RD.evaluate(con::OffsetLineSegmentConstraint, X::RD.DataVector)
  x, y, θ = X[con.inds]
  x0 = x + con.l * cos(θ)
  y0 = y + con.l * sin(θ)
  x1 = con.points[:, 1]
  y1 = con.points[:, 2]
  x2 = con.points[:, 3]
  y2 = con.points[:, 4]
  len = con.points[:, 5]
  if con.side == :lower
    @. ((x2 - x1) * (y0 - y1) - (y2 - y1) * (x0 - x1)) / len + con.dist
  else
    @. ((y2 - y1) * (x0 - x1) - (x2 - x1) * (y0 - y1)) / len + con.dist
  end
end

function RD.evaluate!(con::OffsetLineSegmentConstraint{P}, c, X::RD.DataVector) where {P}
  x, y, θ = X[con.inds]
  x0 = x + con.l * cos(θ)
  y0 = y + con.l * sin(θ)
  for i = 1:P
    x1 = con.points[i, 1]
    y1 = con.points[i, 2]
    x2 = con.points[i, 3]
    y2 = con.points[i, 4]
    len = con.points[i, 5]
    if con.side == :lower
      c[i] = ((x2 - x1) * (y0 - y1) - (y2 - y1) * (x0 - x1)) / len + con.dist
    else
      c[i] = ((y2 - y1) * (x0 - x1) - (x2 - x1) * (y0 - y1)) / len + con.dist
    end
  end
  return nothing
end

function RD.jacobian!(con::OffsetLineSegmentConstraint{P}, ∇c, c, z::RD.AbstractKnotPoint) where {P}
  l = con.l
  θ = z[con.inds[end]]
  x1 = con.points[:, 1]
  y1 = con.points[:, 2]
  x2 = con.points[:, 3]
  y2 = con.points[:, 4]
  len = con.points[:, 5]
  if con.side == :lower
    @. ∇c[:, 1] = (y1 - y2) / len
    @. ∇c[:, 2] = (x2 - x1) / len
    @. ∇c[:, 3] = l * ((x2 - x1) * cos(θ) + (y2 - y1) * sin(θ)) / len
  else
    @. ∇c[:, 1] = (y2 - y1) / len
    @. ∇c[:, 2] = (x1 - x2) / len
    @. ∇c[:, 3] = l * ((x1 - x2) * cos(θ) + (y1 - y2) * sin(θ)) / len
  end

  return nothing
end

#########################################################
struct LineSegmentConstraint{P,W,T} <: StageConstraint
  n::Int
  m::Int
  points::SizedMatrix{P,W,T,2,Matrix{T}}
  dist::Float64
  inds::SVector{2,Int}
  side::Symbol
  function LineSegmentConstraint(
    n::Int,
    m::Int,
    points::StaticMatrix{P,W,T},
    inds = 1:2;
    dist::Float64 = 1.1,
    side::Symbol = :lower,
  ) where {P,W,T}
    @assert size(points, 2) == W == 5 "Length of segment points coord error"
    inds = SVector{2}(inds)
    new{P,W,T}(n, m, points, dist, inds, side)
  end
end

function LineSegmentConstraint(
  n::Int,
  m::Int,
  points::AbstractMatrix,
  inds = 1:2;
  dist = 1.0,
  side = :lower,
)
  @assert size(points, 2) == 5 "Points coord size error"
  @assert side == :lower || side == :upper
  p, q = size(points)
  T = promote_type(eltype(points))
  points = SizedMatrix{p,q,T}(points)
  LineSegmentConstraint(n, m, points, inds, dist = dist, side = side)
end

@inline TO.sense(::LineSegmentConstraint) = Inequality()
@inline RD.output_dim(::LineSegmentConstraint{P}) where {P} = P
@inline RD.state_dim(con::LineSegmentConstraint) = con.n
@inline RD.control_dim(con::LineSegmentConstraint) = con.m
RD.functioninputs(::LineSegmentConstraint) = RD.StateOnly()

function RD.evaluate(con::LineSegmentConstraint, X::RD.DataVector)
  x, y = X[con.inds]
  x1 = con.points[:, 1]
  y1 = con.points[:, 2]
  x2 = con.points[:, 3]
  y2 = con.points[:, 4]
  len = con.points[:, 5]
  if con.side == :lower
    @. ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)) / len + con.dist
  else
    @. ((y2 - y1) * (x - x1) - (x2 - x1) * (y - y1)) / len + con.dist
  end
end

function RD.evaluate!(con::LineSegmentConstraint{P}, c, X::RD.DataVector) where {P}
  x, y = X[con.inds]
  for i = 1:P
    x1 = con.points[i, 1]
    y1 = con.points[i, 2]
    x2 = con.points[i, 3]
    y2 = con.points[i, 4]
    len = con.points[i, 5]
    if con.side == :lower
      c[i] = ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)) / len + con.dist
    else
      c[i] = ((y2 - y1) * (x - x1) - (x2 - x1) * (y - y1)) / len + con.dist
    end
  end
  return nothing
end

function RD.jacobian!(con::LineSegmentConstraint{P}, ∇c, c, z::RD.AbstractKnotPoint) where {P}
  x1 = con.points[:, 1]
  y1 = con.points[:, 2]
  x2 = con.points[:, 3]
  y2 = con.points[:, 4]
  len = con.points[:, 5]
  if con.side == :lower
    @. ∇c[:, 1] = (y1 - y2) / len
    @. ∇c[:, 2] = (x2 - x1) / len
  else
    @. ∇c[:, 1] = (y2 - y1) / len
    @. ∇c[:, 2] = (x1 - x2) / len
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
  function OffsetLinearConstraint(
    n::Int,
    m::Int,
    A::StaticMatrix{P,W,T},
    b::StaticVector{P,T},
    sense::ConstraintSense,
    l::Float64,
    inds = 1:3,
  ) where {P,W,T}
    @assert size(A, 1) == length(b) "Length of A, b must be equal"
    @assert W == 2 "Only support two dimension now"
    inds = SVector{3}(inds)
    new{typeof(sense),P,W,T}(n, m, A, b, sense, l, inds)
  end
end

function OffsetLinearConstraint(
  n::Int,
  m::Int,
  A::AbstractMatrix,
  b::AbstractVector,
  sense::S,
  l::Float64,
  inds = 1:3,
) where {S<:ConstraintSense}
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
  for i = 1:P
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
  for i = 1:P
    c[i] = dot(A[i, 1:2], (x0, y0))
  end
  c .-= b
  return nothing
end

function RD.jacobian!(
  con::OffsetLinearConstraint{<:Any,P},
  ∇c,
  c,
  z::RD.AbstractKnotPoint,
) where {P}
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

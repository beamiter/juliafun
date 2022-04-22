
RD.@autodiff struct MyBicycleModel <: RD.ContinuousDynamics
  ref::Symbol
  L::Float64
  lr::Float64
  function MyBicycleModel(; ref::Symbol = :rear, L::Real = 2.7, lr::Real = 1.5)
    @assert ref ∈ (:cg, :front, :rear)
    @assert L > 0 "($L) must be greater than 0"
    @assert L > lr "($L) must be greater than ($lr)"
    new(ref, L, lr)
  end
end

body = quote
  da = u[1]
  ϕ = u[2]

  θ = x[3]
  δ = x[4]
  v = x[5]
  a = x[6]
  if model.ref == :cg
    β = atan(model.lr * δ, model.L)
    s, c = sincos(θ + β)
    ω = v * cos(β) * tan(δ) / model.L
  elseif model.ref == :rear
    s, c = sincos(θ)
    ω = v * tan(δ) / model.L
  elseif model.ref == :front
    s, c = sincos(θ + δ)
    ω = v * sin(δ) / model.L
  end
  ẋ = v * c
  ẏ = v * s
end

@eval function RD.dynamics(model::MyBicycleModel, x, u)
  $body
  return SA[ẋ, ẏ, ω, ϕ, a, da]
end
@eval function RD.dynamics!(model::MyBicycleModel, xdot, x, u)
  $body
  xdot[1] = ẋ
  xdot[2] = ẏ
  xdot[3] = ω
  xdot[4] = ϕ
  xdot[5] = a
  xdot[6] = da
  return nothing
end
RD.state_dim(::MyBicycleModel) = 6
RD.control_dim(::MyBicycleModel) = 2

function BicycleCar(x0, xf, N, tf, ; ref = :cg, U0 = SA_F64[0, 0])
  # 1. define model dynamics
  model = MyBicycleModel(; ref = ref)

  # 2. define solver options
  opts = SolverOptions(penalty_initial = 1e4, verbose = 0, cost_tolerance_intermediate = 1e-1)

  # 3. define Q, R matrix
  # x, y, theta, delta, v, a
  Q = Diagonal(SA_F64[10, 10, 60, 1, 1, 1])
  # jerk, phi
  ρ = 1.0
  R = ρ * Diagonal(SA_F64[1, 1])
  Qf = Diagonal(SA_F64[10, 10, 60, 1, 1, 1])
  obj = LQRObjective(Q, R, Qf, xf, N)


  # 4. define problem
  prob = Problem(model, obj, x0, tf, xf = xf)

  # 5. init, can use warm start here
  initial_controls!(prob, U0)
  rollout!(prob)

  return prob, opts
end

function BicycleCar(
  x0,
  xf,
  N,
  tf,
  cons::ConstraintList,
  ;
  ref = :cg,
  X0 = SA_F64[0, 0, 0, 0, 0, 0],
  U0 = SA_F64[0, 0],
)
  # 1. define model dynamics
  model = MyBicycleModel(; ref = ref)

  # 2. define solver options
  opts = SolverOptions(penalty_initial = 1e4, verbose = 0, cost_tolerance_intermediate = 1e-1)

  # 3. define Q, R matrix
  # x, y, theta, delta, v, a
  Q = Diagonal(SA_F64[10, 10, 60, 1, 1, 1])
  # jerk, phi
  ρ = 1.0
  R = ρ * Diagonal(SA_F64[1, 1])
  Qf = Diagonal(SA_F64[10, 10, 60, 1, 1, 1])
  obj = LQRObjective(Q, R, Qf, xf, N)


  # 4. define problem
  prob = Problem(model, obj, x0, tf, xf = xf, constraints = cons)

  # 5. init, can use warm start here
  initial_controls!(prob, U0)
  initial_states!(prob, X0)
  # rollout!(prob)

  return prob, opts
end

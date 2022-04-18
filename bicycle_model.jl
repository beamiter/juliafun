
RD.@autodiff struct MyBicycleModel <: RD.ContinuousDynamics
    ref::Symbol
    L::Float64
    lr::Float64
    function MyBicycleModel(; ref::Symbol=:rear, L::Real=2.7, lr::Real=1.5)
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

function BicycleCar(scenario=:parallel_park, ; N=51, x0)
    if scenario == :parallel_park
    end
    model = MyBicycleModel(; ref=:cg)
    n, m = size(model)

    opts = SolverOptions(penalty_initial=1e4, verbose=0,
        cost_tolerance_intermediate=1e-1)

    tf = 5.0
    #dt = tf / (N - 1)

    xf = SA[13, -2.2, deg2rad(0), 0, 0.1, 0]

    # x, y, theta, delta, v, a
    Q = Diagonal(SA_F64[10, 10, 60, 1, 1, 1])
    # jerk, phi
    ρ = 1.0
    R = ρ * Diagonal(SA_F64[1, 1])
    Qf = Diagonal(SA_F64[10, 10, 60, 1, 1, 1])
    obj = LQRObjective(Q, R, Qf, xf, N)

    cons = ConstraintList(n, m, N)
    bnd_x_l = [-1, -2.4, Inf, -deg2rad(45), 0.0, -3]
    bnd_x_u = [20, 2.4, Inf, deg2rad(45), 6.0, 3]
    bnd_u_l = [-3, -deg2rad(45)]
    bnd_u_u = [3, deg2rad(45)]

    bnd = BoundConstraint(n, m, x_min=bnd_x_l, x_max=bnd_x_u, u_min=bnd_u_l,
        u_max=bnd_u_u)

    p = 2
    A = zeros(p, m + n)
    A[1, 5] = 1.0
    A[2, 5] = -1.0
    A = SMatrix{p,m + n,Float64}(A)
    b = zeros(p)
    b[1] = 4.0
    b[2] = -0.01
    b = SVector{p,Float64}(b)
    lin = LinearConstraint(n, m, A, b, Inequality())

    xc = SA[7]
    yc = SA[0.5]
    r = SA[1.0]
    cir = CircleConstraint(n, xc, yc, r)

    goal = GoalConstraint(xf)

    p = 1
    A = zeros(p, 2)
    b = zeros(p)
    A[2] = -1.0
    b[1] = 2.4
    off = OffsetLinearConstraint(n, m, A, b, Inequality(), 1.0)

    add_constraint!(cons, bnd, 1:N-1)
    # add_constraint!(cons, goal, N)
    add_constraint!(cons, lin, 2:N)
    # add_constraint!(cons, cir, 1:N)
    add_constraint!(cons, off, 1:N)

    prob = Problem(model, obj, x0, tf, xf=xf, constraints=cons)
    initial_controls!(prob, SA[0.0, 0.0])
    rollout!(prob)


    return prob, opts
end

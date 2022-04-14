using Altro
using RobotZoo
using StaticArrays
using LinearAlgebra
using TrajectoryOptimization
const TO = TrajectoryOptimization
using BenchmarkTools
using Plots
using ForwardDiff, FiniteDiff
using RobotDynamics
const RD = RobotDynamics
import TrajectoryOptimization: StageConstraint, ConstraintSense
using Plots: plot, plot!

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

function loop_for_gif()
    x0 = SA_F64[0, 0, 0, 0, 4, 0]
    plt = plot([0, 20, 20, 0, 0], 2.0 * [-1.2, -1.2, 1.2, 1.2, -1.2])
    anim = @animate for i in UnitRange(1, 40)
        @show x0
        bicycle = BicycleCar(:parallel_park, x0=x0)
        solver = ALTROSolver(bicycle...)
        solve!(solver)
        X = states(solver)
        U = controls(solver)
        @show size(X)
        @show size(U)
        x0 = X[2]
        p = plot(plt, [x[1] for x in X], [x[2] for x in X])
        display(p)
        sleep(5)
    end
    gif(anim, "anim_fps15.gif", fps=1)
end
function loop_for_display()
    gr()
    x0 = SA_F64[0, 0, 0, 0, 4, 0]
    plt = plot([0, 20], [-2.4, -2.4])
    plot!(plt, [0, 20], [-2.0, -2.0])
    plot!(plt, [0, 20], [2.4, 2.4])
    his_x = []
    his_y = []
    his_x_f = []
    his_y_f = []
    for i in UnitRange(1, 40)
        @show x0
        bicycle = BicycleCar(:parallel_park, x0=x0)
        solver = ALTROSolver(bicycle...)
        solve!(solver)
        X = states(solver)
        U = controls(solver)
        @show size(X)
        @show size(U)
        push!(his_x, x0[1])
        push!(his_y, x0[2])
        l = 1.0
        push!(his_x_f, x0[1] + l * cos(x0[3]))
        push!(his_y_f, x0[2] + l * sin(x0[3]))
        x0 = X[2]
        p = plot(plt, [x[1] for x in X], [x[2] for x in X])
        r = 1.0
        theta = (-pi:0.1:pi)
        plot!(p, [r * cos(i) + 7.0 for i in theta], [r * sin(i) + 0.5 for i in theta])
        r = 1.0
        plot!(p, [r * cos(i) + his_x[end] for i in theta], [r * sin(i) + his_y[end] for i in theta])
        plot!(p, [r * cos(i) + his_x_f[end] for i in theta], [r * sin(i) + his_y_f[end] for i in theta])
        scatter!(p, his_x, his_y)
        scatter!(p, his_x_f, his_y_f)
        display(p)
        sleep(0.5)
        readline()
    end
end

#loop_for_display()

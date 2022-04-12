using TrajectoryOptimization
using Altro
import RobotZoo.Cartpole
using StaticArrays, LinearAlgebra
using RobotDynamics

model = Cartpole()
n, m = size(model)

N = 101
tf = 5.0
dt = tf / (N - 1)

x0 = @SVector zeros(n)
xf = @SVector [0, pi, 0, 0]

Q = 1.0e-2 * Diagonal(@SVector ones(n))
Qf = 100.0 * Diagonal(@SVector ones(n))
R = 1.0e-1 * Diagonal(@SVector ones(m))

obj = LQRObjective(Q, R, Qf, xf, N)

conSet = ConstraintList(n, m, N)
u_bnd = 3.0
bnd = BoundConstraint(n, m, u_min=-u_bnd, u_max=u_bnd)
goal = GoalConstraint(xf)
add_constraint!(conSet, bnd, 1:N-1)
add_constraint!(conSet, goal, N)

u0 = @SVector fill(0.01, m)
U0 = [u0 for k = 1:N-1]

prob = Problem(model, obj, xf, tf, x0=x0, constraints=conSet)
initial_controls!(prob, U0)

opts = SolverOptions(
  cost_tolerance_intermediate=1e-2,
  penalty_scaling=10.0,
  penalty_initial=1.0
)

altro = ALTROSolver(prob, opts)
solve!(altro)

println(max_violation(altro), ",", cost(altro), ",", iterations(altro))

X = states(altro)
U = controls(altro)

stats = Altro.stats(altro)
println(stats.iterations, ",", stats.iterations_outer, ",", stats.iterations_pn, ",",
  stats.cost[end], ",", stats.c_max[end], ",", stats.gradient[end])

dstats = Dict(stats)


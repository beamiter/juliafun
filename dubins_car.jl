using RobotZoo
using Altro
using StaticArrays
using LinearAlgebra
using TrajectoryOptimization
using RobotDynamics

function DubinsCar(scenario=:three_obstacles; N=101)
  opts = SolverOptions(
    cost_tolerance_intermediate=1e-2,
    penalty_scaling=10.0,
    penalty_initial=10.0,
  )

  model = RobotZoo.DubinsCar()
  n, m = size(model)

  N = 101
  tf = 5.0
  dt = tf / (N - 1)
  x0 = @SVector [0.0, 0.0, 0.0]
  xf = @SVector [3.0, 3.0, 0.0]

  Q = Diagonal(@SVector [1.0, 1.0, 1.0])
  R = Diagonal(@SVector [0.5, 0.5])
  Qf = 10.0 * Diagonal(@SVector ones(n))
  obj = LQRObjective(Q, R, Qf, xf, N)

  r_circle_3obs = 0.25

  circle_x = 3 * @SVector [0.25, 0.5, 0.75]
  circle_y = 3 * @SVector [0.25, 0.5, 0.75]
  circle_r = @SVector fill(r_circle_3obs + model.radius, 3)

  obs = CircleConstraint(n, circle_x, circle_y, circle_r)
  bnd = BoundConstraint(n, m, u_min=[0, -3], u_max=[3, 3])
  goal = GoalConstraint(xf)

  conSet = ConstraintList(n, m, N)
  add_constraint!(conSet, obs, 2:N-1)
  add_constraint!(conSet, bnd, 1:N-1)
  add_constraint!(conSet, goal, N:N)

  U = [@SVector fill(0.01, m) for _ = 1:N-1]
  car_3obs_static = Problem(model, obj, xf, tf, constraints=conSet, x0=x0)
  initial_controls!(car_3obs_static, U)
  rollout!(car_3obs_static)
  return car_3obs_static, opts
end


solver = ALTROSolver(DubinsCar(:three_obstacles)...)
solve!(solver)
# b = benchmark_solve!(solver)
# @show minimum(b).time/1e6


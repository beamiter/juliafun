using BenchmarkTools
using Altro
using TrajectoryOptimization
using Test
using StaticArrays
const TO = TrajectoryOptimization

solver = ALTROSolver(Problems.DubinsCar(:parallel_park)..., show_summary=true)
b = benchmark_solve!(solver)
@test minimum(b).time / 1e6 < 8
@test max_violation(solver) < 1e-6
@test iterations(solver) == 13
@test solver.stats.gradient[end] < 1e-3
@test status(solver) == Altro.SOLVE_SUCCEEDED
@show median(b), mean(b)

solver = ALTROSolver(Problems.DubinsCar(:three_obstacles)...)
b = benchmark_solve!(solver)
@test minimum(b).time /1e6 < 6 
@test max_violation(solver) < 1e-6
@test iterations(solver) == 20 # 20
@test solver.stats.gradient[end] < 1e-2
@test status(solver) == Altro.SOLVE_SUCCEEDED 


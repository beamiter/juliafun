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

include("bicycle_model.jl")
include("bicycle_constraints.jl")
include("constraints_sampler.jl")
include("bicycle_sim.jl")

## CairoMakie demo
# using CairoMakie
# x = range(0, 10, length=100)
# y1 = sin.(x)
# y2 = cos.(x)
# f, = lines(x, y1)
# lines!(x, y2)
# save("figure.png", f)

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
include("bicycle_sim.jl")

using RobotZoo, RobotDynamics
# const RD = RobotDynamics

using TrajectoryOptimization
# const TO = TrajectoryOptimization

using Altro
using Test

using StaticArrays, LinearAlgebra

import RobotZoo.Cartpole

include("bicycle.jl")

model = Cartpole()
n, m = RD.dims(model)
x, u = rand(model)
t, h = 1.1, 0.1
z = KnotPoint(x, u, t, h)

#################################################

p = 5
A = @SMatrix rand(p, n+m)
b = @SVector rand(p)
lin = LinearConstraint(n, m, A, b, Inequality())
@show RD.evaluate(lin, z)
@show RD.output_dim(lin)

#################################################

xc = SA[1,1,1]
yc = SA[1,2,3]
r = SA[1,1,1]
cir = CircleConstraint(n, xc, yc, r)
c = zeros(3)
@show RD.evaluate(cir,z)
@show RD.output_dim(cir)

#################################################

p = 5
A = @SVector rand(p)
b = @SVector rand(p)
lx = 0.1
ly = 0.1
off = OffsetLinearConstraint(n, m, A, b, Inequality(), lx, ly)
@show RD.evaluate(off, z)
@show RD.output_dim(off)
c = zeros(RD.output_dim(off))
RD.evaluate!(off, c, z)
@test c ≈ RD.evaluate(off, z)
∇c = zeros(2 * p, n)
RD.jacobian!(off, ∇c, c, z)
@show ∇c

#################################################

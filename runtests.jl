using Test
import RobotZoo.Cartpole

include("bicycle.jl")

@testset "unit tests" begin

  model = Cartpole()
  n, m = RD.dims(model)
  x, u = rand(model)
  t, h = 1.1, 0.1
  z = KnotPoint(x, u, t, h)

  #################################################

  p = 5
  A = @SMatrix rand(p, n + m)
  b = @SVector rand(p)
  lin = LinearConstraint(n, m, A, b, Inequality())
  @show RD.evaluate(lin, z)
  @show RD.output_dim(lin)
  @show typeof(lin)

  #################################################

  xc = SA[1, 1, 1]
  yc = SA[1, 2, 3]
  r = SA[1, 1, 1]
  cir = CircleConstraint(n, xc, yc, r)
  c = zeros(3)
  @show RD.evaluate(cir, z)
  @show RD.output_dim(cir)
  ∇c = TO.gen_jacobian(cir)
  RD.jacobian!(cir, ∇c, c, z)
  @show ∇c
  @show size(∇c)
  @show typeof(cir)

  #################################################

  l = 1.0
  off_cir = OffsetCircleConstraint(n, m, xc, yc, r, l)
  @show RD.evaluate(off_cir, z)
  @show RD.output_dim(off_cir)
  c = zeros(RD.output_dim(off_cir))
  RD.evaluate!(off_cir, c, z)
  @test c ≈ RD.evaluate(off_cir, z)
  ∇c = zeros(p, n)
  RD.jacobian!(off_cir, ∇c, c, z)
  @show ∇c
  @show size(∇c)
  @show typeof(off_cir)

  #################################################

  p = 5
  A = @SMatrix rand(p, 2)
  b = @SVector rand(p)
  off_lin = OffsetLinearConstraint(n, m, A, b, Inequality(), l)
  @show RD.evaluate(off_lin, z)
  @show RD.output_dim(off_lin)
  c = zeros(RD.output_dim(off_lin))
  RD.evaluate!(off_lin, c, z)
  @test c ≈ RD.evaluate(off_lin, z)
  ∇c = zeros(p, n)
  RD.jacobian!(off_lin, ∇c, c, z)
  @show ∇c
  @show typeof(off_lin)
  if off_lin isa OffsetLinearConstraint
  end

  #################################################
  p = 1
  xc = [9]
  yc = [1]
  a = [3]
  b = [1]
  θ = [deg2rad(30.0)]
  elli = EllipseConstraint(n, m, xc, yc, a, b, θ)
  @show RD.evaluate(elli, z)
  @show RD.output_dim(elli)
  c = zeros(RD.output_dim(elli))
  RD.evaluate!(elli, c, z)
  @test c ≈ RD.evaluate(elli, z)
  ∇c = zeros(p, n)
  RD.jacobian!(elli, ∇c, c, z)
  @show ∇c
  @show typeof(elli)

  #################################################
  p = 1
  xc = [9]
  yc = [1]
  a = [3]
  b = [1]
  θ = [deg2rad(30.0)]
  l = 2.0
  off_elli = OffsetEllipseConstraint(n, m, xc, yc, a, b, θ, l)
  @show RD.evaluate(off_elli, z)
  @show RD.output_dim(off_elli)
  c = zeros(RD.output_dim(off_elli))
  RD.evaluate!(off_elli, c, z)
  @test c ≈ RD.evaluate(off_elli, z)
  ∇c = zeros(p, n)
  RD.jacobian!(off_elli, ∇c, c, z)
  @show ∇c
  @show typeof(off_elli)

end

loop_for_display()

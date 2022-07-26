using Test

include("bicycle.jl")

# @testset "unit tests" begin

#   model = MyBicycleModel()
#   n, m = RD.dims(model)
#   x = zeros(n)
#   u = zeros(m)
#   t, h = 5, 0.1
#   z = KnotPoint(x, u, t, h)
#   println("-----------------------\n")

#   #################################################

#   p = 5
#   A = @SMatrix rand(p, n + m)
#   b = @SVector rand(p)
#   lin = LinearConstraint(n, m, A, b, Inequality())
#   @show RD.evaluate(lin, z)
#   @show RD.output_dim(lin)
#   @show typeof(lin)
#   println("-----------------------\n")

#   #################################################

#   xc = SA[1, 1, 1]
#   yc = SA[1, 2, 3]
#   r = SA[1, 1, 1]
#   cir = CircleConstraint(n, xc, yc, r)
#   c = zeros(3)
#   @show RD.evaluate(cir, z)
#   @show RD.output_dim(cir)
#   ∇c = TO.gen_jacobian(cir)
#   RD.jacobian!(cir, ∇c, c, z)
#   @show ∇c
#   @show size(∇c)
#   @show typeof(cir)
#   println("-----------------------\n")

#   #################################################

#   l = 2.0
#   off_cir = OffsetCircleConstraint(n, m, xc, yc, r, l)
#   @show RD.evaluate(off_cir, z)
#   @show RD.output_dim(off_cir)
#   c = zeros(RD.output_dim(off_cir))
#   RD.evaluate!(off_cir, c, z)
#   @test c ≈ RD.evaluate(off_cir, z)
#   ∇c = zeros(p, n)
#   RD.jacobian!(off_cir, ∇c, c, z)
#   @show ∇c
#   @show size(∇c)
#   @show typeof(off_cir)
#   println("-----------------------\n")

#   #################################################

#   points = [0 1.2 5 1.2 5;]
#   line_seg = LineSegmentConstraint(n, m, points, side = :lower)
#   @show RD.evaluate(line_seg, z)
#   @show RD.output_dim(line_seg)
#   c = zeros(RD.output_dim(line_seg))
#   RD.evaluate!(line_seg, c, z)
#   @test c ≈ RD.evaluate(line_seg, z)
#   ∇c = zeros(RD.output_dim(line_seg), n)
#   RD.jacobian!(line_seg, ∇c, c, z)
#   @show ∇c
#   @show typeof(line_seg)
#   println("-----------------------\n")

#   points = [0 -1.2 25 -1.2 25;]
#   line_seg = LineSegmentConstraint(n, m, points, side = :upper)
#   @show RD.evaluate(line_seg, z)
#   @show RD.output_dim(line_seg)
#   c = zeros(RD.output_dim(line_seg))
#   RD.evaluate!(line_seg, c, z)
#   @test c ≈ RD.evaluate(line_seg, z)
#   ∇c = zeros(RD.output_dim(line_seg), n)
#   RD.jacobian!(line_seg, ∇c, c, z)
#   @show ∇c
#   @show typeof(line_seg)
#   println("-----------------------\n")

#   #################################################

#   points = [0 1.2 5 1.2 5;]
#   off_line_seg = OffsetLineSegmentConstraint(n, m, points, l, side = :lower)
#   @show RD.evaluate(off_line_seg, z)
#   @show RD.output_dim(off_line_seg)
#   c = zeros(RD.output_dim(off_line_seg))
#   RD.evaluate!(off_line_seg, c, z)
#   @test c ≈ RD.evaluate(off_line_seg, z)
#   ∇c = zeros(RD.output_dim(off_line_seg), n)
#   RD.jacobian!(off_line_seg, ∇c, c, z)
#   @show ∇c
#   @show typeof(off_line_seg)
#   println("-----------------------\n")

#   points = [0 -1.2 25 -1.2 25;]
#   off_line_seg = OffsetLineSegmentConstraint(n, m, points, l, side = :upper)
#   @show RD.evaluate(off_line_seg, z)
#   @show RD.output_dim(off_line_seg)
#   c = zeros(RD.output_dim(off_line_seg))
#   RD.evaluate!(off_line_seg, c, z)
#   @test c ≈ RD.evaluate(off_line_seg, z)
#   ∇c = zeros(RD.output_dim(off_line_seg), n)
#   RD.jacobian!(off_line_seg, ∇c, c, z)
#   @show ∇c
#   @show typeof(off_line_seg)
#   println("-----------------------\n")

#   #################################################

#   p = 5
#   A = @SMatrix rand(p, 2)
#   b = @SVector rand(p)
#   off_lin = OffsetLinearConstraint(n, m, A, b, Inequality(), l)
#   @show RD.evaluate(off_lin, z)
#   @show RD.output_dim(off_lin)
#   c = zeros(RD.output_dim(off_lin))
#   RD.evaluate!(off_lin, c, z)
#   @test c ≈ RD.evaluate(off_lin, z)
#   ∇c = zeros(p, n)
#   RD.jacobian!(off_lin, ∇c, c, z)
#   @show ∇c
#   @show typeof(off_lin)
#   println("-----------------------\n")

#   #################################################
#   p = 1
#   xc = [9]
#   yc = [1]
#   a = [3]
#   b = [1]
#   θ = [deg2rad(30.0)]
#   elli = EllipseConstraint(n, m, xc, yc, a, b, θ)
#   @show RD.evaluate(elli, z)
#   @show RD.output_dim(elli)
#   c = zeros(RD.output_dim(elli))
#   RD.evaluate!(elli, c, z)
#   @test c ≈ RD.evaluate(elli, z)
#   ∇c = zeros(p, n)
#   RD.jacobian!(elli, ∇c, c, z)
#   @show ∇c
#   @show typeof(elli)
#   println("-----------------------\n")

#   #################################################
#   p = 1
#   xc = [9]
#   yc = [1]
#   a = [3]
#   b = [1]
#   θ = [deg2rad(30.0)]
#   l = 2.0
#   off_elli = OffsetEllipseConstraint(n, m, xc, yc, a, b, θ, l)
#   @show RD.evaluate(off_elli, z)
#   @show RD.output_dim(off_elli)
#   c = zeros(RD.output_dim(off_elli))
#   RD.evaluate!(off_elli, c, z)
#   @test c ≈ RD.evaluate(off_elli, z)
#   ∇c = zeros(p, n)
#   RD.jacobian!(off_elli, ∇c, c, z)
#   @show ∇c
#   @show typeof(off_elli)
#   println("-----------------------\n")

# end

loop_for_display()

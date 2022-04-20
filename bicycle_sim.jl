function loop_for_gif()
  x0 = SA_F64[0, 0, 0, 0, 4, 0]
  xf = SA[13, -0.5, deg2rad(0), 0, 1.0, 0]
  plt = plot([0, 20, 20, 0, 0], 2.0 * [-1.2, -1.2, 1.2, 1.2, -1.2])
  anim = @animate for i in UnitRange(1, 40)
    @show x0
    bicycle = BicycleCar(:cg, x0=x0, xf=xf, tf=5.0, i=i)
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

# (h, k): new center point, a: semimajor axes length, b: semiminor axes length,
# ψ: rotation angle
function ellipse!(h::Float64, k::Float64, a::Float64, b::Float64, ψ::Float64, plt::P) where {P}
  theta = (0.0:0.1:2*pi+0.1)
  s, c = sincos(ψ)
  xt = t -> cos(t)
  yt = t -> sin(t)
  fx = t -> c * a * xt(t) - s * b * yt(t) + h
  fy = t -> s * a * xt(t) + c * b * yt(t) + k
  points = [2.5 1.0
    -2.5 1.0
    -2.5 -1.0
    2.5 -1.0]
  ma = [c -s
    s c]
  points = points * ma
  points[:, 1] .+= h
  points[:, 2] .+= k
  points = vcat(points, points[1, :]')
  plot!(plt, points[:, 1], points[:, 2])
  # @testset "ellipse" begin
  #   f = (x, y) -> (x / a)^2 + (y / b)^2
  #   for i in theta
  #     x = fx(i)
  #     y = fy(i)
  #     # x′ = dot((c, s), (x, y)) + dot((c, s), (-h, -k))
  #     # y′ = dot((-s, c), (x, y)) + dot((-s, c), (-h, -k))
  #     x′ = c * x + s * y + (-c * h - s * k)
  #     y′ = -s * x + c * y + (s * h - c * k)
  #     @test f(x′, y′) ≈ 1.0
  #   end
  # end
  plot!(plt, [fx(i) for i in theta], [fy(i) for i in theta])
  # X, Y = [a * cos(i) for i in theta], [b * sin(i) for i in theta]
  # plot!(plt, c .* X .- s .* Y .+ x, s .* X .+ c .* Y .+ y)
end

function plot_ellipse_con!(c::C, plt::P) where {C<:Union{EllipseConstraint,OffsetEllipseConstraint},P}
  for i in 1:RD.output_dim(c)
    h = c.x[i]
    k = c.y[i]
    a = c.a[i]
    b = c.b[i]
    ψ = c.ψ[i]
    ellipse!(h, k, a, b, ψ, plt)
    half_width = 1.0
    if b > half_width
      ellipse!(h, k, a - half_width, b - half_width, ψ, plt)
    end
  end
end

function circle!(x::Float64, y::Float64, r::Float64, plt::P) where {P}
  theta = (0.0:0.1:2*pi+0.1)
  plot!(plt, [r * cos(i) + x for i in theta], [r * sin(i) + y for i in theta])
end

function plot_circle_con!(c::C, plt::P) where {C<:Union{CircleConstraint,OffsetCircleConstraint},P}
  for i in 1:RD.output_dim(c)
    x = c.x[i]
    y = c.y[i]
    r = c.radius[i]
    circle!(x, y, r, plt)
  end
end

function loop_for_display()
  gr()
  x0 = SA_F64[0, 0, 0, 0, 4, 0]
  xf = SA[15, -1.2, deg2rad(0), 0, 1.0, 0]
  plt = plot([0, 20], [-2.4, -2.4], aspect_ratio=:equal)
  plot!(plt, [0, 20], [-2.0, -2.0])
  plot!(plt, [0, 20], [2.4, 2.4])
  scatter!(plt, [xf[1]], [xf[2]], marker_size=2, shape=:star5)
  his_x = []
  his_y = []
  his_x_f = []
  his_y_f = []
  for i in UnitRange(1, 40)
    @show x0
    bicycle = BicycleCar(:cg, x0=x0, xf=xf, tf=5.0, i=i)
    solver = ALTROSolver(bicycle...)
    solve!(solver)
    X = states(solver)
    U = controls(solver)
    p = plot(plt, [x[1] for x in X], [x[2] for x in X])
    cons = get_constraints(bicycle[1])
    l = 2.0
    for con in cons
      @show typeof(con)
      if con isa LinearConstraint
      elseif con isa BoundConstraint
      elseif con isa CircleConstraint
        plot_circle_con!(con, p)
      elseif con isa OffsetLinearConstraint
        l = con.l
      elseif con isa OffsetCircleConstraint
        plot_circle_con!(con, p)
      elseif con isa EllipseConstraint
        plot_ellipse_con!(con, p)
      elseif con isa OffsetEllipseConstraint
        plot_ellipse_con!(con, p)
      end
    end
    @show size(X)
    @show size(U)
    push!(his_x, x0[1])
    push!(his_y, x0[2])
    push!(his_x_f, x0[1] + l * cos(x0[3]))
    push!(his_y_f, x0[2] + l * sin(x0[3]))

    r = 1.0
    circle!(his_x[end], his_y[end], r, p)
    circle!(his_x_f[end], his_y_f[end], r, p)
    scatter!(p, his_x, his_y)
    scatter!(p, his_x_f, his_y_f)
    # savefig(p, "pics/$(i)_haha.png")
    display(p)
    readline()
    x0 = X[2]
  end
end

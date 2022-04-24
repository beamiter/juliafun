function loop_for_gif()
  x0 = SA_F64[0, 0, 0, 0, 4, 0]
  xf = SA[13, -0.5, deg2rad(0), 0, 1.0, 0]
  plt = plot([0, 20, 20, 0, 0], 2.0 * [-1.2, -1.2, 1.2, 1.2, -1.2])
  anim = @animate for i in UnitRange(1, 40)
    @show x0
    bicycle = BicycleCar(:cg, x0 = x0, xf = xf, tf = 5.0, i = i)
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
  gif(anim, "anim_fps15.gif", fps = 1)
end

# (h, k): new center point, a: semimajor axes length, b: semiminor axes length,
# ψ: rotation angle
function ellipse!(
  h::Float64,
  k::Float64,
  a::Float64,
  b::Float64,
  ψ::Float64,
  plt::P;
  alpha = 1.0,
  color = :red,
  linewidth = 1.0,
) where {P}
  if a < 0.1 || b < 0.1
    return
  end
  theta = (0.0:0.1:2*pi+0.1)
  s, c = sincos(ψ)
  xt = t -> cos(t)
  yt = t -> sin(t)
  fx = t -> c * a * xt(t) - s * b * yt(t) + h
  fy = t -> s * a * xt(t) + c * b * yt(t) + k
  points = [
    2.5 1.0
    -2.5 1.0
    -2.5 -1.0
    2.5 -1.0
  ]
  ma = [
    c -s
    s c
  ]
  points = points * ma
  points[:, 1] .+= h
  points[:, 2] .+= k
  points = vcat(points, points[1, :]')
  plot!(plt, points[:, 1], points[:, 2], alpha = alpha, color = color, linewidth = linewidth)
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
  plot!(
    plt,
    [fx(i) for i in theta],
    [fy(i) for i in theta],
    alpha = alpha,
    color = color,
    linewidth = linewidth,
  )
  # X, Y = [a * cos(i) for i in theta], [b * sin(i) for i in theta]
  # plot!(plt, c .* X .- s .* Y .+ x, s .* X .+ c .* Y .+ y)
end

function plot_line_segment_con!(
  c::C,
  plt::P;
  pred = false,
) where {C<:Union{LineSegmentConstraint},P}
  for i = 1:RD.output_dim(c)
    x1 = c.points[i, 1]
    y1 = c.points[i, 2]
    x2 = c.points[i, 3]
    y2 = c.points[i, 4]
    if pred
      plot!(plt, [x1, x2], [y1, y2], alpha = 0.3, color = :red)
    else
      plot!(plt, [x1, x2], [y1, y2], linewidth = 1.5, color = :red)
    end
  end
end

function plot_ellipse_con!(
  c::C,
  plt::P;
  pred = false,
) where {C<:Union{EllipseConstraint,OffsetEllipseConstraint},P}
  for i = 1:RD.output_dim(c)
    h = c.x[i]
    k = c.y[i]
    a = c.a[i]
    b = c.b[i]
    ψ = c.ψ[i]
    if pred
      ellipse!(h, k, a, b, ψ, plt, alpha = 0.3, color = :azure3)
      ellipse!(h, k, a - 1.0, b - 1.0, ψ, plt, alpha = 0.3, color = :azure3)
    else
      ellipse!(h, k, a, b, ψ, plt, linewidth = 1.5)
      ellipse!(h, k, a - 1.0, b - 1.0, ψ, plt, linewidth = 1.5)
    end
  end
end

function circle!(x::Float64, y::Float64, r::Float64, plt::P) where {P}
  theta = (0.0:0.1:2*pi+0.1)
  plot!(plt, [r * cos(i) + x for i in theta], [r * sin(i) + y for i in theta])
end

function plot_circle_con!(c::C, plt::P) where {C<:Union{CircleConstraint,OffsetCircleConstraint},P}
  for i = 1:RD.output_dim(c)
    x = c.x[i]
    y = c.y[i]
    r = c.radius[i]
    circle!(x, y, r, plt)
  end
end

function loop_for_display()
  gr()
  x0 = SA_F64[0, 0, 0, 0, 4, 0]
  xf = SA[15, -1.0, deg2rad(0), 0, 1.0, 0]
  tf = 5.0
  N = 51
  plt = plot()
  plot!(plt, [0, 20], [-2.4, -2.4], aspect_ratio = :equal, alpha = 0.4)
  plot!(plt, [0, 20], [-2.0, -2.0], alpha = 0.4)
  plot!(plt, [0, 20], [2.4, 2.4], aspect_ratio = :equal, alpha = 0.4)
  line_segments_x = Vector{Float64}([0.0, 5.0, 10.0, 15.0])
  line_segments_y = Vector{Float64}([1.2, 1.2, 0.0, 0.0])
  @assert length(line_segments_x) == length(line_segments_y) "x, y length should be equal"
  line_segments = collect(zip(line_segments_x, line_segments_y))
  plot!(plt, line_segments_x, line_segments_y)
  scatter!(plt, [xf[1]], [xf[2]], marker_size = 2, shape = :star5)
  his_x = []
  his_y = []
  his_x_f = []
  his_y_f = []
  warm_up = true
  for i in UnitRange(1, 40)
    @show x0
    p = plot(plt)

    X0 = SA_F64[]
    U0 = SA_F64[0, 0]
    line_segment_points = []
    if warm_up
      # 1. use ilqr to warm up
      bicycle = BicycleCar(x0, xf, N, tf, i = i, constrained = false)
      solver = Altro.iLQRSolver(bicycle..., use_static = Val(true))
      solve!(solver)
      X = states(solver)
      plot!(
        p,
        [x[1] for x in X],
        [x[2] for x in X],
        linewidth = 1.5,
        linestyle = :dashdot,
        color = :blue,
      )
      line_segment_points = map(X) do x
        local dis = []
        local x0, y0 = x[1:2]
        for (x1, y1) in line_segments[1:end-1]
          push!(dis, hypot(x1 - x0, y1 - y0))
        end
        local id = argmin(dis)
        x1, y1 = line_segments[id]
        x2, y2 = line_segments[id+1]
        len = hypot(x2 - x1, y2 - y1)
        [x1, y1, x2, y2, len]
      end
      # 2. warm up
      X0 = states(solver)
      U0 = controls(solver)
    end
    bicycle = BicycleCar(
      x0,
      xf,
      N,
      tf,
      X0 = X0,
      U0 = U0,
      i = i,
      constrained = true,
      line_segment_points = line_segment_points,
    )
    # 3. use altro to solve
    solver = ALTROSolver(bicycle...)
    solve!(solver)
    X = states(solver)
    plot!(p, [x[1] for x in X], [x[2] for x in X], linewidth = 1.5, color = :green)

    l = 2.0
    con_list = get_constraints(bicycle[1])
    for (j, con) in enumerate(con_list.constraints)
      ids = con_list.inds[j]
      @assert length(ids) != 0
      # @show typeof(con)
      if con isa LinearConstraint
      elseif con isa BoundConstraint
      elseif con isa CircleConstraint
        plot_circle_con!(con, p)
      elseif con isa OffsetLinearConstraint
        l = con.l
      elseif con isa OffsetCircleConstraint
        plot_circle_con!(con, p)
      elseif con isa EllipseConstraint
        if ids[1] == 1
          plot_ellipse_con!(con, p)
        else
          plot_ellipse_con!(con, p, pred = true)
        end
      elseif con isa OffsetEllipseConstraint
        if ids[1] == 1
          plot_ellipse_con!(con, p)
        else
          plot_ellipse_con!(con, p, pred = true)
        end
      elseif con isa LineSegmentConstraint
        if ids[1] == 1
          plot_line_segment_con!(con, p)
        end
      end
    end
    push!(his_x, x0[1])
    push!(his_y, x0[2])
    push!(his_x_f, x0[1] + l * cos(x0[3]))
    push!(his_y_f, x0[2] + l * sin(x0[3]))

    r = 1.0
    circle!(his_x[end], his_y[end], r, p)
    circle!(his_x_f[end], his_y_f[end], r, p)
    scatter!(p, his_x, his_y, markersize = 2)
    scatter!(p, his_x_f, his_y_f, markersize = 2)
    # savefig(p, "pics/$(i)_haha.png")
    display(p)
    readline()
    x0 = X[2]
  end
end

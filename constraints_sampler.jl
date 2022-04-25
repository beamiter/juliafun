function revise_id(
  line_segments_x::Vector{T},
  line_segments_y::Vector{T},
  id::Int,
  x0::Float64,
  y0::Float64,
) where {T}
  x_prev, y_prev = line_segments_x[id-1], line_segments_y[id-1]
  x_cur, y_cur = line_segments_x[id], line_segments_y[id]
  x_next, y_next = line_segments_x[id+1], line_segments_y[id+1]
  x0_vec = (x0 - x_cur, y0 - y_cur)
  prev_vec = (x_prev - x_cur, y_prev - y_cur)
  next_vec = (x_next - x_cur, y_next - y_cur)
  cosine1 = dot(x0_vec, prev_vec) / (prev_vec[1]^2 + prev_vec[2]^2)
  cosine2 = dot(x0_vec, next_vec) / (next_vec[1]^2 + next_vec[2]^2)
  if cosine1 > cosine2
    id = id - 1
  end
  id
end

function get_line_segments_points(
  line_segments_x::Vector{T},
  line_segments_y::Vector{T},
  X,
  ;
  l = 2.0,
) where {T}
  @assert length(line_segments_y) == length(line_segments_x)
  seg_inds = map(X) do x
    dis = []
    dis_l = []
    x0, y0, θ = x[1:3]
    xl = x0 + l * cos(θ)
    yl = y0 + l * sin(θ)
    for i = 1:(length(line_segments_x)-1)
      x1 = line_segments_x[i]
      y1 = line_segments_y[i]
      push!(dis, (x1 - x0)^2 + (y1 - y0)^2)
      push!(dis_l, (x1 - xl)^2 + (y1 - yl)^2)
    end
    argmin(dis), argmin(dis_l)
  end
  line_segments_points = []
  line_segments_points_l = []
  for (i, x) in enumerate(X)
    x0, y0, θ = x[1:3]
    xl = x0 + l * cos(θ)
    yl = y0 + l * sin(θ)
    id, id_l = seg_inds[i]

    if id != 1 && id != length(X)
      id = revise_id(line_segments_x, line_segments_y, id, x0, y0)
    end
    x_1, y_1 = line_segments_x[id], line_segments_y[id]
    x_2, y_2 = line_segments_x[id+1], line_segments_y[id+1]
    len = hypot(x_1 - x_2, y_1 - y_2)
    push!(line_segments_points, [x_1, y_1, x_2, y_2, len])

    if id_l != 1 && id_l != length(X)
      id_l = revise_id(line_segments_x, line_segments_y, id_l, xl, yl)
    end
    x_1, y_1 = line_segments_x[id_l], line_segments_y[id_l]
    x_2, y_2 = line_segments_x[id_l+1], line_segments_y[id_l+1]
    len = hypot(x_1 - x_2, y_1 - y_2)
    push!(line_segments_points_l, [x_1, y_1, x_2, y_2, len])
  end
  line_segments_points, line_segments_points_l
end

function generate_constraints(n, m, N, tf, xf; i::Int = 1, X0 = SA_F64[])
  cons = ConstraintList(n, m, N)
  l = 2.0

  ######## bound constraint
  bnd_x_l = [-1, -2.4, Inf, -deg2rad(45), 0.0, -3]
  bnd_x_u = [20, 2.4, Inf, deg2rad(45), 6.0, 3]
  bnd_u_l = [-3, -deg2rad(45)]
  bnd_u_u = [3, deg2rad(45)]
  bnd = BoundConstraint(n, m, x_min = bnd_x_l, x_max = bnd_x_u, u_min = bnd_u_l, u_max = bnd_u_u)
  add_constraint!(cons, bnd, 1:N-1)

  ######## line segment constraint
  line_segments_x = Vector{Float64}([0.0, 5.0, 10.0, 15.0])
  line_segments_y = Vector{Float64}([1.2, 1.2, 0.0, 0.0])
  line_segment_points, line_segment_points_l =
    get_line_segments_points(line_segments_x, line_segments_y, X0, l = l)

  if length(line_segment_points) > 0
    @assert length(line_segment_points) == N
    # @show line_segment_points
    for (j, line_segment) in enumerate(line_segment_points)
      line_segment = SMatrix{1,5}(line_segment)
      line_seg = LineSegmentConstraint(n, m, line_segment, dist = 1.1)
      add_constraint!(cons, line_seg, j)
    end
  end
  points = [0 -2.4 25 -2.4 25;]
  line_seg = LineSegmentConstraint(n, m, points, side = :upper)
  add_constraint!(cons, line_seg, 1:N)

  if length(line_segment_points_l) > 0
    @assert length(line_segment_points_l) == N
    # @show line_segment_points_l
    for (j, line_segment) in enumerate(line_segment_points_l)
      line_segment = SMatrix{1,5}(line_segment)
      off_line_seg = OffsetLineSegmentConstraint(n, m, line_segment, l, dist = 1.1)
      add_constraint!(cons, off_line_seg, j)
    end
  end
  points = [0 -2.4 25 -2.4 25;]
  off_line_seg = OffsetLineSegmentConstraint(n, m, points, l, side = :upper, dist = 1.1)
  add_constraint!(cons, off_line_seg, 1:N)

  ######## linear constraint
  p = 2
  A = zeros(p, m + n)
  A[1, 5] = 1.0
  A[2, 5] = -1.0
  A = SMatrix{p,m + n,Float64}(A)
  b = zeros(p)
  b[1] = 4.0
  b[2] = -0.01
  b = SVector{p,Float64}(b)
  lin = LinearConstraint(n, m, A, b, Inequality())
  # add_constraint!(cons, lin, 2:N)

  p = 2
  A = zeros(p, 2)
  b = zeros(p)
  A[2] = -1.0
  b[1] = 2.4
  l = 2.0
  off_lin = OffsetLinearConstraint(n, m, A, b, Inequality(), l)
  # add_constraint!(cons, off_lin, 1:N)

  ######## ellipse constraint
  dt = tf / (N - 1)
  dist = 15.0 - (i - 1) * 2.0 * dt
  for j = 1:N
    xc = SA[dist]
    yc = SA[1.2]
    a = SA[4.0]
    b = SA[2.0]
    ψ = SA[deg2rad(0.0)]
    elli = EllipseConstraint(n, m, xc, yc, a, b, ψ)
    add_constraint!(cons, elli, j:j)
    dist -= 2.0 * dt
    off_elli = OffsetEllipseConstraint(n, m, xc, yc, a, b, ψ, l)
    add_constraint!(cons, off_elli, j:j)
  end

  ######## circle constraint
  xc = SA[7]
  yc = SA[0.5]
  r = SA[1.0]
  cir = CircleConstraint(n, xc, yc, r)
  # add_constraint!(cons, cir, 1:N)

  l = 2.0
  off_cir = OffsetCircleConstraint(n, m, xc, yc, r, l)
  # add_constraint!(cons, off_cir, 1:N)

  ######## goal constraint
  goal = GoalConstraint(xf)
  # add_constraint!(cons, goal, N)

  cons
end

function generate_constraints(n, m, N, tf, xf; i::Int = 1)
  cons = ConstraintList(n, m, N)

  ######## bound constraint
  bnd_x_l = [-1, -2.4, Inf, -deg2rad(45), 0.0, -3]
  bnd_x_u = [20, 2.4, Inf, deg2rad(45), 6.0, 3]
  bnd_u_l = [-3, -deg2rad(45)]
  bnd_u_u = [3, deg2rad(45)]
  bnd = BoundConstraint(n, m, x_min=bnd_x_l, x_max=bnd_x_u, u_min=bnd_u_l,
    u_max=bnd_u_u)

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

  p = 2
  A = zeros(p, 2)
  b = zeros(p)
  A[2] = -1.0
  b[1] = 2.4
  l = 2.0
  off_lin = OffsetLinearConstraint(n, m, A, b, Inequality(), l)

  ######## ellipse constraint
  p = 1
  dt = tf / (N - 1)
  dist = 15.0 - (i - 1) * 2.0 * dt
  xc = SA[dist, dist - 5.0]
  yc = SA[1.2, 1.2]
  a = SA[4.0, 3.0]
  b = SA[2.0, 1.5]
  ψ = SA[deg2rad(0.0), deg2rad(0.0)]
  elli = EllipseConstraint(n, m, xc, yc, a, b, ψ)

  l = 2.0
  off_elli = OffsetEllipseConstraint(n, m, xc, yc, a, b, ψ, l)

  ######## circle constraint
  xc = SA[7]
  yc = SA[0.5]
  r = SA[1.0]
  cir = CircleConstraint(n, xc, yc, r)

  l = 2.0
  off_cir = OffsetCircleConstraint(n, m, xc, yc, r, l)


  ######## goal constraint
  goal = GoalConstraint(xf)


  # add_constraint!(cons, goal, N)
  add_constraint!(cons, bnd, 1:N-1)
  add_constraint!(cons, lin, 2:N)
  # add_constraint!(cons, off_lin, 1:N)
  # add_constraint!(cons, cir, 1:N)
  # add_constraint!(cons, off_cir, 1:N)
  add_constraint!(cons, elli, 1:N)
  add_constraint!(cons, off_elli, 1:N)

  cons
end

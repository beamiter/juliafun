
function loop_for_gif()
    x0 = SA_F64[0, 0, 0, 0, 4, 0]
    plt = plot([0, 20, 20, 0, 0], 2.0 * [-1.2, -1.2, 1.2, 1.2, -1.2])
    anim = @animate for i in UnitRange(1, 40)
        @show x0
        bicycle = BicycleCar(:parallel_park, x0=x0)
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

function circle!(x::Float64, y::Float64, r::Float64, plt::P) where {P}
    theta = (0.0:0.1:2*pi+0.1)
    plot!(plt, [r * cos(i) + x for i in theta], [r * sin(i) + y for i in theta])
end
function plot_circle_con!(c::C, plt::P) where {C<:CircleConstraint,P}
    for i in 1:RD.output_dim(c)
        x = c.x[i]
        y = c.y[i]
        r = c.radius[i]
        println(x, ",", y, ",", r)
        circle!(x, y, r, plt)
    end
end

function loop_for_display()
    gr()
    x0 = SA_F64[0, 0, 0, 0, 4, 0]
    xf = SA[13, 0.6, deg2rad(0), 0, 0.1, 0]
    plt = plot([0, 20], [-2.4, -2.4])
    plot!(plt, [0, 20], [-2.0, -2.0])
    plot!(plt, [0, 20], [2.4, 2.4])
    scatter!(plt, [xf[1]], [xf[2]], marker_size=2)
    his_x = []
    his_y = []
    his_x_f = []
    his_y_f = []
    for i in UnitRange(1, 40)
        @show x0
        bicycle = BicycleCar(:parallel_park, x0=x0, xf=xf)
        solver = ALTROSolver(bicycle...)
        solve!(solver)
        X = states(solver)
        U = controls(solver)
        p = plot(plt, [x[1] for x in X], [x[2] for x in X])
        cons = get_constraints(bicycle[1])
        l = 1.0
        for con in cons
            @show typeof(con)
            if con isa LinearConstraint
            elseif con isa BoundConstraint
            elseif con isa CircleConstraint
                plot_circle_con!(con, p)
            elseif con isa OffsetLinearConstraint
                l = con.l
            else
            end
        end
        @show size(X)
        @show size(U)
        push!(his_x, x0[1])
        push!(his_y, x0[2])
        push!(his_x_f, x0[1] + l * cos(x0[3]))
        push!(his_y_f, x0[2] + l * sin(x0[3]))

        # r = 1.0
        # circle!(his_x[end], his_y[end], r, p)
        # circle!(his_x_f[end], his_y_f[end], r, p)
        scatter!(p, his_x, his_y)
        scatter!(p, his_x_f, his_y_f)
        display(p)
        sleep(0.1)
        readline()
        x0 = X[2]
    end
end


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
function loop_for_display()
    gr()
    x0 = SA_F64[0, 0, 0, 0, 4, 0]
    plt = plot([0, 20], [-2.4, -2.4])
    plot!(plt, [0, 20], [-2.0, -2.0])
    plot!(plt, [0, 20], [2.4, 2.4])
    his_x = []
    his_y = []
    his_x_f = []
    his_y_f = []
    for i in UnitRange(1, 40)
        @show x0
        bicycle = BicycleCar(:parallel_park, x0=x0)
        solver = ALTROSolver(bicycle...)
        solve!(solver)
        X = states(solver)
        U = controls(solver)
        cons = get_constraints(bicycle[1])
        for con in cons
            # @show typeof(con)
            if con isa LinearConstraint
                println("linear")
            elseif con isa BoundConstraint
                println("bound")
            elseif con isa CircleConstraint
                println("circle")
            elseif con isa OffsetLinearConstraint
                println("offset")
            else
                println(typeof(con))
            end
        end
        @show size(X)
        @show size(U)
        push!(his_x, x0[1])
        push!(his_y, x0[2])
        l = 1.0
        push!(his_x_f, x0[1] + l * cos(x0[3]))
        push!(his_y_f, x0[2] + l * sin(x0[3]))
        x0 = X[2]
        p = plot(plt, [x[1] for x in X], [x[2] for x in X])
        r = 1.0
        theta = (-pi:0.1:pi)
        plot!(p, [r * cos(i) + 7.0 for i in theta], [r * sin(i) + 0.5 for i in theta])
        r = 1.0
        plot!(p, [r * cos(i) + his_x[end] for i in theta], [r * sin(i) + his_y[end] for i in theta])
        plot!(p, [r * cos(i) + his_x_f[end] for i in theta], [r * sin(i) + his_y_f[end] for i in theta])
        scatter!(p, his_x, his_y)
        scatter!(p, his_x_f, his_y_f)
        # display(p)
        sleep(0.1)
        readline()
    end
end

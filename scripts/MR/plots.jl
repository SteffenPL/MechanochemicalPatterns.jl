function init_plot(s, p; show_frustration = false, show_heads = !show_frustration)

    fig = Figure(resolution = (1024, 768))

    state_obs = Observable(s)

    X_node = @lift $state_obs.X
    t_node = @lift @sprintf "Time = %.1fmin" $state_obs.t
    c_node = @lift [ atan(P[2] / P[1] ) for P in $state_obs.P ]
    H_node = @lift [$state_obs.X[seg[1]] for seg in $state_obs.colony]
    Hc_node = @lift [ atan($state_obs.P[seg[1]][2] / $state_obs.P[seg[1]][1] ) for seg in $state_obs.colony]
    frust_node = lift(state_obs) do s_
        F = zeros(length(s_.X))
        for (k, seg) in enumerate(s_.colony)
            for i in seg
                F[i] = s_.frustration[k]
            end
        end
        return F
    end 
    

    rm = p.cells.reversal_mechanism
    R = hasproperty(rm, :directional_density) ? rm.directional_density.R : p.cells.R_hard

    # create plot
    ax = Axis(fig[1,1], title = t_node)
    ax.aspect = DataAspect()

    xlims!(ax, p.env.domain.min[1], p.env.domain.max[1])
    ylims!(ax, p.env.domain.min[2], p.env.domain.max[2])


    if show_heads 
        scatter!(X_node, 
            markersize = 2 * p.cells.R_hard, markerspace = :data, 
            color = c_node, marker = Makie.Circle, colormap = :cyclic_mygbm_30_95_c78_n256)

            
        scatter!(X_node, 
            markersize = 2 * p.cells.R_hard, markerspace = :data, 
            color = (:black, 0.3), marker = Makie.Circle)


        scatter!(H_node, 
            markersize = 2 * R, markerspace = :data, 
            color = (:white, 0.5), marker = Makie.Circle, colormap = :cyclic_mygbm_30_95_c78_n256)

        scatter!(H_node, 
            markersize = 2 * p.cells.R_hard, markerspace = :data, 
            color = Hc_node, marker = Makie.Circle, colormap = :cyclic_mygbm_30_95_c78_n256)
    elseif show_frustration
        scatter!(H_node, 
            markersize = 2 * R, markerspace = :data, 
            color = (:white, 0.5), marker = Makie.Circle)

        sc = scatter!(X_node, 
            markersize = 2 * p.cells.R_hard, markerspace = :data, 
            color = frust_node, marker = Makie.Circle)

        Colorbar(fig[1,2], sc)
    else 
        scatter!(X_node, 
            markersize = 2 * p.cells.R_hard, markerspace = :data, 
            color = c_node, marker = Makie.Circle, colormap = :buda)

    end

    return fig, state_obs
end

function update_plot!(fig, state_obs, s, p)
    state_obs[] = s
end

function update_plot_callback!(fig, s_obs, dt = 0.1, timeout = 0.01)
    last_time = Ref(0.0)
    return function (s, p, cache)
        cur_time = time()
        if cur_time - last_time[] > dt
            last_time[] = cur_time
            update_plot!(fig, s_obs, s, p)
            sleep(timeout)
        end
    end
end

function add_slider!(fig, state_obs, states, p)
    i_slider = Slider(fig[end+1,1], range = 1:length(states), startvalue = 1)
    on(i_slider.value) do i
        update_plot!(fig, state_obs, states[i], p)
    end
    return i_slider
end

function play_animation!(fig, state_obs, states, skip = 10)
    for i in 1:skip:length(states)
        state_obs[] = states[i]
        sleep(0.01)
    end
end
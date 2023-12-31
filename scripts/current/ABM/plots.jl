function init_plot(s, p;
                    colors = [:magenta, :lightgreen],
                    bond_color = :darkorange,
                    show_concentration = false,
                    show_bonds = true,
                    show_polarities = true,
                    show_soft_spheres = false,
                    transparency = false,
                    bottom_plots = true,
                    alpha = clamp(200 / length(s.X), 0.01, 0.2)
                    )

    transparent_colors = transparency ? (c -> (c,alpha)).(colors) : colors 
    bond_color = transparency ? (bond_color, 0.5) : bond_color

    fig = Figure(resolution = (1024, 768))

    state_obs = Observable(s)

    SVecD = svec(p)

    X_node = @lift $state_obs.X
    ct = @lift $state_obs.cell_type
    E_node = @lift Tuple{SVecD,SVecD}[($state_obs.X[src(e)], $state_obs.X[dst(e)]) for e in edges($state_obs.adh_bonds) ]
    u_node = @lift $state_obs.u
    t_node = @lift @sprintf "Time = %.1fh" $state_obs.t
    R_node = @lift cache.R_hard[1:length($state_obs.X)]

    # create plot
    if dim(p) == 3
        ax = Axis3(fig[1,1], title = t_node, aspect = :data) # LScene(fig[1, 1])
    elseif dim(p) == 2 
        ax = Axis(fig[1,1], title = t_node) # LScene(fig[1, 1])
    end

    xlims!(ax, p.env.domain.min[1], p.env.domain.max[1])
    ylims!(ax, p.env.domain.min[2], p.env.domain.max[2])
    dim(p) == 3 && zlims!(ax, p.env.domain.min[3], p.env.domain.max[3])

    meshscatter!(X_node, 
                markersize = R_node, space = :data, 
                color = ct, colormap = colors; transparency)

    if show_soft_spheres
        meshscatter!(X_node, 
                    markersize = p.cells.R_soft, space = :data, 
                    color = ct, colormap = transparent_colors; transparency)
    end 

    if show_bonds
        linesegments!(E_node, color = bond_color, linewidth = 2; transparency)
    end 

    if show_polarities
        P_n = lift(state_obs) do s 
            Tuple{SVecD,SVecD}[(s.X[j], s.X[j] + p.cells.R_soft * s.P[j]) for j in eachindex(s.X) ]
        end
        
        linesegments!(P_n, color = :black; transparency)
    end

    if bottom_plots && dim(p) == 3
        dm = (p.env.domain.max[1], p.env.domain.max[2], p.env.domain.min[3])
        X_xy_node = @lift [ Point3f(x[1], x[2],dm[3]) for x in $state_obs.X]
        X_xz_node = @lift [ Point3f(x[1], dm[2],x[3]) for x in $state_obs.X]
        X_yz_node = @lift [ Point3f(dm[1],x[2], x[3]) for x in $state_obs.X]

        scatter!(X_xy_node, color = ct, colormap = colors, markersize = 5.0; transparency)
        scatter!(X_xz_node, color = ct, colormap = colors, markersize = 5.0; transparency)
        scatter!(X_yz_node, color = ct, colormap = colors, markersize = 5.0; transparency)
    end

    # create signal plot
    if show_concentration && dim(p) == 3
        ax2 = Axis3(fig[1,2], xlabel = "x", ylabel = "y", title = "u", aspect = :data)
        volume!(u_node, colorrange = (0,1.5))
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
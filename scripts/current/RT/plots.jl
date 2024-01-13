function init_plot(s, p, cache;
                    colors = [:magenta, :lightgreen],
                    bond_color = :darkorange,
                    show_concentration = false,
                    show_v = true,
                    n_peaks = 0,
                    show_bonds = true,
                    show_polarities = true,
                    show_soft_spheres = false,
                    transparency = false,
                    bottom_plots = true,
                    alpha = dim(p) == 2 ? 0.4 : clamp(200 / length(s.X), 0.01, 0.2)
                    )

    transparency = transparency || dim(p) == 2

    transparent_colors = transparency ? (c -> (c,alpha)).(colors) : colors 
    bond_color = transparency ? (bond_color, 0.5) : bond_color

    fig = Figure(resolution = (show_concentration ? 1500 : 1024 , 768))

    state_obs = Observable(s)

    SVecD = svec(p)

    function get_bond(state, edge)
        i, j = src(edge), dst(edge)
        bnd = (state.X[i], state.X[j])
        if dist(bnd[1], bnd[2]) > 10*p.cells.R_adh
            return (bnd[1], bnd[1])
        else
            return bnd
        end
    end

    X_node = @lift $state_obs.X
    ct = @lift $state_obs.cell_type
    E_node = @lift [get_bond($state_obs,e) for e in edges($state_obs.adh_bonds) ]
    u_node = @lift $state_obs.u
    v_node = @lift $state_obs.v
    t_node = @lift @sprintf "Time = %.1fh" $state_obs.t
    R_node = @lift 2*cache.R_hard[1:length($state_obs.X)]
    Rs_node = @lift 2*cache.R_soft[1:length($state_obs.X)]

    # create plot
    if dim(p) == 3
        ax = Axis3(fig[1,1], title = t_node, aspect = :data) # LScene(fig[1, 1])
    elseif dim(p) == 2 
        ax = Axis(fig[1,1], title = t_node) # LScene(fig[1, 1])
        ax.aspect = DataAspect()
    end

    xlims!(ax, p.env.domain.min[1], p.env.domain.max[1])
    ylims!(ax, p.env.domain.min[2], p.env.domain.max[2])
    dim(p) == 3 && zlims!(ax, p.env.domain.min[3], p.env.domain.max[3])

    if dim(p) == 3
        meshscatter!(X_node, 
                    markersize = R_node, space = :data, 
                    color = ct, colormap = colors; transparency)

        if show_soft_spheres
            meshscatter!(X_node, 
                        markersize = Rs_node, space = :data, 
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
    end

    if dim(p) == 2 

        if show_concentration

            xs, ys = LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid )
            heatmap!(xs, ys, u_node, colorrange = (0,1.5), colormap = :thermal, alpha = 0.5)
        end

        scatter!(X_node, 
                    markersize = R_node, markerspace = :data, marker = Makie.Circle, 
                    color = ct, colormap = colors)

        scatter!(X_node, 
                    markersize = Rs_node, markerspace = :data, marker = Makie.Circle,
                    color = ct, colormap = transparent_colors)


        if show_bonds
            linesegments!(E_node, color = bond_color, linewidth = 2; transparency)
        end 

        if show_polarities
            P_n = lift(state_obs) do s 
                Tuple{SVecD,SVecD}[(s.X[j], s.X[j] + p.cells.R_soft * s.P[j]) for j in eachindex(s.X) ]
            end
            
            linesegments!(P_n, color = :black; transparency)
        end

        if show_concentration

            xs, ys = LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid )
            ax2 = Axis(fig[1,2], xlabel = "x", ylabel = "y", title = "u")
            
            ax2.aspect = DataAspect()
            heatmap!(ax2, xs, ys, u_node, colorrange = (0,1.5), colormap = :thermal, interpolate = true)


            if n_peaks > 0
                peaks = @lift filteredpeaks($state_obs, p, n_peaks)
                pos = @lift @. Point2f(getproperty($peaks, :pos))
                scatter!(pos, color = :red, markersize = 10.0)
            end

            if show_v 
                ax3 = Axis(fig[1,3], xlabel = "x", ylabel = "y", title = "v")
                ax3.aspect = DataAspect()

                heatmap!(ax3, xs, ys, v_node, colorrange = (0,1.5), colormap = :thermal, interpolate = true)
            end

        end

        linkaxes!(ax, ax2)

        if show_v 
            linkaxes!(ax, ax3)
        end
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

begin
    # figure configuration
    colors = [:magenta, :lightgreen]
    bond_color = :black
    show_bonds = true
    show_soft_spheres = false
    transparency = false
    alpha = clamp(200 / length(states[end].X), 0.01, 0.2)

    transparent_colors = transparency ? (c -> (c,alpha)).(colors) : colors 
    bond_color = transparency ? (bond_color, 0.5) : bond_color

    fig = Figure(resolution = (1024, 768))

    i_slider = Slider(fig[2, 1], range = 1:length(states), startvalue = 1)
    
    X_node = @lift states[$(i_slider.value)].X
    ct = @lift states[$(i_slider.value)].cell_type

    E_node = lift(i_slider.value) do i 
        bonds = states[i].adh_bonds
        @show ne(bonds)
        Tuple{SVecD,SVecD}[ (states[i].X[src(e)], states[i].X[dst(e)]) for e in edges(bonds) ]
    end

    t_node = @lift @sprintf "Time = %.1fh" states[$(i_slider.value)].t

    ax = Axis3(fig[1,1], title = t_node, aspect = :data) # LScene(fig[1, 1])

    xlims!(ax, p.env.domain.min[1], p.env.domain.max[1])
    ylims!(ax, p.env.domain.min[2], p.env.domain.max[2])
    zlims!(ax, p.env.domain.min[3], p.env.domain.max[3])

    meshscatter!(X_node, 
                markersize = p.cells.R_hard, space = :data, 
                color = ct, colormap = colors; transparency)

    if show_soft_spheres
        meshscatter!(X_node, 
                    markersize = p.cells.R_soft, space = :data, 
                    color = ct, colormap = transparent_colors; transparency)
    end 

    if show_bonds
        linesegments!(E_node, color = bond_color, linewidth = 2; transparency)
    end 

    display(fig)

    for i in 1:10:length(states)
        i_slider.value[] = i
        sleep(0.01)
    end
end

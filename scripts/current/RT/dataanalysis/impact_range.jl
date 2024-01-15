function compute_kernel_density(s, p, ct, xs, ys)

    x = s.X[s.cell_type .== ct]
    X = getindex.(x, 1)
    Y = getindex.(x, 2)

    dom = p.env.domain
    boundary = ((dom.min[1],dom.max[1]),(dom.min[2],dom.max[2]))
    kd = kde( (X, Y), boundary = boundary, npoints=(length(xs), length(ys)))

    return kd
end




s = states[end]
peaks = filteredpeaks(s, p, 10)

p_ = @set p.plot.alpha_soft = 0.2

fig, s_obs = init_plot(s, p_, cache; show_polarities = true, bottom_plots = false, show_concentration = false)

xs = LinRange(p_.env.domain.min[1], p_.env.domain.max[1], p_.signals.grid[1])
ys = LinRange(p_.env.domain.min[2], p_.env.domain.max[2], p_.signals.grid[2])
data_smooth = [lift(s -> compute_kernel_density(s, p, i, xs, ys).density, s_obs) for i in 1:2]

ct_cols = [:magenta, :green]

for ct in 1:2
    Axis(fig[1,1+ct], xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(xs, ys, data_smooth[ct], colormap = cgrad([:black,ct_cols[ct]]), colorrange = (0,5e-6))
end 
display(fig)

add_slider!(fig, s_obs, states[2:end], p_)
using KernelDensity

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
p_ = @set p.plot.alpha_soft = 0.2
p_ = (p_..., signals = (;grid = (100, 100), min = (-250,-250)))

fig, s_obs = init_plot(s, p_, cache; show_polarities = true, bottom_plots = false, show_concentration = false)

xs = LinRange(p_.env.domain.min[1], p_.env.domain.max[1], p_.signals.grid[1])
ys = LinRange(p_.env.domain.min[2], p_.env.domain.max[2], p_.signals.grid[2])
data_smooth = [lift(s -> compute_kernel_density(s, p, i, xs, ys).density, s_obs) for i in 1:2]
peaks_pos = [lift(d -> [pk.pos for pk in filteredpeaks((;u=d), p_, 15)], data_smooth[i]) for i in 1:2]

ct_cols = [:magenta, :green]

for ct in 1:2
    Axis(fig[1,1+ct], xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(xs, ys, data_smooth[ct], colormap = cgrad([:black,ct_cols[ct]]), colorrange = (0,5e-6))
    scatter!(peaks_pos[ct], color = :red, markersize = 10.0)
end 
display(fig)

add_slider!(fig, s_obs, states[2:end], p_)


function compute_distances(states, p, ct, i1, i2, peak = missing)
    xs = LinRange(p.env.domain.min[1], p.env.domain.max[1], p.signals.grid[1])
    ys = LinRange(p.env.domain.min[2], p.env.domain.max[2], p.signals.grid[2])

    s1 = states[i1]
    s2 = states[i2]
    data = compute_kernel_density(s2, p, ct, xs, ys).density
    peaks = filteredpeaks((;u = data), p, 10)
    @show length(peaks)

    if !ismissing(peak)
        peaks = filteredpeaks((;u=data), p, 15)[peak:peak]
        @show peaks
    end

    pk_idx = [(i, argmin(dist²(pk.pos, s2.X[i]) for pk in peaks)) for i in eachindex(s2.X) if s2.cell_type[i] == ct && i <= length(s1.X)]
    
    dists2 = [dist(p, s2.X[i], peaks[j].pos) for (i, j) in pk_idx]
    dists1 = [dist(p, s1.X[i], peaks[j].pos) for (i, j) in pk_idx]

    return dists1, dists2
end

function digitalize(x, y, nbins)
    xmin = minimum(x) - eps(x[1])
    xmax = maximum(x)
    dx = (xmax - xmin) / nbins
    bins = LinRange(xmin + dx/2, xmax - dx/2, nbins)
    i_d = round.(Int, (x .- xmin) ./ dx) .+ 1
    y_hist = [mean(y[i_d .== i]) for i in 1:nbins]
    return bins, y_hist
end    

timeindex(states, t) = findfirst(s -> s.t >= t, states)

ts = (0.1, 30.0)

ct = 1

dist1, dist2 = compute_distances(states, p_, 2, timeindex(states, ts[1]), timeindex(states, ts[2]))

id = "chemo_2"
mkpath("scripts/current/RT/outputs/$(id)")
logscale = true
nbins = 40 
begin 
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "initial distance", ylabel = "distance delta towards peaks", yscale = logscale ? log10 : identity)
    logscale && ylims!(ax,(2.0^-4, 2^4))
    hlines!([logscale ? 1.0 : 0.0], color = :orange, linewidth = 5)
    data = logscale ? dist2 ./ dist1 : dist2 .- dist1
    scatter!(dist1, data, color = :gray, markersize = 5.0)

    if !logscale
        bins, data_means = digitalize(dist1, dist2 .- dist1, nbins)
        lines!(bins, data_means, color = :red, linewidth = 5.0)
    end
    save("""scripts/current/RT/outputs/$(id)/distance_$(logscale ? "log" : "lin")_delta_$(id)_$(ts[1])_$(ts[2]).png""", fig)
    display(fig)
end

vel(states, p, k, l, i) = wrap(p, states[k].X[i] - states[l].X[i]) / (states[k].t - states[l].t)

function compute_distances(states, p, ct, steps, i2, peaks; k_peaks = 10)
    s2 = states[i2]

    pk_idxs = [(i, argmin(dist²(pk.pos, s2.X[i]) for pk in peaks)) for i in eachindex(s2.X) if s2.cell_type[i] == ct && i <= length(states[steps[1]].X)]
    
    df = DataFrame(k = Int64[], i = Int64[], di = Float64[], align = Float64[], dir = Float64[])

    for k_step in eachindex(steps)[2:end]
        k = steps[k_step]
        s = states[k]
        for (i, pk_idx) in pk_idxs
            pk = peaks[pk_idx].pos

            Vi = vel(states, p, k, steps[k_step-1], i)
            Δpi = pk - s.X[i]
            di = dist(p, s.X[i], pk)

            align = dot(Vi, Δpi) / norm(Vi) / norm(Δpi)

            dir = atan(cross_2d(Vi / norm(Vi), Δpi / norm(Δpi)), align)
            push!(df, (k, i, di, align, dir))
        end 
    end

    return df
end

timeindex(states, t) = findfirst(s -> s.t >= t, states)
p = (p..., signals = (; grid = SVector{2,Int64}(100, 100)))

ct = 1
i2 = timeindex(states, 25.0)
Δt = timeindex(states, 0.1) - timeindex(states, 0.0) 

xs = LinRange(p.env.domain.min[1], p.env.domain.max[1], p.signals.grid[1])
ys = LinRange(p.env.domain.min[2], p.env.domain.max[2], p.signals.grid[2])

s2 = states[i2]
data = compute_kernel_density(s2, p, ct, xs, ys).density
peaks = filteredpeaks((;u = data), p, 10)[1:1]

begin 
    fig_data = Figure()
    Axis(fig_data[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())

    heatmap!(xs, ys, data, colormap = cgrad([:black, [:magenta,:green][ct]]), colorrange = (0,5e-6))
    scatter!(getproperty.(peaks, :pos), color = :red, markersize = 10.0)
    save("scripts/current/RT/outputs/align/$(id)/data_$(ct_name)_$(id).png", fig_data)
    fig_data
end

df = compute_distances(states, p, ct, 1:Δt:i2, i2, peaks; k_peaks = 1)
di_org = copy(df.di)
df.di = round.(df.di ./ 5, digits = 0) .* 5

df_mean = combine(groupby(df, :di), :dir => x -> cos(atan(mean( svec(p)(sincos(a)) for a in x)...)))

id = "pure_chemo_final4"
mkpath("scripts/current/RT/outputs/align/$(id)")
ct_name = ["Proximal", "Distal"][ct]


begin 
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Distance to peak (μm)", ylabel = "cos(θ)", title = "$(ct_name)")

    scatter!(di_org, df.align, color = :paleturquoise1, markersize = 2.0, alpha = 0.2)
    scatter!(df_mean.di, df_mean.dir_function, color = :red, markersize = 20.0)

    # violin!(round.(df.di ./ 40) .* 40, acos.(df.align), color = :red, width = 30)

    save("scripts/current/RT/outputs/align/$(id)/alignment_$(ct_name)_$(id).png", fig)
    fig
end
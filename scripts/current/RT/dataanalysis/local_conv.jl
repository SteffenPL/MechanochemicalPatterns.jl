include("../base.jl")
using KernelDensity, Dates

Random.seed!(1)
fn = "scripts/current/RT/inputs/parameters_2D.toml"
p = load_parameters(fn)
p = @set p.sim.t_end = 80.0

id = Dates.format(now(), "yyyy-mm-dd_HH-MM")

s = init_state(p)
cache = init_cache(p, s)

states = [deepcopy(s)]
simulate(s, p, cache; states = states)

fig, s_obs = init_plot(states[end], p, cache; show_polarities = true, show_concentration = true, n_peaks = 10)

#add_slider!(fig, s_obs, states, p)

# find main peak
s_end = states[end]
main_peak = filteredpeaks(, p, 20)
display(fig)

mkpath("scripts/current/RT/outputs/growth/$(id)")

for ct in 1:2
    ct_name = ["proximal", "distal"][ct]
    ct_col = [:magenta, :green][ct]

    x = s_end.X[s_end.cell_type .== ct]
    dom = p.env.domain
    xs = LinRange(dom.min[1], dom.max[1], p_.signals.grid[1])
    ys = LinRange(dom.min[2], dom.max[2], p_.signals.grid[2])
    kd = kde( (getindex.(x, 1), getindex.(x, 2)); boundary = ((dom.min[1],dom.max[1]),(dom.min[2],dom.max[2])), npoints=tuple(p_.signals.grid...))

    # compute number of cells within radius 
    R = 100.0
    fig = Figure()
    Axis(fig[1,1], xlabel = "x", ylabel = "y")
    current_axis().aspect = DataAspect()
    heatmap!(xs, ys, kd.density, colormap = cgrad([:black,ct_col]), colorrange = (0,5e-6))
    main_peak = filteredpeaks((;u=kd.density), p_, 10)[1]
    scatter!(main_peak.pos, color = :red, markersize = 10.0)
    arc!(main_peak.pos, R, 0.0, 2π, color = :red)
    save("scripts/current/RT/outputs/growth/$(id)/state_$(ct_name)_$(id)_R=$(R).png", fig)
    fig 


    # compute number of cells within radius
    function count_cells_within_radius(s, p, main_peak, R)
        return count(dist²(p, s.X[i], main_peak.pos) < R^2 for i in eachindex(s.X) if s.cell_type[i] == ct)
    end

    cell_count = [count_cells_within_radius(s, p, main_peak, R) for s in states]
    ts = [s.t for s in states]

    with_theme(Theme()) do
        fig2 = Figure()
        Axis(fig2[1,1], xlabel = "time", ylabel = "cell count", title = ct_name)
        lines!(ts, cell_count)
        save("scripts/current/RT/outputs/growth/$(id)/growth_$(ct_name)_$(id)_R=$(R).png", fig2)
        display(fig2)
    end
end


# repeat the same with multiple simulations 

n_reps = 10 


cell_counts = zeros(n_reps, length(states))

prog = Progress(n_reps, 1, "Repetitions...")
Threads.@threads for k_rep in 1:n_reps
    let cell_cpunts = cell_counts, p = p_m, R = R, k_rep = k_rep
        Random.seed!(k_rep)
        s = init_state(p)
        cache = init_cache(p, s)
        states = [deepcopy(s)]
        simulate(s, p, cache; states = states, show_prog = k_rep == 1)

        # find main peak
        s_end = states[end]
        main_peak = filteredpeaks(s_end, p, 1)[1]

        # count cells within radius
        for (k, s) in enumerate(states)
            cell_counts[k_rep, k] = count_cells_within_radius(s, p, main_peak, R)
        end
    end
    next!(prog)
end

rel_counts = cell_counts ./ [ cell_counts[:,end]; ]

cc_mean = mean(rel_counts, dims = 1)[:]
cc_std = std(rel_counts, dims = 1)[:]
ts = [s.t for s in states]

with_theme(Theme()) do
    fig2 = Figure()
    Axis(fig2[1,1], xlabel = "time (h)", ylabel = "cell count")
    lines!(ts, cc_mean, color = :black)
    band!(ts, cc_mean - cc_std, cc_mean + cc_std, color = :black, alpha = 0.2)
    save("scripts/current/RT/outputs/local_conv_avg_normalized_$(id).png", fig2)
    display(fig2)
end
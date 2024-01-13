include("../base.jl")

Random.seed!(1)
fn = "scripts/current/RT/inputs/parameters_2D.toml"
p = load_parameters(fn)
p = @set p.sim.t_end = 80.0

id = "blub"

s = init_state(p)
cache = init_cache(p, s)

states = [deepcopy(s)]
simulate(s, p, cache; states = states)

fig, s_obs = init_plot(states[end], p, cache; show_polarities = true, show_concentration = true, n_peaks = 10)

#add_slider!(fig, s_obs, states, p)

# find main peak
s_end = states[end]
main_peak = filteredpeaks(s_end, p, 10)[4]
display(fig)


# compute number of cells within radius 
R = 100.0

# plot circle
ax1 = content(fig[1,1])
ax2 = content(fig[1,3])
arc!(ax1, main_peak.pos, R, 0.0, 2π, color = :red)
arc!(ax2, main_peak.pos, R, 0.0, 2π, color = :red)
save("scripts/current/RT/outputs/local_conv_state_$(id)_R=$(R).png", fig)
display(fig)

# compute number of cells within radius
function count_cells_within_radius(s, p, main_peak, R)
    return count(dist²(p, s.X[i], main_peak.pos) < R^2 for i in eachindex(s.X) if s.cell_type[i] == 2)
end

cell_count = [count_cells_within_radius(s, p, main_peak, R) for s in states]
ts = [s.t for s in states]

with_theme(Theme()) do
    fig2 = Figure()
    Axis(fig2[1,1], xlabel = "time", ylabel = "cell count")
    lines!(ts, cell_count)
    save("scripts/current/RT/outputs/local_conv_single_$(id)_R=$(R).png", fig2)
    display(fig2)
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
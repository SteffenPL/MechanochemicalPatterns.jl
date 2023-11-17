
function time_step!(s, p, cache) 


    # update cache (in case of cell division)
    update_cache!(s, p, cache)
    updateboxes!(cache.st, s.X)

    # cell events
    add_bonds!(s, p, cache)
    remove_bonds!(s, p, cache)

    # reset forces
    reset_forces!(s, p, cache)


    # update forces
    compute_adhesive_forces!(s, p, cache)
    compute_interaction_forces!(s, p, cache)
    compute_gravity_forces!(s, p, cache)

    # add noise 
    add_random_forces!(s, p, cache)

    # add forces
    for i in eachindex(s.X)
        s.X[i] += cache.F[i] * p.sim.dt / p.env.damping
    end
    
    # detect instability 
    if any(x -> any(isnan(xk) for xk in x), s.X)
        error("Instability at time t = $(s.t)")
    end

    # detect out of bounds issues
    ma = p.sim.collision_detection.margin
    dom = p.env.domain 
    for i in eachindex(s.X)
        if !all(dom.min .- ma .* dom.size .<= s.X[i] .<= dom.max .+ ma .* dom.size)
            error("Out of bounds at time t = $(s.t) at position $(s.X[i])")
        end
    end

    # deal with constraints
    project_non_overlap!(s, p, cache)
    project_onto_domain!(s, p, cache)

end




function simulate(s, p, cache, callbacks = Function[])

    (; dt ) = p.sim

    s = deepcopy(s)
    states = [deepcopy(s)]

    resize_cache!(s, p, cache)

    t_unsaved = 0.0
    n_steps = Int(round(p.sim.t_end / dt))

    prog = Progress(n_steps, 1, "Simulating... ")
    for k_step in 1:n_steps 
        time_step!(s, p, cache)
        s.t += dt 
        if t_unsaved > p.sim.saveat
            push!(states, deepcopy(s))
            t_unsaved = 0.0
        end
        t_unsaved += dt
        next!(prog)
    end

    return states
end

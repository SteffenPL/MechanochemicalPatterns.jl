
function time_step!(s, p, cache) 

    # update cache (in case of cell division)
    update_cache!(s, p, cache)
    updatetable!(cache.st, s.X)


    # reset forces
    reset_forces!(s, p, cache)

    # update forces
    # compute_interaction_forces!(s, p, cache)

    # add noise 
    # add_random_forces!(s, p, cache)

    # polarity dynamics
    flip_bacteria!(s, p, cache)
    pull_polarities!(s, p, cache)

    # compute forces
    compute_bending_forces!(s, p, cache)
    add_self_prop!(s, p, cache)
    
    # add forces
    for i in eachindex(s.X)
        s.X[i] += cache.F[i] * p.sim.dt / p.env.damping
    end

    # deal with constraints
    project_non_overlap!(s, p, cache)
    project_bonds!(s, p, cache)
    project_onto_domain!(s, p, cache)
end

function simulate(s, p, cache; callbacks = Function[], states = [deepcopy(s)])

    (; dt ) = p.sim

    s = deepcopy(s)

    resize_cache!(s, p, cache)
    update_cache!(s, p, cache)
    project_onto_domain!(s, p, cache)

    t_unsaved = 0.0
    t_lazy = 0.0
    n_steps = Int(round( (p.sim.t_end-s.t) / dt))

    prog = Progress(n_steps, 1, "Simulating... ")
    for k_step in 1:n_steps 
        time_step!(s, p, cache)
        s.t += dt 

        # update solution vector
        if t_unsaved > p.sim.saveat
            push!(states, deepcopy(s))
            t_unsaved = 0.0

            # call callbacks as well
            foreach(f -> f(s, p, cache), callbacks)
        end
        t_unsaved += dt
        t_lazy += dt
        
        next!(prog)
    end

    return states
end

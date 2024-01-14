
function time_step!(s, p, cache) 

    cell_divisions!(s, p, cache)

    # update cache (in case of cell division)
    update_cache!(s, p, cache)
    updatetable!(cache.st, s.X)

    # track velocity 
    for i in eachindex(s.X)
        cache.dX[i] = s.X[i]
    end

    # cell events
    add_bonds!(s, p, cache)
    remove_bonds!(s, p, cache)

    # medium interaction 
    compute_neighbourhood!(s, p, cache)

    # do diffusion integration 
    add_source!(s, p, cache)
    ode_time_step!(s, p, cache)
    follow_source!(s, p, cache)

    # reset forces
    reset_forces!(s, p, cache)

    # update forces
    compute_adhesive_forces!(s, p, cache)
    compute_interaction_forces!(s, p, cache)
    # compute_gravity_forces!(s, p, cache)
    compute_medium_forces!(s, p, cache)

    # add noise 
    add_random_forces!(s, p, cache)
    add_self_prop!(s, p, cache)

    # add forces
    for i in eachindex(s.X)
        s.X[i] += cache.F[i] * p.sim.dt / p.env.damping
    end

    # deal with constraints
    project_non_overlap!(s, p, cache)
    project_onto_domain!(s, p, cache)

    # update velocity
    for i in eachindex(s.X)
        cache.dX[i] = (wrap(p, s.X[i] - cache.dX[i])) / p.sim.dt
    end

    # polarity dynamics
    update_polarity!(s, p, cache)

end

function simulate(s, p, cache; callbacks = Function[], states = [deepcopy(s)], show_prog = true)

    (; dt ) = p.sim

    s = deepcopy(s)

    resize_cache!(s, p, cache)

    t_unsaved = 0.0
    t_lazy = 0.0
    n_steps = Int(round( (p.sim.t_end-s.t) / dt))

    project_onto_domain!(s, p, cache)

    prog = Progress(n_steps, 1, "Simulating... ")
    for k_step in 1:n_steps 
        time_step!(s, p, cache)
        s.t += dt 

        # update solution vector
        if t_unsaved > p.sim.saveat
            push!(states, deepcopy(s))
            #push!(states, partialcopy(s, t_lazy > p.signals.saveat))
            t_unsaved = 0.0
            t_lazy = t_lazy > p.signals.saveat ? 0.0 : t_lazy

            # call callbacks as well
            foreach(f -> f(s, p, cache), callbacks)
        end
        t_unsaved += dt
        t_lazy += dt
        
        if show_prog
            next!(prog)
        end
    end

    return states
end

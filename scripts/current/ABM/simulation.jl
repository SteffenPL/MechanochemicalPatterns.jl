
function time_step!(s, p, cache) 


    # update cache (in case of cell division)
    update_cache!(s, p, cache)
    updatetable!(cache.st, s.X)

    # cell events
    add_bonds!(s, p, cache)
    remove_bonds!(s, p, cache)


    # do diffusion integration 
    add_source!(s, p, cache)
    ode_time_step!(s, p, cache)

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
    t_lazy = 0.0
    n_steps = Int(round(p.sim.t_end / dt))

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
        end
        t_unsaved += dt
        t_lazy += dt
        
        next!(prog)
    end

    return states
end
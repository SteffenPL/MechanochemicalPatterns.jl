using OrdinaryDiffEq.SciMLBase 

function ode_time_step!(s, p, cache)
    set_u!(cache.ode_integrator, s.u)
    step!(cache.ode_integrator, p.sim.dt)
    s.u .= cache.ode_integrator.u
    return nothing
end


function add_source!(s, p, cache)
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        if ct == 2 
            rel = (s.X[i] - p.env.domain.min) ./ p.env.domain.size 
            I = round.(Int64, rel .* p.signals.grid) 
            I = clamp.(I, (1,1,1), p.signals.grid)
            s.u[I[1],I[2],I[3]] += 0.8 * p.sim.dt
        end
    end
end

using OrdinaryDiffEq.SciMLBase 

function ode_time_step!(s, p, cache)
    set_u!(cache.ode_integrator, s.u)
    step!(cache.ode_integrator, p.sim.dt)
    s.u .= cache.ode_integrator.u
    return nothing
end

function indexat(s, p, cache, X, offset = 0)
    rel = (X - p.env.domain.min) ./ p.env.domain.size 
    I = round.(Int64, rel .* p.signals.grid) 
    I = clamp.(I, one.(p.signals.grid) .+ offset, p.signals.grid .- offset)
    return I
end

function add_source!(s, p, cache)
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        signal_emission = get_param(p, ct, :signal_emission, 0.0)
    
        I = indexat(s, p, cache, s.X[i])
        s.u[I...] += signal_emission * p.sim.dt
    end
end

cross_2d(a, b) = a[1]*b[2] - a[2]*b[1]

function follow_source!(s, p, cache)
    inv_dV =  p.signals.grid ./ p.env.domain.size
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        chemotaxis_strength = get_param(p, ct, :chemotaxis_strength, 0.0)

        I = indexat(s, p, cache, s.X[i], 1)  # only interior indices

        # get gradient
        grad = @SVector[ s.u[(I .+ (1,0))...] - s.u[(I .+ (-1,0))...], 
                         s.u[(I .+ (0,1))...] - s.u[(I .+ (0,-1))...]]  .* inv_dV
        
        grad_l = sqrt(sum( z->z^2, grad))

        if grad_l > 0.0
            grad_n = grad ./ grad_l

            angle_cos = dot(grad_n, s.P[i])
            angle_sign = cross_2d(grad_n, s.P[i]) > 0 ? -1 : 1

            P_tangent = @SVector[ -s.P[i][2], s.P[i][1] ]

            # update polarity
            s.P[i] += chemotaxis_strength * p.sim.dt * grad_l * (1 - angle_cos) * angle_sign * P_tangent
            s.P[i] = normalize(s.P[i])
        end
    end
end
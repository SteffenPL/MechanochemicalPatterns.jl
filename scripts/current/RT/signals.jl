using OrdinaryDiffEq.SciMLBase 

using ComponentArrays: ComponentArray

function ode_time_step!(s, p, cache)
    set_u!(cache.ode_integrator, ComponentArray( u = s.u, v = s.v) )
    step!(cache.ode_integrator, p.sim.dt)
    s.u .= cache.ode_integrator.u.u 
    s.v .= cache.ode_integrator.u.v
    return nothing
end

function indexat(s, p, cache, X, offset = 0)
    rel = (X - p.env.domain.min) ./ p.env.domain.size 
    I = round.(Int64, rel .* p.signals.grid) 
    I = clamp.(I, one.(p.signals.grid) .+ offset, p.signals.grid .- offset)
    return I
end

function indexpos(s, p, I)
    return p.env.domain.min .+ (I .- 0.5) ./ p.signals.grid .* p.env.domain.size
end

function add_source!(s, p, cache)
    inv_dvol =  1.0 / prod(p.env.domain.size ./ p.signals.grid) 
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        signal_emission = get_param(p, ct, :signal_emission, 0.0)
    
        I = indexat(s, p, cache, s.X[i])
        if ct == 1 
            s.u[I...] += signal_emission * inv_dvol * p.sim.dt
        else 
            s.u[I...] += signal_emission * inv_dvol * p.sim.dt
        end
    end
end

cross_2d(a, b) = a[1]*b[2] - a[2]*b[1]

@inline function fd_grad(u::AbstractArray{Float64,2}, I, inv_dV)
    return @SVector[ u[(I .+ (1,0))...] - u[(I .+ (-1,0))...], 
                     u[(I .+ (0,1))...] - u[(I .+ (0,-1))...]]  .* inv_dV
end

@inline function fd_grad(u::AbstractArray{Float64,3}, I, inv_dV)
    return @SVector[ u[(I .+ (1,0,0))...] - u[(I .+ (-1,0,0))...], 
                     u[(I .+ (0,1,0))...] - u[(I .+ (0,-1,0))...],
                     u[(I .+ (0,0,1))...] - u[(I .+ (0,0,-1))...]]  .* inv_dV
end

function follow_source!(s, p, cache)
    inv_dV =  p.signals.grid ./ p.env.domain.size
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        chemotaxis_strength = get_param(p, ct, :chemotaxis_strength, 0.0)

        I = indexat(s, p, cache, s.X[i], 1)  # only interior indices

        # get gradient
        grad_u = fd_grad(s.u, I, inv_dV)
        # grad_v = fd_grad(s.v, I, inv_dV)
                        
        grad = grad_u
        grad_l = sqrt(sum( z->z^2, grad))

        if grad_l > 0.0
            grad_n = grad ./ grad_l

            angle_cos = dot(grad_n, s.P[i])
            P_tangent = normalize(grad_n - s.P[i] * angle_cos)

            # update polarity
            s.P[i] += chemotaxis_strength * p.sim.dt * grad_l * (1 - angle_cos) * P_tangent
            s.P[i] = normalize(s.P[i])
        end
    end
end
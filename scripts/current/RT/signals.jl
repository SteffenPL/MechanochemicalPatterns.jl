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
    I = ceil.(Int64, rel .* p.signals.grid) 
    #I = clamp.(I, one.(p.signals.grid) .+ offset, p.signals.grid .- offset)
    return (mod(I[1], 1:p.signals.grid[1]), mod(I[2], 1:p.signals.grid[2]))
end


function indexat(s, p, cache, X::SVector{3,Float64}, offset = 0)
    rel = (X - p.env.domain.min) ./ p.env.domain.size 
    I = ceil.(Int64, rel .* p.signals.grid) 
    #I = clamp.(I, one.(p.signals.grid) .+ offset, p.signals.grid .- offset)
    return (mod(I[1], 1:p.signals.grid[1]), mod(I[2], 1:p.signals.grid[2]), mod(I[3], 1:p.signals.grid[3]))
end

function indexpos(s, p, I)
    return p.env.domain.min .+ (I .- 0.5) .* p.env.domain.size ./ p.signals.grid
end

function add_source!(s, p, cache)
    inv_dvol =  1.0 / prod(p.env.domain.size ./ p.signals.grid) 
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        signal_emission = get_param(p, ct, :signal_emission, 0.0)
    
        I = indexat(s, p, cache, s.X[i])
        if !hasproperty(p.signals.types, :v) || ct == 2 
            s.u[I...] += signal_emission * inv_dvol * p.sim.dt
        else 
            s.v[I...] += signal_emission * inv_dvol * p.sim.dt
        end
    end
end

cross_2d(a, b) = a[1]*b[2] - a[2]*b[1]

@inline function fd_grad(u::AbstractArray{Float64,2}, I, inv_dV, p)
    fix(i) = (mod(i[1], axes(u,1)), mod(i[2], axes(u,2)))
    
    return @SVector[ u[fix(I .+ (1,0))...] - u[fix(I .+ (-1,0))...], 
                     u[fix(I .+ (0,1))...] - u[fix(I .+ (0,-1))...]]  .* inv_dV  
    #TODO: fix /2 later
end

@inline function fd_grad(u::AbstractArray{Float64,3}, I, inv_dV, p)
    fix(i) = (mod(i[1], axes(u,1)), mod(i[2], axes(u,2)), mod(i[3], axes(u,3)))

    return @SVector[ u[fix(I .+ (1,0,0))...] - u[fix(I .+ (-1,0,0))...], 
                     u[fix(I .+ (0,1,0))...] - u[fix(I .+ (0,-1,0))...],
                     u[fix(I .+ (0,0,1))...] - u[fix(I .+ (0,0,-1))...]]  .* inv_dV
end

function follow_source!(s, p, cache)
    inv_dV =  p.signals.grid ./ p.env.domain.size
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        chemotaxis_strength = get_param(p, ct, :chemotaxis_strength, 0.0)

        I = indexat(s, p, cache, s.X[i], 0)
        if !p.env.periodic
            @error "Implement"
        end

        # get gradient
                      
        
        if !hasproperty(p.signals.types, :v) || ct == 2 
            grad = fd_grad(s.u, I, inv_dV, p)
        else
            grad = fd_grad(s.v, I, inv_dV, p)
        end

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
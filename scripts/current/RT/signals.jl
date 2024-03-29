using OrdinaryDiffEq.SciMLBase 

using ComponentArrays: ComponentArray

function ode_time_step!(s, p, cache)
    set_u!(cache.ode_integrator, s.U)
    step!(cache.ode_integrator, p.sim.dt)
    s.U .= cache.ode_integrator.u
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
    n_signals = length(p.signals.types) 

    cond_s = p.cells.cond_signal_emission

    for i in eachindex(s.X)
        ct = s.cell_type[i]
        signal_emission = get_param(p, ct, :signal_emission, 0.0)
    
        I = CartesianIndex(indexat(s, p, cache, s.X[i]))

        for k in 1:n_signals     
            s.U.x[k][I] += signal_emission[clamp(k, axes(signal_emission,1))] * inv_dvol * p.sim.dt
        end

        for (t,src,des,f) in cond_s
            if ct == t
                s.U.x[des][I] += f * s.U.x[src][I] * inv_dvol * p.sim.dt
            end
        end
    end
end

cross_2d(a, b) = a[1]*b[2] - a[2]*b[1]

@inline function fd_grad(u::AbstractArray{Float64,2}, I, inv_dV, p)
    fix(i) = p.env.periodic ? mod.(i, axes(u)) : clamp.(i, axes(u))
    
    return @SVector[ u[fix(I .+ (1,0))...] - u[fix(I .+ (-1,0))...], 
                     u[fix(I .+ (0,1))...] - u[fix(I .+ (0,-1))...]]  .* inv_dV  
    #TODO: fix /2 later
end

@inline function fd_grad(u::AbstractArray{Float64,3}, I, inv_dV, p)
    fix(i) = p.env.periodic ? mod.(i, axes(u)) : clamp.(i, axes(u))

    return @SVector[ u[fix(I .+ (1,0,0))...] - u[fix(I .+ (-1,0,0))...], 
                     u[fix(I .+ (0,1,0))...] - u[fix(I .+ (0,-1,0))...],
                     u[fix(I .+ (0,0,1))...] - u[fix(I .+ (0,0,-1))...]]  .* inv_dV
end

function update_gradients!(s, p, cache)
    if p.signals.main isa SVector 

        inv_dV =  p.signals.grid ./ p.env.domain.size
        for i in eachindex(s.X)
            I = indexat(s, p, cache, s.X[i], 0)            
            grad = zero(svec(p)) 
            for (k, factor) in enumerate(p.signals.main)     
                grad += fd_grad(s.U.x[k], I, inv_dV, p) * factor
            end
            cache.data.grad[i] = grad
        end
    else
        i_main = get(p.signals, :main, 1)
        inv_dV =  p.signals.grid ./ p.env.domain.size
        for i in eachindex(s.X)
            I = indexat(s, p, cache, s.X[i], 0)                   
            cache.data.grad[i] = fd_grad(s.U.x[i_main], I, inv_dV, p)
            # s.sox9[i] = s.U.x[i_main][I...]
        end
    end
end

function follow_source!(s, p, cache)
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        chemotaxis_strength = get_param(p, ct, :chemotaxis_strength, 0.0)
        grad = cache.data.grad[i]

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
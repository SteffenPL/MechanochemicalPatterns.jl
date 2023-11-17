
function reset_forces!(s, p, cache)
    for i in eachindex(cache.F)
        cache.F[i] = zero(SVecD)
    end
end

function add_random_forces!(s, p, cache)
    inv_sqrt_dt = 1/sqrt(p.sim.dt)
    for i in eachindex(s.X)
        cache.F[i] += randn(SVecD) * p.cells.sigma * inv_sqrt_dt
    end
end

function add_bonds!(s, p, cache, i, j, Xi, Xj, dij)
    if dij < 2*p.cells.R_adh && !has_edge(s.adh_bonds, i, j)
        if rand() < 1.0 - exp(-p.cells.new_adh_rate * p.sim.dt)
            add_edge!(s.adh_bonds, i, j, 0.0)
        end
    end
end

add_bonds!(s,p,cache) = apply_interaction_kernel!(s, p, cache, add_bonds!, 2*p.cells.R_adh)

function apply_interaction_kernel!(s, p, cache, fnc, R)
    for i in eachindex(s.X)
        Xi = s.X[i]
        for j in neighbours(cache.st, Xi, R_max)
            if i < j 
                Xj = s.X[j]
                dij = dist(Xi, Xj)
                fnc(s, p, cache, i, j, Xi, Xj, dij)
            end
        end
    end
end

function remove_bonds!(s, p, cache)
    bonds = s.adh_bonds
    for e in edges(bonds)
        i, j = src(e), dst(e)
        dij = dist(s.X[i], s.X[j])

        if bonds[i,j] > p.cells.adh_duration # || dij > 4*p.cells.R_adh
            rem_edge!(bonds, i, j)
        else
            bonds[i,j] += p.sim.dt
        end
    end
end

function compute_adhesive_forces!(s, p, cache)
    bonds = s.adh_bonds
    for e in edges(bonds)
        i, j = src(e), dst(e)
        Xi, Xj = s.X[i], s.X[j]
        k = cache.adhesion_stiffness[i] + cache.adhesion_stiffness[j]
        cache.F[i] += k * (Xj - Xi)
        cache.F[j] -= k * (Xj - Xi)
    end
end

function compute_gravity_forces!(s, p, cache)
    for i in eachindex(s.X)
        cache.F[i] += p.env.gravity
    end
end

function repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    Rij = cache.R_soft[i] + cache.R_soft[j] 
    if 0.0 < dij < Rij 
        k = cache.repulsion_stiffness[i] + cache.repulsion_stiffness[j]
        cache.F[i] += (Rij - dij) * k / dij * (Xi - Xj) 
        cache.F[j] -= (Rij - dij) * k / dij * (Xi - Xj) 
    end
end

function attraction_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    Rij = cache.R_attract[i] + cache.R_attract[j] 
    if 0.0 < dij < Rij
        k = cache.attraction_stiffness[i] + cache.attraction_stiffness[j]
        cache.F[i] -= (Rij - dij) * k / dij * (Xi - Xj) 
        cache.F[j] += (Rij - dij) * k / dij * (Xi - Xj) 
    end
end


function interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    attraction_kernel!(s, p, cache, i, j, Xi, Xj, dij)
end

function compute_interaction_forces!(s, p, cache)
    R_int = 2*p.cells.R_interact
    for i in eachindex(s.X)
        Xi = s.X[i]
        for j in neighbours(cache.st, Xi, R_max)  # 1:i-1
            if i < j 
                Xj = s.X[j]
                dij² = dist²(Xi, Xj)
                if dij² < R_int^2
                    interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, sqrt(dij²))
                end
            end
        end
    end
end


function project_non_overlap!(s, p, cache)
    for i in eachindex(s.X)
        Ri = cache.R_hard[i]
        for j in neighbours(cache.st, s.X[i], R_max) # 1:i-1
            if i < j 
                dij² = dist²(s.X[i], s.X[j])
                Rij = Ri + cache.R_hard[j]
                if 0.0 < dij² < Rij^2
                    dij = sqrt(dij²)
                    Xij = s.X[j] - s.X[i]
                    s.X[i] -= 0.5 * (Rij - dij) / dij * Xij
                    s.X[j] += 0.5 * (Rij - dij) / dij * Xij
                end
            end
        end
    end
end


function project_onto_domain!(s, p, cache)
    for i in eachindex(s.X)
        R_hard = cache.R_hard[i]
        s.X[i] = clamp.(s.X[i], p.env.domain.min .+ R_hard, p.env.domain.max .- R_hard)
    end     
end

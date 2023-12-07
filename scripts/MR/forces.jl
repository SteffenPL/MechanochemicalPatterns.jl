function reset_forces!(s, p, cache)
    for i in eachindex(cache.F)
        cache.F[i] = zero(SVecD)
    end
end

function add_self_prop!(s, p, cache)
    dt_factor = p.cells.v_self_prop * p.env.damping
    for i in eachindex(s.X)
        cache.F[i] += s.P[i] * dt_factor
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
        k = get_hetero_param(s, p, cache, i, j, :attraction_stiffness)
        cache.F[i] -= (Rij - dij) * k / dij * (Xi - Xj) 
        cache.F[j] += (Rij - dij) * k / dij * (Xi - Xj) 
    end
end


function interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    attraction_kernel!(s, p, cache, i, j, Xi, Xj, dij)
end

function compute_interaction_forces!(s, p, cache)
    R_int = p.cells.R_interact
    for i in eachindex(s.X)
        Xi = s.X[i]
        for (j, o) in neighbours_bc(cache.st, Xi, R_int)  # 1:i-1
            if i < j 
                Xj = s.X[j] + o
                dij² = dist²(Xi, Xj)
                if dij² < R_int^2
                    interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, sqrt(dij²))
                end
            end
        end
    end
end

function compute_bending_forces!(s, p, cache)
    n_disks = p.cells.n_disks
    R = 2*p.cells.R_hard
    k = p.cells.bending_stiffness

    for i_head in 1:n_disks:length(s.X)
        for i in i_head+1:i_head+n_disks-2

            e⁺ = wrap(p, s.X[i+1] - s.X[i])
            l⁺ = norm(e⁺)
            e⁺ /= l⁺

            e⁻ = wrap(p, s.X[i-1] - s.X[i])
            l⁻ = norm(e⁻)
            e⁻ /= l⁻

            if l⁺ > 0 && l⁻ > 0
                e⁺⁻ = dot(e⁺, e⁻)
                v⁻ = -k / R * (e⁺ - e⁺⁻ * e⁻)
                v⁺ = -k / R * (e⁻ - e⁺⁻ * e⁺)

                cache.F[i-1] += v⁻
                cache.F[i+1] += v⁺
                cache.F[i]   += -v⁺ - v⁻
            end
        end
    end
end

function pull_polarities!(s, p, cache)
    for e in edges(s.bonds)
        j, i = src(e), dst(e)  # follower (j) -> leader (i)
        Xij = wrap(p, s.X[i] - s.X[j])
        s.P[j] = safe_normalize(Xij, s.P[i])

        # set head polaritites
        if outdegree(s.bonds, i) == 0
            s.P[i] = safe_normalize(Xij, s.P[i])
        end
    end
end

function project_non_overlap!(s, p, cache)
    Rij = 2*p.cells.R_hard
    for i in eachindex(s.X)
        for (j, o) in neighbours_bc(p, cache, s.X[i], Rij) # 1:i-1
            if i < j 
                Xi = s.X[i]
                Xj = s.X[j] + o
                dij² = dist²(Xi, Xj)
                if 0.0 < dij² < Rij^2
                    dij = sqrt(dij²)
                    Xij = Xj - Xi
                    s.X[i] -= 0.5 * (Rij - dij) / dij * Xij
                    s.X[j] += 0.5 * (Rij - dij) / dij * Xij
                end
            end
        end
    end
end

function project_bonds!(s, p, cache)
    Rij = 2*p.cells.R_hard
    for e in edges(s.bonds)
        i, j = src(e), dst(e)
        Xij = wrap(p, s.X[j] - s.X[i])
        dij = dist(Xij)     
        if dij > 0   
            s.X[i] -= 0.5 * (Rij - dij) / dij * Xij
            s.X[j] += 0.5 * (Rij - dij) / dij * Xij
        end
    end
end

function project_onto_domain!(s, p, cache)
    if p.env.periodic
        for i in eachindex(s.X)
            s.X[i] = @. mod(s.X[i] - p.env.domain.min, p.env.domain.size) + p.env.domain.min
        end
    else
        R_hard = p.cells.R_hard
        for i in eachindex(s.X)
            s.X[i] = clamp.(s.X[i], p.env.domain.min .+ R_hard, p.env.domain.max .- R_hard)
        end     
    end
end

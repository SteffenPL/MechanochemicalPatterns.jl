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

function compute_bending_forces!(s, p, cache)
    n_disks = p.cells.n_disks
    R = 2*p.cells.R_hard
    k = p.cells.bending_stiffness

    for seg in s.colony
        for k in 2:length(seg)-1
            i = seg[k]
            i⁻ = seg[k-1]
            i⁺ = seg[k+1]

            e⁺ = wrap(p, s.X[i⁺] - s.X[i])
            l⁺ = norm(e⁺)
            e⁺ /= l⁺

            e⁻ = wrap(p, s.X[i⁻] - s.X[i])
            l⁻ = norm(e⁻)
            e⁻ /= l⁻

            if l⁺ > 0 && l⁻ > 0
                e⁺⁻ = dot(e⁺, e⁻)
                v⁻ = -k / R * (e⁺ - e⁺⁻ * e⁻)
                v⁺ = -k / R * (e⁻ - e⁺⁻ * e⁺)

                cache.F[i⁻] += v⁻
                cache.F[i⁺] += v⁺
                cache.F[i]   += -v⁺ - v⁻
            end
        end
    end
end

function pull_polarities!(s, p, cache)
    for seg in s.colony
        for (i, j) in segments(seg)
            Xji = wrap(p, s.X[i] - s.X[j])  # points from j to i (i.e. forward)

            # polarities of followers
            s.P[j] = safe_normalize(Xji, s.P[j])

            # set head polaritites
            if i == seg[1]
                s.P[i] = s.P[j]
            end
        end
    end
end

function project_non_overlap!(s, p, cache)
    Rij = 2*p.cells.R_hard
    for i in eachindex(s.X)
        for (j, o) in neighbours_bc(p, cache.st, s.X[i], Rij) # 1:i-1
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

    for bact in s.colony
        for (i,j) in segments(bact)
            Xij = wrap(p, s.X[j] - s.X[i])
            dij = dist(Xij)   
            if dij > 0   
                s.X[i] -= 0.5 * (Rij - dij) / dij * Xij
                s.X[j] += 0.5 * (Rij - dij) / dij * Xij
            end  
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

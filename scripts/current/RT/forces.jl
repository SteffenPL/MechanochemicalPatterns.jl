@inline function get_hetero_param(s, p, cache, i, j, sym)
    if s.cell_type[i] == s.cell_type[j]
        return getproperty(cache.data, sym)[i]
    else
        return getproperty(p.cells.interaction, sym)
    end
end

function reset_forces!(s, p, cache)
    for i in eachindex(s.F)
        s.F[i] = zero(svec(p))
    end
end

function add_center_gravity!(s, p, cache)
    for i in eachindex(s.X)
        s.F[i] += (s.X[i] - p.env.domain.center) * get(p.env, :gravity, 0.0)
    end
end

function add_random_forces!(s, p, cache)
    inv_sqrt_dt = 1/sqrt(p.sim.dt)
    for i in eachindex(s.X)
        s.F[i] += randn(svec(p)) * p.cells.sigma * inv_sqrt_dt
    end
end

function rotate(v, angle, forward)

end

# polarity dynamics 
function update_polarity!(s, p, cache)
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        s.P[i] += randn(svec(p)) * p.cells.sigma_p * sqrt(p.sim.dt)
        s.P[i] = safe_normalize(s.P[i])

        nF = norm(cache.data.dX[i])
        cil = get_param(p, ct, :CIL, 0.0)
        if cil > 0 && nF > 0 && cache.data.neighbour_count[i] < 5
            dFP = dot(s.P[i], cache.data.dX[i])/nF
            PT = nF * (1 - dFP) * safe_normalize(cache.data.dX[i] - dFP*nF * s.P[i])
            s.P[i] += cil * p.sim.dt * PT
            s.P[i] = safe_normalize(s.P[i])
        end

        # if p.cells.medium_active && p.cells.plithotaxis > 0 && cache.data.neighbour_count[i] > 6 
        #     nf = sqrt(sum(z -> z^2, s.F[i]))
        #     fp = @. s.F[i]/nf - s.P[i]
        #     align = dot(s.P[i], s.F[i]) / nf
        #     fp = fp
        #     fp /= norm(fp)

        #     s.P[i] += fp * (1-align) * p.cells.plithotaxis * p.sim.dt
        #     s.P[i] /= norm(s.P[i])
        # end
    end
end

function add_self_prop!(s, p, cache)
    dt_factor = p.cells.v_self_prop * p.env.damping
    for i in eachindex(s.X)
        if !p.cells.medium_active || cache.data.neighbour_count[i] > 6
            s.F[i] += s.P[i] * dt_factor
        else 
            s.F[i] += s.P[i] * dt_factor * p.cells.medium_slowdown
        end
    end
end

function add_bonds!(s, p, cache, i, j, Xi, Xj, dij)
    if dij < 2*p.cells.R_adh && !has_edge(s.adh_bonds, i, j)
        rate = s.cell_type[i] == s.cell_type[j] ? cache.data.new_adh_rate[i] : p.cells.interaction.new_adh_rate
        if rand() < 1.0 - exp(-rate * p.sim.dt)
            add_edge!(s.adh_bonds, i, j, 0.0)
        end
    end
end

add_bonds!(s,p,cache) = apply_interaction_kernel!(s, p, cache, add_bonds!, 2*p.cells.R_adh)

function apply_interaction_kernel!(s, p, cache, fnc, R)
    R_int = R
    for i in eachindex(s.X)
        Xi = s.X[i]
    
        for (j, o) in neighbours_bc(p, cache.st, Xi, R)
            if i < j 
                Xj = s.X[j] - o
                dij = dist(p, Xi, Xj)
                fnc(s, p, cache, i, j, Xi, Xj, dij)
            end
        end
    end
end

function remove_bonds!(s, p, cache)
    bonds = s.adh_bonds
    for e in edges(bonds)
        i, j = src(e), dst(e)
        dij = dist(p, s.X[i], s.X[j])
        rate = get_hetero_param(s, p, cache, i, j, :break_adh_rate)
        
        if rand() < 1.0 - exp(-rate * p.sim.dt) || dij > 2*p.cells.R_adh
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
        
        Xij = wrap(p, s.X[j] - s.X[i])
        l = dist(p, s.X[i], s.X[j])
        # k = get_hetero_param(s, p, cache, i, j, :adhesion_stiffness)

        # 1. we compute the attraction willingness of both cells,
        # that depends on the gradient signal and bias 

        if l > 0
            # chem = s.U.x[p.signals.main]
            # cti = s.cell_type[i]
            # ctj = s.cell_type[j]

            # Ai = chem[indexat(s, p, cache, s.X[i], 0)]  
            # Aj = chem[indexat(s, p, cache, s.X[j], 0)]

            # Gi = cache.data.grad[i]
            # Gj = cache.data.grad[j]

            # Bi = norm(Gi)
            # Bj = norm(Gj)

            # ai = cache.data.attraction_chemo_base[i]
            # aj = cache.data.attraction_chemo_base[j]

            # bi = cache.data.attraction_chemo_bias[i]
            # bj = cache.data.attraction_chemo_bias[j]

            # f1 = Bi > 0 ? dot(Gi,  Xij) / (l*Bi) : 0.0
            # f2 = Bj > 0 ? dot(Gj, -Xij) / (l*Bj) : 0.0
            # # note -1 ≤ f1, f2 ≤ 1

            # k = get_hetero_param(s, p, cache, i, j, :adhesion_stiffness)
            # k_extr = ai * Ai + aj * Aj + bi * f1 * (0.1 + Bi) + bj * f2 * (0.1 + Bj)
            # if k_extr > 0.1
            #     @show i, j, k_extr
            # end
            # print(k_extr, " ")

            # k = k + clamp(k_extr, 0, k)

            # s.F[i] += k * Xij
            # s.F[j] -= k * Xij
        end
    end
end

function compute_gravity_forces!(s, p, cache)
    for i in eachindex(s.X)
        s.F[i] += p.env.gravity
    end
end

function repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    Rij = cache.data.R_soft[i] + cache.data.R_soft[j] 
    if 0.0 < dij < Rij 
        k = cache.data.repulsion_stiffness[i] + cache.data.repulsion_stiffness[j]
        Xji = wrap(p, s.X[i] - s.X[j])
        s.F[i] += (Rij - dij) * k / dij * Xji
        s.F[j] -= (Rij - dij) * k / dij * Xji
    end
end

function attraction_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    Rij = cache.data.R_attract[i] + cache.data.R_attract[j] 
    if 0.0 < dij < Rij
        Xij = wrap(p, s.X[j] - s.X[i])
        l = sqrt(sum(x->x^2, Xij))

        cti = s.cell_type[i]
        ctj = s.cell_type[j]

        Ai = s.sox9[i]
        Aj = s.sox9[j]

        Gi = cache.data.grad[i]
        Gj = cache.data.grad[j]

        Bi = norm(Gi)
        Bj = norm(Gj)

        ai = cache.data.attraction_chemo_base[i]
        aj = cache.data.attraction_chemo_base[j]

        bi = cache.data.attraction_chemo_bias[i]
        bj = cache.data.attraction_chemo_bias[j]

        fi = dot(s.P[i],  Xij) / l
        fj = dot(s.P[j], -Xij) / l
        # note -1 ≤ f1, f2 ≤ 1

        # adhesion potential per cell 
        C1 = clamp(1 + ai + fi*bi*Bi, 0.0, 1.5)
        C2 = clamp(1 + aj + fj*bj*Bj, 0.0, 1.5)        
        f = C1 * C2

        k = get_hetero_param(s, p, cache, i, j, :attraction_stiffness)
        k = ( f < 1 ) ? (f * k) : (k + (f-1) * (p.cells.chemo_max-k))

        s.F[i] += k / dij * Xij
        s.F[j] -= k / dij * Xij
    end
end


function interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    repulsion_kernel!(s, p, cache, i, j, Xi, Xj, dij)
    attraction_kernel!(s, p, cache, i, j, Xi, Xj, dij)
end

function compute_interaction_forces!(s, p, cache)
    R_int = 2 * max(p.cells.R_soft, p.cells.R_adh)

    for i in eachindex(s.X)
        Xi = s.X[i]
        for (j, o) in neighbours_bc(p, cache.st, Xi, R_int)  # 1:i-1 
            if i < j
                Xj = s.X[j] - o
                dij² = dist²(p, Xi, Xj)
                if 0.0 < dij² < R_int^2
                    interaction_force_kernel!(s, p, cache, i, j, Xi, Xj, sqrt(dij²))
                end
            end
        end
    end
end

function compute_neighbourhood!(s, p, cache)
    R_int = 2*p.cells.R_interact
    cache.data.neighbour_count .= 0

    for i in eachindex(cache.data.neighbour_avg)
        cache.data.neighbour_avg[i] = zero(svec(p))
    end

    if !p.cells.medium_active 
        return 
    end

    for i in eachindex(s.X)
        Ri = cache.data.R_hard[i]
        for (j, o) in neighbours_bc(p, cache.st, s.X[i], R_int) # 1:i-1
            if i < j 
                Xij = s.X[j] - o - s.X[i]
                dij = sqrt(sum(z -> z^2, Xij))
                Rij = Ri + cache.data.R_hard[j]
                if 0 < dij < 2*R_int
                    factor = Rij^p.cells.medium_alpha / dij^(1+p.cells.medium_alpha)
                    cache.data.neighbour_avg[i] += Xij * factor
                    cache.data.neighbour_avg[j] -= Xij * factor

                    cache.data.neighbour_count[i] += 1
                    cache.data.neighbour_count[j] += 1
                end
            end
        end
    end

    # normalize
    for i in eachindex(cache.data.neighbour_avg)
        nc = cache.data.neighbour_count[i]
        cache.data.neighbour_avg[i] *= (nc > 0 ? 1/nc : 0.0)
    end

    return nothing
end

function compute_medium_forces!(s, p, cache)
    for i in eachindex(s.X)
        nl = sqrt(sum(x->x^2, cache.data.neighbour_avg[i]))
        s.F[i] += p.cells.medium_repulsion * cache.data.neighbour_avg[i] * nl
    end
end

function project_non_overlap!(s, p, cache)
    Rij = 2*p.cells.R_hard
    for i in eachindex(s.X)
        Xi = s.X[i]
    
        for (j, o) in neighbours_bc(p, cache.st, s.X[i], Rij) # 1:i-1
            if i < j 
                Xi = s.X[i]
                Xj = s.X[j] - o
                dij² = dist²(p, Xi, Xj)
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

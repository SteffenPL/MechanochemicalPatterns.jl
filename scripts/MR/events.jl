function add_disk!(s, p, cache, i_bact, Xj, Pj)
    # add disk to the end of the colony
    push!(s.X, Xj)
    push!(s.P, Pj)

    # add bond to the end of the colony
    push!(s.colony[i_bact], length(s.X))

    cache.outdated = true
end



switch(p, clock, age) = false
switch(p, clock::@NamedTuple{constant::Float64}, age) = age > clock.constant
switch(p, clock::@NamedTuple{random_rate::Float64}, age) = rand() < 1 - exp(-clock.random_rate*p.sim.dt)

function flip_bacteria!(s, p, cache, k)
    s.colony[k] = reverse(s.colony[k])
    for j in s.colony[k]
        s.P[j] = -s.P[j]
    end
end 

function compute_head_neighbours!(s, p, cache, k)
    R = 2*p.cells.reversal_mechanism.directional_density.R
    count = 0
    polarity = 0.0
    for (j, o) in neighbours_bc(p, cache.st_heads, cache.Heads[k], R)
        djk² = dist²(p, cache.Heads[j], cache.Heads[k])
        if 0 < djk² < R^2

            Pj = s.P[s.colony[j][1]]
            Pk = s.P[s.colony[k][1]]

            count += 1
            polarity += dot(Pk, Pj) / (norm(Pj)*norm(Pk))

        end
    end
    polarity = count > 0 ? polarity/count : 1.0
    return polarity 
end

function flip_bacteria!(s, p, cache)
    for (i, seg) in enumerate(s.colony)
        cache.Heads[i] = s.X[seg[1]]
    end
    updatetable!(cache.st_heads, cache.Heads)

    cache.flipping .= false

    rm = p.cells.reversal_mechanism
    for k in eachindex(s.colony)
        if hasproperty(rm, :directional_density)
            dd = rm.directional_density
            polarity = compute_head_neighbours!(s, p, cache, k)
            
            if polarity < dd.threshold
                cache.flipping[k] = true
            end
        elseif hasproperty(rm, :clock)
            if s.tsr[k] > p.cells.T_clock.constant
                cache.flipping[k] = true
                s.tsr[k] = 0.0
            else 
                s.tsr[k] += p.sim.dt
            end
        elseif hasproperty(rm, :random)
            dd=rm.random
            if rand() < 1 - exp(-dd.random_rate*p.sim.dt)
                cache.flipping[k] = true
            end
        end
    end

    for (k, flip) in enumerate(cache.flipping)
        if flip
            flip_bacteria!(s, p, cache, k)
        end
    end
end
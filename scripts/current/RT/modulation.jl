function modulate_sym!(s, p, cache, sym, p_mod)

    if hasproperty(p_mod, :start) && s.t < p_mod.start
        return nothing  
    elseif hasproperty(p_mod, :end) && s.t > p_mod.end
        return nothing
    end

    if hasproperty(p_mod, :speed)
        X = getproperty(cache, sym)
        speed = p_mod.speed
        range = p_mod.range

        for i in eachindex(s.X)
            ct = s.cell_type[i]
            X[i] = clamp(X[i] + speed[ct] * p.sim.dt, range[1], range[2])
        end
    elseif hasproperty(p_mod, :random_rate)
        X = getproperty(s, sym)
        rate = p_mod.random_rate

        for i in eachindex(s.X)
            ct = s.cell_type[i]
            if rand() < 1.0 - exp(-rate[ct] * p.sim.dt)
                X[i] = p_mod.state[ct]
            end
        end
    end
end

function modulate_parameters!(s, p, cache)
    if hasproperty(p, :mods) 
        map(sym -> modulate_sym!(s,p,cache,sym,getproperty(p.mods,sym)), keys(p.mods))
    end
end
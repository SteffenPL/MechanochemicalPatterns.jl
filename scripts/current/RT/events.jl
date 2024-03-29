function cell_divisions!(s, p, cache)
    s.cell_age .+= p.sim.dt

    for i in eachindex(s.X)
        if s.cell_age[i] > p.cells.lifespan
            if dim(p) == 2 
                X = s.X[i]
                P = s.P[i]

                rd = random_direction(dim(p))
                s.X[i] += rd * p.cells.R_hard / 2
                s.P[i] = random_direction(dim(p))
                s.cell_age[i] = 0.0

                # create new cell
                ct = s.cell_type[i]

                push!(s.X, X - rd * p.cells.R_hard / 2)
                push!(s.P, random_direction(dim(p)))
                push!(s.cell_type, ct)
                push!(s.cell_age, 0.0)
                push!(s.sox9, s.sox9[i])
                push!(s.F, zero(svec(p)))

                add_vertex!(s.adh_bonds)

                cache.outdated = true
            elseif dim(p) == 3 
                X = s.X[i]
                P = s.P[i]

                rd = random_direction(dim(p))
                s.X[i] += rd * p.cells.R_hard / 2
                s.P[i] = random_direction(dim(p))
                s.cell_age[i] = 0.0

                # create new cell
                ct = s.cell_type[i]

                push!(s.X, X - rd * p.cells.R_hard / 2)
                push!(s.P, random_direction(dim(p)))
                push!(s.cell_type, ct)
                push!(s.cell_age, 0.0)
                push!(s.sox9, s.sox9[i])
                push!(s.F, zero(svec(p)))

                add_vertex!(s.adh_bonds)

                cache.outdated = true
            end
        end
    end

    if cache.outdated

        resize_cache!(s, p, cache)
    end
end


function update_cell_internals!(s, p, cache)
    sox_limit = 1.0
    for i in eachindex(s.X)
        ct = s.cell_type[i]
        k_sox9 = cache.data.k_sox9[i]
        if s.sox9[i] < sox_limit
            s.sox9[i] += k_sox9 * p.sim.dt
        end
    end
end
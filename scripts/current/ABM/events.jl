function cell_divisions!(s, p, cache)
    s.cell_age .+= p.sim.dt

    for i in eachindex(s.X)
        if s.cell_age[i] > p.cells.lifespan
            if Dim == 2 
                error("Cell division not implemented for 2D")
            elseif Dim == 3 
                X = s.X[i]
                P = s.P[i]

                rd = random_direction(Dim)
                s.X[i] += rd * p.cells.R_hard / 2
                s.P[i] = random_direction(Dim)
                s.cell_age[i] = 0.0

                # create new cell
                ct = s.cell_type[i]

                push!(s.X, X - rd * p.cells.R_hard / 2)
                push!(s.P, random_direction(Dim))
                push!(s.cell_type, ct)
                push!(s.cell_age, 0.0)

                add_vertex!(s.adh_bonds)

                cache.outdated = true
            end
        end
    end

    if cache.outdated

        resize_cache!(s, p, cache)
    end
end
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

function flip_bacteria!(s, p, cache)
    n_disks = p.cells.n_disks
    for i in 1:n_disks:length(s.X)
        if rand() < 1 - exp(-p.cells.flip_rate * p.sim.dt)
            # detect direction 
            if indegree(s.bonds, i) == 1 
                # i is the head 
                for j in i+1:i+n_disks-1
                    rem_edge!(s.bonds, j, j-1)
                    add_edge!(s.bonds, j-1, j)
                end
            else
                # i is the tail, i + n_disks - 1 is the head
                for j in i+1:i+n_disks-1
                    rem_edge!(s.bonds, j-1, j)
                    add_edge!(s.bonds, j, j-1)
                end
            end

            for j in i:i+n_disks-1
                s.P[j] = -s.P[j]
            end
        end
    end
end
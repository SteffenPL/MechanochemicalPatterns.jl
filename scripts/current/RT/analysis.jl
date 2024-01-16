using Statistics 

function avg_location(s, ct = 1:2)
    return mean(s.X[i] for i in eachindex(s.X) if s.cell_type[i] in ct)
end

function avg_locations(states, ct = 1:2)
    return [avg_location(s, ct) for s in states]
end

# begin 
#     Figure()
#     yspan = (-80, 20)
#     Axis(current_figure()[1,1])
#     scatter!(map(p -> p[1], avg_locations(states, 1:2)), color = :white)
#     scatter!(map(p -> p[1], avg_locations(states, 1)), color = :magenta)
#     scatter!(map(p -> p[1], avg_locations(states, 2)), color = :lightgreen)
#     ylims!(yspan...)

#     Axis(current_figure()[1,2])
#     scatter!(map(p -> p[2], avg_locations(states, 1:2)), color = :white)
#     scatter!(map(p -> p[2], avg_locations(states, 1)), color = :magenta)
#     scatter!(map(p -> p[2], avg_locations(states, 2)), color = :lightgreen)
#     ylims!(yspan...)

#     display(current_figure())
# end




# peaks 
import Images 

function atboundary(u, I)
    return any( (i == 1 || i == size(u,k)) for (k, i) in enumerate(I))
end


function filteredpeaks(s, scale, k, scalemin = 0.1, varmin = 0.8)
    p = (; env = (;domain = (;min = @SVector[0.0, 0.0], size = SVector{2,Float64}(scale .* size(s.u)))), signals = (;grid = SVector{2,Float64}(size(s.u))))
    return filteredpeaks(s, p, k, scalemin, varmin)
end

function filteredpeaks(s, p::NamedTuple, k, scalemin = 0.1, varmin = 0.8)
    data = s.u
    peaks = Images.findlocalmaxima(data)
    filter!(pk -> !atboundary(data, pk.I), peaks)

    scale = maximum(data)
    u_values = s.u[peaks]
    peaks = peaks[sortperm(u_values, rev=true)]
    peaks = first(peaks, k)


    unique_pks = eltype(peaks)[]
    for pk in peaks
        new_pk = true
        
        # scale test 
        u_val = data[pk.I...]
        if u_val < scalemin * scale
            continue
        end

        # line test
        for qk in unique_pks
            all_above = true
            for x in LinRange(0,1,5)
                Ix = round.(Int64, (1-x) .* pk.I .+ x .* qk.I)
                Ix = mod.(Ix, (axes(data,1), axes(data,2)))
                if data[Ix...] < varmin * u_val
                    all_above = false
                    break
                end
            end
            if all_above
                new_pk = false
                break
            end
        end

        if new_pk
            push!(unique_pks, pk)
        end
    end


    return [ (I = pk.I, pos = indexpos(s,p, pk.I), u = s.u[pk]) for pk in unique_pks]
end

# filteredpeaks(states[end], p, 10)
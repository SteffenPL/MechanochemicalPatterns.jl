using Statistics 

function avg_location(s, ct = 1:2)
    return mean(s.X[i] for i in eachindex(s.X) if s.cell_type[i] in ct)
end

function avg_locations(states, ct = 1:2)
    return [avg_location(s, ct) for s in states]
end

scatter(map(p -> p[1], avg_locations(states, 1:2)))
scatter!(map(p -> p[1], avg_locations(states, 1)))
scatter!(map(p -> p[1], avg_locations(states, 2)))
display(current_figure())
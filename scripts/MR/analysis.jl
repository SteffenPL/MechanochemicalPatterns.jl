using Statistics 

function avg_location(s, ct = 1:2)
    return mean(s.X[i] for i in eachindex(s.X) if s.cell_type[i] in ct)
end

function avg_locations(states, ct = 1:2)
    return [avg_location(s, ct) for s in states]
end

begin 
    Figure()
    yspan = (-80, 20)
    Axis(current_figure()[1,1])
    scatter!(map(p -> p[1], avg_locations(states, 1:2)), color = :white)
    scatter!(map(p -> p[1], avg_locations(states, 1)), color = :magenta)
    scatter!(map(p -> p[1], avg_locations(states, 2)), color = :lightgreen)
    ylims!(yspan...)

    Axis(current_figure()[1,2])
    scatter!(map(p -> p[2], avg_locations(states, 1:2)), color = :white)
    scatter!(map(p -> p[2], avg_locations(states, 1)), color = :magenta)
    scatter!(map(p -> p[2], avg_locations(states, 2)), color = :lightgreen)
    ylims!(yspan...)

    display(current_figure())
end
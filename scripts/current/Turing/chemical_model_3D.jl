using ComponentArrays

p = load_parameters()
s = init_state(p)

function rhs!(dz, z, p, t)
    dz .= 0.0
    laplace!(dz, z, p_ode.D, p_ode.dV, 1.0)
end

xs, ys, zs = LinRange.( p.env.domain.min, p.env.domain.max, p.signals.grid )
p_ode = (; D = p.signals.types.u.D, dV = (xs[2]-xs[1], ys[2]-ys[1], zs[2]-zs[1]))

p_ode = @set p_ode.D = 0.5

ode_prob = ODEProblem(rhs!, s.u, (0.0, p.sim.t_end), p_ode)
ode_integrator = init(ode_prob, ROCK2(); save_everystep=false)

sol = solve(ode_prob, ROCK2())

begin 
    fig = Figure()
    ax = Axis3(fig[1,1], aspect = :data)
    sl = Slider(fig[2,1], range = 1:length(sol), startvalue = 1)
    u_node = @lift sol[$(sl.value)]


    volume!(u_node, colorrange = (0.0, 1.0))
    fig
end
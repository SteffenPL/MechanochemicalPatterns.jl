using GLMakie

N = 100
s = (
    X = rand(N) .+ 0.5,
    fgf = zeros(N),
    wnt = zeros(N)
)

p = (
    sim = (
        L = 2.0,
        dt = 0.0001,
        dx = 2.0/N,
        t_end = 1.0,
    ),
    fgf = (
        D = 0.01,
    ),
    cells = (
        follow = (1.0,),
    )
)

cache = (; fgf = zeros(N))

function step_!(s, p, cache)
    dt = p.sim.dt 
    dx = p.sim.dx 
    (;X, fgf, wnt) = s

    # boundary conditions
    fgf[1] = 1.0
    fgf[end] = 1.0

    # diffusion
    for i in 2:length(fgf)-1
        cache.fgf[i] = p.fgf.D * (s.fgf[i-1] - 2*s.fgf[i] + s.fgf[i+1]) / dx^2
    end

    for i in 2:length(fgf)-1
        s.fgf[i] += dt * cache.fgf[i]

        if isnan(s.fgf[i])
            error("NaN at $i")
        end
    end

    

end


function simulate_!(s, p, cache)
    t = 0.0

    s_ = deepcopy(s)

    sol = [deepcopy(s)]
    while t < p.sim.t_end
        step_!(s_, p, cache)
        t += p.sim.dt
        push!(sol, deepcopy(s_))
    end
    return sol
end



sol = simulate_!(s, p, cache)

begin
    fig = Figure()
    ax = Axis(fig[1, 1])
    sl = Slider(fig[2, 1], range = 1:length(sol), startvalue = 1)
    xgrid = LinRange(0, p.sim.L, N)
    s_obs = @lift sol[$(sl.value)]

    X = @lift $s_obs.X
    FGF = @lift $s_obs.fgf

    scatter!(ax, X, fill(-0.1,N), color = :orange)
    lines!(ax, xgrid, FGF, color = :white)

    ylims!(ax, -0.2, 1.1)
    fig
end
include("base.jl")

rng_seeds = MersenneTwister()

set_theme!(theme_dark())

begin 
    last_seed = rand(rng_seeds, 1:1000)
    Random.seed!(last_seed)

    printstyled("Setting seed to $last_seed\n", color = :green)

    p = load_parameters("$(@__DIR__)/inputs/parameters.toml")
    includet(joinpath("inputs", p.signals.model))

    s = init_state(p)
    cache = init_cache(p, s)

    fig, s_obs = init_plot(s, p, cache; show_polarities = false, bottom_plots = false, show_signals = 1:3, force_scale = 0.3)
    display(fig)
    
    states = [deepcopy(s)]
    s = states[end]
    #
    # @profview simulate(s, p, cache)
    simulate(s, p, cache; callbacks = (update_plot_callback!(fig, s_obs, 0.05),), states = states)
    add_slider!(fig, s_obs, states, p) #view(states,1:600), p)
end
display(fig)
play_animation!(fig, s_obs, states, 10)

# include("analysis.jl")

# fig, s_obs = init_plot(s, p, cache; show_polarities = false, bottom_plots = true)

fig, s_obs = init_plot(states[end], p, cache; show_polarities = false, bottom_plots = false, show_signals = 2:3, force_scale = 0.5)
display(fig)

record(fig, "3D_large_crazy.mp4", 1:1:length(states); framerate = 30) do i
     update_plot!(fig, s_obs, states[i], p)
end




# GLMakie.activate!(ssao=true)
# GLMakie.closeall() # close any open screen

# set_theme!(theme_light())

# begin 

#     fig = Figure(size= ( 1200, 500))
#     ssao = Makie.SSAO(radius = 20.0, blur = 5)
#     cm = (x -> (x, 1.0)).([:lightgreen, :magenta])
#     limits = Rect3( (0,0,0),(1000,1000,800))
#     cm2 = (:thermal, 0.5) #[(:white,0.1), (:lightblue,1.0), (:blue,1.0)]
#     scenekw = (ssao=ssao, )
#     ax0 = LScene(fig[1,1]; scenekw)

#     s = states[5]
#     meshscatter!(Point3f.(s.X), markersize = 6, markerspace = :data, color = s.cell_type, colormap = cm, ssao=true)
#     surface!(LinRange(-500,500,50), LinRange(-500,500,50), fill(-300,50,50), color = sum(s.U.x[2], dims=3)[:,:,1], colormap = cm2)


#     ax1 = LScene(fig[1,2]; scenekw)

#     s = states[30]
#     meshscatter!(Point3f.(s.X), markersize = 6, markerspace = :data, color = s.cell_type, colormap = cm, ssao=true)
#     surface!(LinRange(-500,500,50), LinRange(-500,500,50), fill(-300,50,50), color = sum(s.U.x[2], dims=3)[:,:,1], colormap = cm2)

#     ax2 = LScene(fig[1,3]; scenekw)
#     s = states[70]
#     meshscatter!(Point3f.(s.X), markersize = 6, markerspace = :data, color = s.cell_type, colormap = cm, ssao=true)
#     surface!(LinRange(-500,500,50), LinRange(-500,500,50), fill(-300,50,50), color = sum(s.U.x[2], dims=3)[:,:,1], colormap = cm2)

#     rotate_cam!(ax0.scene, (10.0, 35.0, 0.0) .* π./180)
#     rotate_cam!(ax1.scene, (10.0, 35.0, 0.0) .* π./180)
#     rotate_cam!(ax2.scene, (10.0, 35.0, 0.0) .* π./180)

#     fig
# end

# begin 
#     fig = Figure() 

#     ax = Axis(fig[1,1])

#     s = states[70]
#     heatmap!(sum(s.U.x[3], dims=3)[:,:,1], colormap = :heat, interpolate = true)

#     fig 
# end
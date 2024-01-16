using CSV, DataFrames, Printf

data_dir = "./input/cell_tracking/split_20231212/split/" 
csv_name(rep, type, postfix = "spots"; case = "FWF") = @sprintf "Automation_out/%s%03d_%s_%s.csv" case rep type postfix 
csv_path(rep, type, postfix = "spots"; case = "FWF") = joinpath(data_dir, csv_name(rep, type, postfix; case = case))

df_s1 = CSV.read(csv_path(3, "Distal", "spots"), DataFrame, header = 1, skipto = 5)
mapcols!(collect, df_s1)

df_s2 = CSV.read(csv_path(3, "Proximal", "spots"), DataFrame, header = 1, skipto = 5)
mapcols!(collect, df_s2)

function get_trajs(df, min_frames = 10)
    trajs = Matrix{Float64}[]
    for sdf in groupby(df, :TRACK_ID)
        if nrow(sdf) < min_frames
            continue
        end
        push!(trajs, sdf[:, [:POSITION_X, :POSITION_Y]] |> Matrix |> transpose)
    end 
    return trajs 
end 

trajs = get_trajs(df_1)

function attime(df, frame)
    sdf = filter(row -> row.FRAME == frame, df)
    sort!(sdf, :MEAN_INTENSITY_CH1, rev = true)
    return df[df.FRAME .== frame, [:POSITION_X, :POSITION_Y]] |> Matrix |> transpose
end

using KernelDensity

xs = LinRange(-200, 1500, 200)
ys = LinRange(-200, 1500, 200)

function smooth(pts, xs, ys)
    kd = kde(pts'; npoints = (length(xs), length(ys)), bandwidth = (80.0, 80.0))
    return kd
end

begin
    fig = Figure() 
    sl = Slider(fig[3,1:2], range = 1:290, startvalue = 1)

    X = @lift attime(df_s1, $(sl.value))
    Y = @lift attime(df_s2, $(sl.value))
    #X = @lift $X[:, 1:min(500, size($X, 2))]
    #Y = @lift $Y[:, 1:min(500, size($Y, 2))]
    ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
    scatter!(X, color = :green, markersize = 10.0)
    ax2 = Axis(fig[1,2], xlabel = "x", ylabel = "y", aspect = DataAspect())
    scatter!(Y, color = :magenta, markersize = 10.0)

    kd_X = @lift smooth($X, xs, ys).density
    kd_Y = @lift smooth($Y, xs, ys).density

    ax3 = Axis(fig[2,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(xs, ys, kd_X, colormap = cgrad([:black, :green]))

    ax4 = Axis(fig[2,2], xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(xs, ys, kd_Y, colormap = cgrad([:black, :magenta]))

    for ax in [ax1, ax2, ax3, ax4]
        limits!(ax, 0, 1300, 0, 1300)
    end

    fig 
end
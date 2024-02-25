using CSV, DataFrames, Printf

data_dir = "./input/cell_tracking/split_20231212/split/" 
csv_name(rep, type, postfix = "spots"; case = "FWF") = @sprintf "Automation_out/%s%03d_%s_%s.csv" case rep type postfix 
csv_path(rep, type, postfix = "spots"; case = "FWF") = joinpath(data_dir, csv_name(rep, type, postfix; case = case))

rep = 6

df_s1 = CSV.read(csv_path(rep, "Distal", "spots"), DataFrame, header = 1, skipto = 5)
mapcols!(collect, df_s1)

df_s2 = CSV.read(csv_path(rep, "Proximal", "spots"), DataFrame, header = 1, skipto = 5)
mapcols!(collect, df_s2)

df_t1 = CSV.read(csv_path(rep, "Distal", "tracks"), DataFrame, header = 1, skipto = 5)
mapcols!(collect, df_t1)
df_t2 = CSV.read(csv_path(rep, "Proximal", "tracks"), DataFrame, header = 1, skipto = 5)
mapcols!(collect, df_t2)


function attime(df, frame)
    sdf = filter(row -> row.FRAME == frame, df)
    sort!(sdf, :MEAN_INTENSITY_CH1, rev = true)
    return df[df.FRAME .== frame, [:POSITION_X, :POSITION_Y]] |> Matrix |> transpose
end


using KernelDensity

xs = LinRange(0, 1300, 100)
ys = LinRange(0, 1300, 100)

function smooth(pts, xs, ys)
    kd = kde(pts'; npoints = (length(xs), length(ys)), bandwidth = (50.0, 50.0))
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

    scale = xs[end] / length(xs)
    peaks_X = @lift [pk.pos for pk in filteredpeaks((;u = $kd_X), scale, 10, 0.1, 0.8)]
    peaks_Y = @lift [pk.pos for pk in filteredpeaks((;u = $kd_Y), scale, 10, 0.1, 0.8)]

    ax3 = Axis(fig[2,1], xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(xs, ys, kd_X, colormap = cgrad([:black, :green]))
    scatter!(peaks_X, color = :red, markersize = 10.0)

    ax4 = Axis(fig[2,2], xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(xs, ys, kd_Y, colormap = cgrad([:black, :magenta]))
    scatter!(peaks_Y, color = :red, markersize = 10.0)

    for ax in [ax1, ax2, ax3, ax4]
        limits!(ax, 0, 1300, 0, 1300)
    end

    fig 
end


# analyse distances 

ct = 1
ct_name = ["proximal", "distal"][ct]
(df_spots, df_tracks) = [(df_s1, df_t1), (df_s2, df_t2)][1]

t_end = maximum(df_s1.FRAME)
X_end = attime(df_s1, 290)
X_smooth = smooth(X_end, xs, ys).density
peaks = filteredpeaks((;u = X_smooth), scale, 10, 0.1, 0.95)

min_duration = 20.0

df_t1_long = df_t1[df_t1.TRACK_DURATION .> min_duration, :]

function get_row(df, track_id, frame)
    data = df[df.TRACK_ID .== track_id .&& df.FRAME .== frame, [:POSITION_X, :POSITION_Y]] |> Matrix
    return @SVector[ data[1], data[2] ]
end

tracks_starts = map(row -> get_row(df_s1, row.TRACK_ID, row.TRACK_START), eachrow(df_t1_long))
tracks_ends = map(row -> get_row(df_s1, row.TRACK_ID, row.TRACK_STOP), eachrow(df_t1_long))


pk_idx = [(i, argmin(dist²(pk.pos, tracks_ends[i]) for pk in peaks)) for i in eachindex(tracks_ends)]
dist2 = [dist(tracks_ends[i], peaks[j].pos) for (i, j) in pk_idx]
dist1 = [dist(tracks_starts[i], peaks[j].pos) for (i, j) in pk_idx]



id = "rep$(rep)"
mkpath("scripts/current/RT/outputs/data/$(id)")
logscale = false 
nbins = 40 
begin 
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "initial distance", ylabel = "distance delta towards peaks", yscale = logscale ? log10 : identity)
    logscale && ylims!(ax,(0.5, 2))
    data = logscale ? dist2 ./ dist1 : dist2 .- dist1
    scatter!(dist1, data, color = :gray, markersize = 5.0)

    if !logscale
        bins, data_means = digitalize(dist1, dist2 .- dist1, nbins)
        lines!(bins, data_means, color = :red, linewidth = 5.0)
    end
    hlines!([logscale ? 1.0 : 0.0], color = :orange, linewidth = 5)
    save("""scripts/current/RT/outputs/data/$(id)/distance_$(logscale ? "log" : "lin")_delta_$(id)_$(ts[1])_$(ts[2]).png""", fig)
    display(fig)
end
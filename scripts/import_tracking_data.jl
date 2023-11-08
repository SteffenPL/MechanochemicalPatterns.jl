using CSV, DataFrames, Glob, TOML
config = TOML.parsefile("scripts/config.toml")


csv_files = readdir(glob"*/FW_Day1_Cycle_01/*.csv", config["data_dir"])
df = CSV.read(csv_files[2], DataFrame, skipto = 5)

# clean up data 
sort!(df, :POSITION_T)
last_frame = maximum(df.FRAME)
filter!( row -> row.FRAME < last_frame, df)

using GLMakie
begin 
    fig = Figure()
    ax = Axis3(fig[1,1])

    # for label in unique(df[:, :ID])
    #     df_label = filter(row -> row.ID == label, df)
    #     lines!(ax, df_label[:, :POSITION_X], df_label[:, :POSITION_Y], df_label[:, :POSITION_Z], label=label)
    # end

    for (label, df_label) in pairs(groupby(df, :TRACK_ID))
        lines!(ax, df_label.POSITION_X, df_label.POSITION_Y, df_label.POSITION_Z, label=label)
    end

    df_init = filter(row -> row.POSITION_T == 0, df)
    scatter!(ax, df_init.POSITION_X, df_init.POSITION_Y, df_init.POSITION_Z, markersize = 50^2, color = :black, label = "init")
    
    fig 
end
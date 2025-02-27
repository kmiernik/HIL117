using CSV
using DataFrames
using GLMakie
using TOML
using StatsBase


pub_theme = Theme(
               Axis = (
                   xlabelsize = 24, 
                   xticklabelsize = 24, 
                   xminorticks = IntervalsBetween(5), 
                   xminorticksvisible = true, 
                   xgridvisible = false,
                   ylabelsize = 24, 
                   yticklabelsize = 24, 
                   ygridvisible = false, 
                   yminorticks = IntervalsBetween(5), 
                   yminorticksvisible = true,
               ),
               Legend = (
                   labelsize = 24,
               )
                        )

function load_data(ref_ch=6000)
    rf = DataFrame(CSV.File("hil117_runs.csv"))

    runs = Int64[]
    ge_ch = Dict{Int64, Vector{Float64}}()
    ge_t = Dict{Int64, Vector{Float64}}()
    for i in 1:2:32
        ge_ch[i] = Float64[]
        ge_t[i] = Float64[]
    end

    files = readdir("../config/", join=true)
    for file in files
        if !startswith(basename(file), "config_") || !endswith(basename(file), ".toml")
            continue
        end
        run_number = parse(Int64, split(basename(file), ['_', '.'])[2])
        push!(runs, run_number)
        config = TOML.parsefile(file)
        for loc in 1:2:32
            cal = config["label"]["$loc"]["cal"]
            push!(ge_ch[loc], ref_ch * cal[2] + cal[1])
            push!(ge_t[loc], config["label"]["$loc"]["dt"])
        end
    end

    df = DataFrame(run=runs)
    for loc in [1:2:32; 17:2:32]
        df[!, "E_$loc"] = ge_ch[loc]
        df[!, "t_$loc"] = ge_t[loc]
    end
    innerjoin(df, rf, on=:run)
end


function plot_shifts(df)
    set_theme!(Theme(
               Axis = (
                   xlabelsize = 20, 
                   xticklabelsize = 20, 
                   xminorticks = IntervalsBetween(5), 
                   xminorticksvisible = true, 
                   xgridvisible = false,
                   ylabelsize = 20, 
                   yticklabelsize = 20, 
                   ygridvisible = false, 
                   yminorticks = IntervalsBetween(5), 
                   yminorticksvisible = true,),
                    Legend = (
                        labelsize = 20,
                  )))
    df = df[(df.run .> 27) .& (df.run .< 200), :]
    #=
    fig = Figure(size=(800, 800))
    for loc in [1:2:14; 17:2:32]
        c = round(Int, loc/2, RoundUp)
        ix = round(Int, (c-1)/4, RoundDown)+1
        iy = (c-1)%4+1
        ax = Axis(fig[ix, iy], title="$loc", xticklabelrotation=45.0) 
        ym = df[df.run .== 28, "t_$loc"]
        scatter!(ax, df.start, df[!, "t_$loc"] .- ym)
        ylims!(ax, -20, 40)
    end
    =#
    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1], xticklabelrotation=45.0, ylabel="ΔT (ns)")
    for loc in [1:2:14; 17:2:32]
        ym = df[df.run .== 28, "t_$loc"]
        scatterlines!(ax, df.start, df[!, "t_$loc"] .- ym, linestyle=:dash)
        ylims!(ax, -20, 40)
    end
    fig

end


function plot_energies(df)
    set_theme!(Theme(
               Axis = (
                   xlabelsize = 20, 
                   xticklabelsize = 20, 
                   xminorticks = IntervalsBetween(5), 
                   xminorticksvisible = true, 
                   xgridvisible = false,
                   ylabelsize = 20, 
                   yticklabelsize = 20, 
                   ygridvisible = false, 
                   yminorticks = IntervalsBetween(5), 
                   yminorticksvisible = true,),
                    Legend = (
                        labelsize = 20,
                  )))
    df = df[(df.run .> 27) .& (df.run .< 200), :]
    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1], xticklabelrotation=45.0, ylabel="ΔE (keV)", 
              title="Ch 6000")
    dys = Float64[]
    for loc in [1:2:14; 17:2:32]
        ym = df[df.run .== 28, "E_$loc"]
        dy = df[!, "E_$loc"] .- ym
        scatterlines!(ax, df.start, dy, linestyle=:dash)
        push!(dys, dy...)
        ylims!(ax, -3, 3)
    end
    h = fit(Histogram, dys, -3:0.1:3)
    ax2 = Axis(fig[1, 2], yaxisposition=:right, xreversed=true, 
               ylabel="ΔE (keV)")
    stairs!(ax2, h.weights, h.edges[1][1:end-1])
    lines!(ax2, cumsum(h.weights) ./ sum(h.weights) .* 400, 
           h.edges[1][1:end-1], color=:black, linestyle=:dash)
    ylims!(ax2, -3, 3)
    xlims!(400, 0)
    colsize!(fig.layout, 2, Relative(1/3))
    fig

end

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
                   xlabelsize = 12, 
                   xticklabelsize = 12, 
                   xminorticks = IntervalsBetween(5), 
                   xminorticksvisible = true, 
                   xgridvisible = false,
                   ylabelsize = 12, 
                   yticklabelsize = 12, 
                   ygridvisible = false, 
                   yminorticks = IntervalsBetween(5), 
                   yminorticksvisible = true,),
                    Legend = (
                        labelsize = 12,
                  )))
    df = df[df.run .> 27, :]
    #=
    fig = Figure(size=(800, 800))
    for loc in 1:2:32
        c = round(Int, loc/2, RoundUp)
        ix = round(Int, (c-1)/4, RoundDown)+1
        iy = (c-1)%4+1
        ax = Axis(fig[ix, iy], title="$loc")
        x1 = findfirst(x->x>=28, runs)
        ym = median(ge_ch[loc][x1:end])
        stem!(ax, runs, ge_ch[loc] .- ym)
        xlims!(ax, 28, 50)
        ylims!(ax, -2, 2)
    end
    fig
    =#
    fig = Figure(size=(800, 800))
    for loc in [1:2:14; 17:2:32]
        c = round(Int, loc/2, RoundUp)
        ix = round(Int, (c-1)/4, RoundDown)+1
        iy = (c-1)%4+1
        ax = Axis(fig[ix, iy], title="$loc", xticklabelrotation=45.0) 
        ym = mean(df[!, "t_$loc"])
        stem!(ax, df.start, df[!, "t_$loc"] .- ym)
        ylims!(ax, -20, 20)
    end
    fig

end


function plot_energies(df)
    set_theme!(Theme(
               Axis = (
                   xlabelsize = 12, 
                   xticklabelsize = 12, 
                   xminorticks = IntervalsBetween(5), 
                   xminorticksvisible = true, 
                   xgridvisible = false,
                   ylabelsize = 12, 
                   yticklabelsize = 12, 
                   ygridvisible = false, 
                   yminorticks = IntervalsBetween(5), 
                   yminorticksvisible = true,),
                    Legend = (
                        labelsize = 12,
                  )))
    df = df[df.run .> 27, :]
    fig = Figure(size=(800, 800))
    for loc in [1:2:14; 17:2:32]
        c = round(Int, loc/2, RoundUp)
        ix = round(Int, (c-1)/4, RoundDown)+1
        iy = (c-1)%4+1
        ax = Axis(fig[ix, iy], title="$loc", xticklabelrotation=45.0) 
        ym = median(df[!, "E_$loc"])
        stem!(ax, df.start, df[!, "E_$loc"] .- ym)
        ylims!(ax, -2, 2)
    end
    fig

end

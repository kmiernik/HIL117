using HDF5
#using GLMakie
using CairoMakie
using DataFrames
using LaTeXStrings


pub_theme = Theme(
        Axis = (
            xlabelsize = 28, 
            xticklabelsize = 24, 
            xminorticks = IntervalsBetween(5), 
            xminorticksvisible = true, 
            xgridvisible = false,
            ylabelsize = 28, 
            yticklabelsize = 24, 
            ygridvisible = false, 
            yminorticks = IntervalsBetween(5), 
            yminorticksvisible = true,
        ),
        Legend = (
            labelsize = 24,
        )
                 )

function plot_gamma_alpha()
    fin = h5open("../run039_edf.h5", "r")
    g = vec(read(fin["g"]))
    ga = vec(read(fin["ga"]))
    g0a = vec(read(fin["g0a"]))
    g2a = vec(read(fin["g2a"]))
    gp0a = vec(read(fin["gp0a"]))

    tf = DataFrame(
                   gammas=Vector{Float64}[],
                label=LaTeXString[],
                marker=Union{Symbol, Char}[],
                color=Symbol[])

    push!(tf, [[286, 416], L"$^{132}_{62}${Sm}", :ltriangle, :green])
    push!(tf, [[163, 316, 418], L"$^{134}_{62}${Sm}", :rtriangle, :blue])
    push!(tf, [[254.9, 431.4], L"$^{136}_{62}${Sm}", :utriangle, :teal])
    push!(tf, [[346.7, 544.4], L"$^{138}_{62}${Sm}", :dtriangle, :brown])
    push!(tf, [[213.2, 397.7, 521.4], L"$^{132}_{60}${Nd}", :star4, :red])
    push!(tf, [[294.2], L"$^{134}_{60}${Nd}", :star6, :slateblue])
    push!(tf, [[207.1, 399.6, 550.4], L"$^{128}_{58}${Ce}", :circle, :orange])

    fig = Figure(size=(1000, 600))
    with_theme(pub_theme) do
        ax1 = Axis(fig[1, 1]; xlabel="E (keV)", ylabel="Counts", 
                   xticks=100:100:600, xminorticks=IntervalsBetween(10))
        norm = sum(ga[100:600, 1]) / sum(g[100:600, 1])
        stairs!(ax1, 1:1:3096, g[:, 1] .* norm, color=:lightsteelblue,
                label="γ")
        stairs!(ax1, 1:1:3096, ga[:, 1], color=:black,
                  label="α-gated γ")
        xlims!(100, 660)
        ylims!(1e4, 2e4)
        for row in eachrow(tf)
            ypos = Float64[]
            for yi in round.(Int64, row.gammas)
                push!(ypos, maximum(ga[yi-3:yi+3]) * 1.01)
            end
            scatter!(ax1, row.gammas, ypos,
                     markersize=14,
                     marker=row.marker, label=row.label, color=row.color)
        end
        axislegend(ax1, position=:rt, fontsize=20)

        #=
        ax2 = Axis(fig[2, 1]; xlabel="E (keV)", ylabel="Counts", 
                  title="0α-p-gated γ")
        stairs!(ax2, 1:1:3096, gp0a[:, 1], color=:black)
        xlims!(100, 700)
        ylims!(0, 400)
        for row in eachrow(tf)
            ypos = Float64[]
            for yi in round.(Int64, row.gammas)
                push!(ypos, maximum(g2a[yi-3:yi+3]) * 1.01)
            end
            scatter!(ax2, row.gammas, ypos,
                     markersize=14,
                     marker=row.marker, label=row.label, color=row.color)
        end
        axislegend(ax2, position=:rt, fontsize=20)
        =#
    end
    save("temp.pdf", fig)
    fig
end

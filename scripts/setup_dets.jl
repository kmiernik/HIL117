"""
    fit_ban(dia_spec, c_lim, E_lim, loc)        
    
    Fits banana shape to diamant at location loc,
    returns params of banana, and saves a plot
"""
function fit_ban(dia_spec, c_lim, E_lim, loc)        
    Ai = Float64[]
    Pi = Float64[]
    Ei = Float64[]
    fig = Figure(size=(800, 800))
    fi = 1
    t = 0.00:0.01:1.27
    b(x, p) = @. p[1] / x^2 + p[2]
    for ie in 1:size(dia_spec)[1]
        if maximum(dia_spec[ie, :]) > c_lim && (100 * (ie - 1) < E_lim)
            ix = round(Int, (fi-1)/6, RoundDown)+1
            iy = (fi-1)%6+1
            fi += 1
            ax = Axis(fig[ix, iy]; title="$(100*ie+50)",
                                   yticklabelsvisible=false)
            tg = t[argmax(dia_spec[ie, :])]
            guess = [sum(dia_spec[ie, :]) * 0.015, tg, 0.04, 
                         sum(dia_spec[ie, :]) * 0.003, tg+0.12, 0.04]
            gf = curve_fit(ngauss, t, dia_spec[ie, :], guess)
            scatter!(ax, t, dia_spec[ie, :], marker=:cross, markersize=5)
            lines!(ax, t, ngauss(t, guess), color=:green, linestyle=:dot)
            lines!(ax, t, ngauss(t, gf.param), color=:red, linestyle=:solid)
            xlims!(ax, 0.3, 1.0)
            ylims!(ax, 0, maximum(dia_spec[:, :]))
            push!(Ei, 100 * ie + 50)
            if gf.param[2] < gf.param[5]
                push!(Ai, gf.param[2])
                push!(Pi, gf.param[5])
            else
                push!(Ai, gf.param[5])
                push!(Pi, gf.param[2])
            end
        else
            continue
        end
    end
    ix = round(Int, fi/6, RoundDown) + 1
    ax = Axis(fig[ix+1:ix+3, 1:6])
    heatmap!(ax, 50:100:100*128+50, 0:0.01:1.28, dia_spec[:, :])
    scatter!(ax, Ei, Ai, marker=:cross, color=:red, markersize=5)
    scatter!(ax, Ei, Pi, marker=:xcross, color=:gray, markersize=5)
    bf = curve_fit(b, Ei, Ai, [30000.0, 0.6])
    lines!(ax, 550:100:E_lim+100, b(550:100:E_lim+100, bf.param), 
           color=:orange, linestyle=:dash)
    xlims!(100, E_lim + 1000)
    ylims!(0, 1.2)
    name = @sprintf("b_%03d.png", loc)
    save(name, fig)
    bf.param
end

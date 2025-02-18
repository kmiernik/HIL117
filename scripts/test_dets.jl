#using HIL117
using GLMakie
using Printf
using MTools
using LsqFit
using TOML


"""
    Calibrate using Eu-152 data (energy calibration)

    Data are loaded from output of scan_raw_ge (16x16384 array)
"""
function cal_energy(filename, configfile, prefix)

    config = TOML.parsefile(configfile)
    spectra = zeros(Int64, 16, 32768)
    fin = open(filename, "r")
    read!(fin, spectra)
    close(fin)

    eu = [
    121.7817 28.37 
    244.6975 7.53 
    295.9390 0.44
    344.2785 26.59
    367.7891 0.859
    411.1165 2.238 
    444.0    3.125
    778.9045 12.97
    867.378  4.214
    964.1    14.63
    1085.836 10.13
    1089.737 1.731 
    1112.074 13.54 
    1212.948 1.412
    1299.140 1.626
    1408.011 20.85
    #1457.643 0.497
    #1460.820 0.000
    2614.511 0.000
    ]

    new_config = copy(config)
    dch = 90
    for ige in 1:16
        label = (ige - 1) * 2 + 1
        if !config["label"]["$label"]["valid"]
            continue
        end
        pinit = config["label"]["$label"]["rcal"]
        
        chfit = Float64[]
        cherr = Float64[]
        ip = 0
        fig = Figure(size=(1000, 1000), title="$label")
        n_peaks = size(eu)[1]
        guess = []
        while ip < n_peaks
            ip += 1
            x0 = round(Int64, (eu[ip, 1] - pinit[1]) / pinit[2])
            x1 = -1
            if eu[ip, 1] == 1085.836 || eu[ip, 1] == 1457.643
                x1 = round(Int64, (eu[ip+1, 1] - pinit[1]) / pinit[2])
                guess = [maximum(spectra[ige, x0-dch:x0+dch]) * 5,
                           argmax(spectra[ige, x0-dch:x0+dch])+x0-dch,
                            4.0, 
                            maximum(spectra[ige, x0-dch:x0+dch]) * 1,
                            argmax(spectra[ige, x0-dch:x0+dch])+x0-dch+20, 
                            4.0,
                            spectra[ige, x0-dch], 0.0 ]
                lower=[0, x0-dch, 1.0,
                       0, x0-dch, 1.0, 
                       -Inf, -Inf]
                upper=[sum(spectra[ige, x0-dch:x0+dch]), x0+dch, 10.0,
                        sum(spectra[ige, x0-dch:x0+dch]), x0+dch, 10.0,
                        Inf, Inf]
                pf = curve_fit(ngausslin, x0-dch:x0+dch, 
                           spectra[ige, x0-dch:x0+dch], 
                           guess, lower=lower, upper=upper)
            else
                guess = [maximum(spectra[ige, x0-dch:x0+dch]) * 5,
                         argmax(spectra[ige, x0-dch:x0+dch])+x0-dch,
                         4.0, 
                         spectra[ige, x0-dch], 0.0 ]
                pf = curve_fit(gausslin, x0-dch:x0+dch, 
                           spectra[ige, x0-dch:x0+dch], 
                           guess)
            end
            pferr = pf.param
            try
                pferr = standard_errors(pf)
            catch
            end
            push!(chfit, pf.param[2])
            push!(cherr, pf.param[3])
            if x1 > 0
                push!(chfit, pf.param[5])
                push!(cherr, pf.param[6])
            end

            ix = round(Int, (ip-1)/5, RoundDown)+1
            iy = (ip-1)%5+1
            ax = Axis(fig[ix, iy]; title="$(eu[ip, 1])")
            stairs!(ax, x0-dch-20:x0+dch+20, 
                    vec(spectra[ige, x0-dch-20:x0+dch+20]))
            if x1 > 0
                lines!(ax, x0-dch:0.1:x0+dch, ngausslin(x0-dch:0.1:x0+dch, 
                                                pf.param), color="orange")
                lines!(ax, x0-dch:0.1:x0+dch, 
                        ngausslin(x0-dch:0.1:x0+dch, guess), color="pink")
            else
                lines!(ax, x0-dch:0.1:x0+dch, gausslin(x0-dch:0.1:x0+dch, 
                                                pf.param), color="red")
                lines!(ax, x0-dch:0.1:x0+dch, gausslin(x0-dch:0.1:x0+dch, 
                                    guess), color="pink")
            end
            if x1 > 0
                ip += 1
            end
        end

        hyplin(x, p) = @. p[1] + p[2] * x + p[3] / x

        sewpoint = Dict( 1 => 5000.0,
                         3 => 4000.0,
                         5 => 5700.0,
                         7 => 4000.0,
                         9 => 4000.0,
                         11 => 4950.0,
                         13 => 4500.0,
                         17 => 3000.0,
                         19 => 5500.0,
                         21 => 4000.0,
                         23 => 5000.0,
                         25 => 5500.0,
                         27 => 4000.0,
                         29 => 4000.0,
                         31 => 3500.0)
        
        lf = curve_fit(lin, chfit, eu[:, 1], pinit[1:2])
        new_config["label"]["$label"]["rcal"] = lf.param
        qf = curve_fit(quad, chfit, eu[:, 1], [pinit; 0.0] )
        llf = curve_fit(linlin, chfit, eu[:, 1], 
                        [pinit; pinit; sewpoint[label]])
        lqf = curve_fit(linquad, chfit, eu[:, 1], 
                        [pinit; pinit; 0.0; sewpoint[label]])
        a7 = ((lqf.param[3] - lqf.param[1]) / lqf.param[6]^2 
              + (lqf.param[4] - lqf.param[2]) / lqf.param[6]
              + lqf.param[5])
        new_config["label"]["$label"]["cal"] = [lqf.param; a7]
        #qqf = curve_fit(quadquad, chfit, eu[:, 1], [pinit; 0.0; pinit; 0.0])
        println(label)
        @printf("\tLin     χ²: %.3f, m: %.3f\n", 
                sum(lf.resid .^ 2 ./ cherr) / (n_peaks - 2), 
                maximum(abs.(lf.resid)))
        @printf("\tQuad    χ²: %.3f, m: %.3f\n", 
                sum(qf.resid .^ 2 ./ cherr) / (n_peaks - 3), 
                maximum(abs.(qf.resid)))
        @printf("\tLinLin  χ²: %.3f, m: %.3f\n", 
                sum(llf.resid .^ 2 ./ cherr) / (n_peaks - 4),
                maximum(abs.(llf.resid)))
        @printf("\tLinQuad χ²: %.3f, m: %.3f\n", 
                sum(lqf.resid .^ 2 ./ cherr) / (n_peaks - 5),
                maximum(abs.(lqf.resid)))
        #@printf("\tQuadQuad χ²: %.3f, m: %.3f\n", 
        #        sum(qqf.resid .^ 2 ./ cherr) / (n_peaks - 6), 
        #        maximum(abs.(qqf.resid)))
        
        axl = Axis(fig[5, 1]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    yticks=-0.5:0.25:0.5,
                    limits=(0, 3000, -0.5, 0.5), title="Lin")
        axq = Axis(fig[5, 2]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    yticks=-0.5:0.25:0.5,
                    limits=(0, 3000, -0.5, 0.5), title="Quad")
        axll = Axis(fig[5, 3]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    yticks=-0.5:0.25:0.5,
                    limits=(0, 3000, -0.5, 0.5), title="LinLin")
        axlq = Axis(fig[5, 4]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    yticks=-0.5:0.25:0.5,
                    limits=(0, 3000, -0.5, 0.5), title="LinQuad")
        axqq = Axis(fig[5, 5]; ylabel="σ (%)", xlabel="E (keV)", 
                    limits=(0, 3000, nothing, nothing))
        stem!(axl, eu[:, 1], lf.resid, color="black", marker=:circle)
        stem!(axq, eu[:, 1], qf.resid, color="blue", marker=:utriangle)
        stem!(axll, eu[:, 1], llf.resid, color="red", marker=:diamond)
        stem!(axlq, eu[:, 1], lqf.resid, color="green", marker=:star5)
        scatter!(axqq, eu[:, 1], 100.0 .* abs.(cherr) ./ eu[:, 1], 
                    color="violet", marker=:star6,
                    label="σ")
        save(@sprintf("%s_%02d.png", prefix, label), fig)
        #=
        fig2 = Figure(size=(800, 600))
        ax = Axis(fig2[1, 1], xticks=0:5000:20000,
                    xminorticks=IntervalsBetween(10), xminorticksvisible=true,
                    xminorgridvisible=true)
        scatter!(ax, chfit, lf.resid, color="black", marker=:circle,
                    label="Lin")
        vlines!(ax, [sewpoint[label]], color="red", linestyle=:dash)
        save(@sprintf("lin_%02d.png", label), fig2)
        =#
    end
    new_config
end


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


function fine_cal_test(E_loc, configfile::String; args...)
    config = TOML.parsefile(configfile)
    fine_cal_test(E_loc, config; args...)
end


function residual(p, input)
    yfit = gausslin(input[1], p)
    return (yfit .- input[2]) .^ 2
end


"""
    fine_cal_test(E_loc, configfile; doplot=false)

    Testing fine cal based on two lines in raw spectrum
"""
function fine_cal_raw(E_loc, config::Dict; prefix="fcal")
    new_config = copy(config)

    fig = Figure(size=(1000, 1000))
    for loc in 1:2:32
        if !config["label"]["$loc"]["valid"]
            continue
        end

        x = collect(1:32768)
        y = E_loc[loc, :]

        x511 = argmax(y[100:5000]) + 100
        ab = 511 / x511
        ch1257 = round(Int, 1257 / ab)

        x1257 = argmax(y[ch1257-100:ch1257+100]) + ch1257 - 100

        ige = round(Int, loc/2, RoundUp)
        ix = round(Int, (ige-1)/4, RoundDown)+1
        iy = (ige-1)%4+1
        ax = Axis(fig[ix, iy]; title="$loc")

        stairs!(ax, x[2000:9000], y[2000:9000])

        pf = curve_fit(gausslin, x[x511-50:x511+50], y[x511-50:x511+50], 
                       [5.0*y[x511], x511, 10.0, y[x511-50], 0.0])
        ch511 = pf.param[2]
        lines!(ax, x[x511-50:x511+50],
                gausslin(x[x511-50:x511+50], pf.param), color="orange")
        guess = [5.0*y[x1257], x1257, 10.0, y[x1257-50], 0.0]
        pf = curve_fit(gausslin, x[x1257-50:x1257+50], y[x1257-50:x1257+50], 
                       guess)
        lines!(ax, x[x1257-50:x1257+50],
                gausslin(x[x1257-50:x1257+50], pf.param), color="orange")
        lines!(ax, x[x1257-50:x1257+50],
                gausslin(x[x1257-50:x1257+50], guess), color="pink")
        ch1257 = pf.param[2]
        ar = (1256.69 - 511.0) / (ch1257 - ch511)
        br = 1256.69 -  ar * ch1257
        println(loc, " ", ch511, " ", ch1257, " ", ar, " ", br)
                        
        
        #=
        Efit = Float64[]
        sfit = Float64[]
        ip = 0
        n_peaks = size(lines)[1]
        guess = []
        fig = Figure(size=(1000, 1000), title="$loc")
        binsize = 1
        #=
        x = rebin(collect(1:8192).-0.5, binsize)[:, 1]
        y = rebin(E_loc[loc, :], binsize)[:, 1]
        =#
        for ip in 1:n_peaks
            E0 = round(Int64, lines[ip, 1] / binsize)
            dE = round(Int64, lines[ip, 2] / binsize)
            xm = argmax(y[E0-dE:E0+dE]) + E0 - dE - 1
            guess = [2 * (maximum(y[E0-dE:E0+dE]) - minimum(y[E0-dE:E0+dE])),
                     x[xm], 1.5, minimum(y[E0-dE:E0+dE]), 0.0]

            pf = curve_fit(gausslin, x[E0-dE:E0+dE], y[E0-dE:E0+dE], 
                        guess, 
                        lower=[0, x[E0-dE-1], 0.8, -Inf, -Inf],
                        upper=[sum(y[E0-dE:E0+dE]), x[E0+dE+1], 3.0, Inf, Inf])
            pferr = pf.param
            try
                pferr = standard_errors(pf)
            catch
            end
            push!(Efit, pf.param[2])
            push!(sfit, pf.param[3])

            ix = round(Int, (ip-1)/4, RoundDown)+1
            iy = (ip-1)%4+1
            ax = Axis(fig[ix, iy]; title="$(lines[ip, 1])")
            stairs!(ax, x[E0-dE-10:E0+dE+10], y[E0-dE-10:E0+dE+10])
            lines!(ax, x[E0-dE:E0+dE],
                   ngausslin(x[E0-dE:E0+dE], pf.param), color="orange")
            lines!(ax, x[E0-dE:E0+dE],
                   ngausslin(x[E0-dE:E0+dE], guess), color="pink")
        end

        lf = curve_fit(lin, Efit, lines[:, 1], [0.0, 1.0] )
        qf = curve_fit(quad, Efit, lines[:, 1], [0.0, 1.0, 0.0] )
        new_config["label"]["$loc"]["fcal"] = qf.param
        
        axl = Axis(fig[3, 1:2]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    yticks=-0.5:0.25:0.5,
                    limits=(0, 3000, -1.0, 1.0), title="Lin")
        axq = Axis(fig[3, 3:4]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    yticks=-0.5:0.25:0.5,
                    limits=(0, 3000, -1.0, 1.0), title="Quad")
        stem!(axl, lines[:, 1], lf.resid, color="black", marker=:circle)
        stem!(axq, lines[:, 1], qf.resid, color="blue", marker=:utriangle)
        =#
        new_config["label"]["$loc"]["rcal"] = [br, ar]
        save(@sprintf("%s_%02d.png", prefix, loc), fig)
    end
    new_config
end

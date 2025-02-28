"""
    cal_energy(filename, configfile, prefix)

    Depraceted since linear online calibration is used

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


"""
    find_shifts(t_loc, config, prefix; dt_min=-1000, dt=1, dt_max=1000)
    
    Based on t_loc table (dt_min:dt:dt_max) calculates t_loc value for
    each channel, returns config Dict with updated values

    saves plots for each channel

"""
function find_shifts(t_loc, config, prefix; dt_min=-1000, dt=1, dt_max=1000,
                                            thr=30)
    last_label = size(t_loc)[1]
    new_config = copy(config)
    for loc in 1:last_label
        if sum(t_loc[loc, :]) == 0 || !config["label"]["$loc"]["valid"]
            continue
        end
        y = vec(t_loc[loc, :])
        t = dt_min:dt:dt_max
        binsize = 1
        if config["label"]["$loc"]["type"] !== "NEDA"
            binsize = 5
            y = vec(rebin(y[2:end], binsize))
            t = vec(rebin(collect(dt_min+1:dt:dt_max), binsize))
        end
        dy = y[2:end] .- y[1:end-1]
        
        tp = -1
        tm = -1
        zc = -1
        thr = abs(minimum(dy) / 4)
        xm = findfirst(x->x>maximum(dy)/2, dy)
        xd = 0
        if config["label"]["$loc"]["type"] == "GE"
            xd = round(Int, 80 / binsize)
        elseif config["label"]["$loc"]["type"] == "BGO"
            xd = round(Int, 80 / binsize)
        elseif config["label"]["$loc"]["type"] == "NEDA"
            xd = round(Int, 10 / binsize)
        elseif config["label"]["$loc"]["type"] == "DIAMANT"
            xd = round(Int, 40 / binsize)
        end
        for i in xm:xm+xd
            if dy[i] > 0
                tp = i
            elseif dy[i] <= 0 && tp > 0
                tm = i
                a = (dy[tm] - dy[tp]) / (t[tm] - t[tp])
                b = dy[tm] - t[tm] * a
                zc = -b/a
                break
            end
        end
        if zc < 0
            zc = t[argmax(y)]
        end

        new_config["label"]["$loc"]["dt"] = zc

        fig = Figure(size=(600, 400))
        ax1 = Axis(fig[1, 1]; title="$loc")
        ax2 = Axis(fig[2, 1]; title="$loc")
        linkxaxes!(ax1, ax2)
        hlines!(ax2, [0.0], color=:black)

        scatter!(ax1, t, y)
        vlines!(ax1, [zc], color=:red, linestyle=:dot)

        scatter!(ax2, t[1:end-1], dy)
        if tp > 0 && tm > 0
            scatterlines!(ax2, [t[tp], t[tm]], [dy[tp], dy[tm]], 
                     marker=:cross, color=:red)
        else
            scatter!(ax1, [zc], [maximum(y)],
                     marker=:utriangle, color=:red)
        end

        vlines!(ax2, [zc], color=:red, linestyle=:dot)
        vlines!(ax2, [t[xm], t[xm+xd]], color=:gray, linestyle=:dot)

        if config["label"]["$loc"]["type"] == "Ge"
            xlims!(ax2, zc-100, zc+100)
        elseif config["label"]["$loc"]["type"] == "BGO"
            xlims!(ax2, zc-100, zc+100)
        elseif config["label"]["$loc"]["type"] == "DIAMANT"
            xlims!(ax2, zc-100, zc+100)
        elseif config["label"]["$loc"]["type"] == "NEDA"
            xlims!(ax2, zc-20, zc+20)
        end
        save(@sprintf("%s_%03d.png", prefix, loc), fig)
    end
    new_config
end


"""
    find_period(t_loc, config, prefix; period_loc=35, dt_min=-1000, dt=1, 
                                        dt_max=1000)

    Find beam period based on NEDA detector timing

    512 ns samples in a "safe" range (wihtout main peak)
    detector next to ref (default 34, first NEDA) 
    shouldn't have scattering "echos"
    finding minima in the safe range will give off-beam period
    
    returns config with updated value, saves plot
"""
function find_period(t_loc, config, prefix; period_loc=35, dt_min=-1000, dt=1, 
        t_safe=200, dt_max=1000)
    
    i1 = round(Int64, (t_safe-dt_min) / dt)
    i2 = i1 + round(Int64, 512 / dt) - 1
    tm = collect(dt_min:dt:dt_max)[i1:i2]
    ym = copy(t_loc[period_loc, i1:i2])
    ym4 = vec(rebin(ym, 4))

    ic = 4
    localmins = Float64[]
    while ic < 124
        ic += 1
        if ((ym4[ic-3] > ym4[ic-2] > ym4[ic-1] > ym4[ic]) &&
            (ym4[ic] < ym4[ic+1] < ym4[ic+2] < ym4[ic+3]))
            fi = ic * 4 - 2
            qf = curve_fit(quad, tm[fi-7:fi+7], ym[fi-7:fi+7], [0.0, 1.0, 0.0])
            push!(localmins, -qf.param[2] / (2 * qf.param[3]))
            ic += 3
        end
    end
    new_config = copy(config)
    new_config["spectra"]["beam_period"] = mean(localmins[2:end] 
                                                - localmins[1:end-1])
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1]; xlabel="t (ns)", 
            title=@sprintf("%.3f ns", new_config["spectra"]["beam_period"]))
    scatter!(ax, tm, ym)
    vlines!(ax, localmins, color=:red, linestyle=:dash)
    save("$prefix.png", fig)

    new_config
end


function fine_cal(E_loc, configfile::String, prefix)
    config = TOML.parsefile(configfile)
    fine_cal(E_loc, config, prefix)
end


"""
    fine_cal(E_loc, configfile, prefix)

    Take E_loc table and perform fine calibration for the run
    Returns new config with fcal fields for Ge

    Depraceted since linear online calibration is used
"""
function fine_cal(E_loc, config::Dict, prefix)
    new_config = copy(config)

    lines = [
             294.17 3
             346.71 5
             373.75 5
             409.20 5
             511.00 5
             552.05 4
             #657.76 4
             #1256.69 5
             1293.56 15
             #1299.91 5
            ]
    for loc in 1:2:32
        if !config["label"]["$loc"]["valid"]
            continue
        end
        
        Efit = Float64[]
        sfit = Float64[]
        ip = 0
        fig = Figure(size=(1000, 1000), title="$loc")
        n_peaks = size(lines)[1]
        guess = []
        for ip in 1:n_peaks
            E0 = round(Int64, lines[ip, 1])
            dE = round(Int64, lines[ip, 2])
            guess = [maximum(E_loc[loc, E0-dE:E0+dE]),
                        argmax(E_loc[loc, E0-dE:E0+dE])+E0-dE,
                        1.0, 
                        E_loc[loc, E0-dE], 0.0 ]
            pf = curve_fit(gausslin, E0-dE:E0+dE, 
                        E_loc[loc, E0-dE:E0+dE], 
                        guess, 
                        lower=[0, E0-dE-1, 0.8, -Inf, -Inf],
                upper=[sum(E_loc[loc, E0-dE:E0+dE]), E0+dE+1, 3.0, Inf, Inf])
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
            stairs!(ax, E0-dE-10:E0+dE+10, 
                    vec(E_loc[loc, E0-dE-10:E0+dE+10]))
            lines!(ax, E0-dE:0.1:E0+dE, gausslin(E0-dE:0.1:E0+dE, 
                                            pf.param), color="red")
        end

        lf = curve_fit(lin, Efit, lines[:, 1], [0.0, 1.0] )
        qf = curve_fit(quad, Efit, lines[:, 1], [0.0, 1.0, 0.0] )
        new_config["label"]["$loc"]["fcal"] = qf.param

        @printf("%4d L: %.3f Q: %.3f\n", loc, maximum(abs.(lf.resid)),
                                maximum(abs.(qf.resid)))
        axl = Axis(fig[3, 1]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    limits=(0, 1500, nothing, nothing), title="Lin")
        axq = Axis(fig[3, 2]; ylabel="ΔE (keV)", xlabel="E (keV)", 
                    limits=(0, 1500, nothing, nothing), title="Quad")
        axs = Axis(fig[3, 3]; ylabel="σ (%)", xlabel="E (keV)", 
                    limits=(0, 1500, nothing, nothing))
        scatter!(axl, lines[:, 1], lf.resid, 
                    color="black", marker=:circle)
        scatter!(axq, lines[:, 1], qf.resid, 
                    color="blue", marker=:utriangle)
        scatter!(axs, lines[:, 1], abs.(sfit),
                color="violet", marker=:star6)
        save(@sprintf("%s_%02d.png", prefix, loc), fig)
    end
    new_config
end

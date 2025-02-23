
"""
    find_shifts(t_loc, config, prefix; dt_min=-1000, dt=1, dt_max=1000)
    
    Based on t_loc table (dt_min:dt:dt_max) calculates t_loc value for
    each channel, returns config Dict with updated values
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
        
        #=
        b = 15
        s = 5
        χ = 0.5

        bs = mean(y[1:b])
        inv = zeros(size(y))
        inv[s+1:end] = y[1:end-s] .- bs
        dy = χ .* (y .- bs) .- inv
        t_lim = argmax(dy)
        t0 = findfirst(x->x<0, dy[t_lim:end]) + t_lim

        a = (dy[t0-1] - dy[t0+1]) / (t[t0-1] - t[t0+1])
        b = dy[t0+1] - t[t0+1] * a
        zc = -b/a
        =#

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

        #=
        new_config["label"]["$loc"]["dt"] = mean(t[xm-7:xm+7], 
                                            weights(y[xm-7:xm+7]))
        =#

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
        #hlines!(ax2, [-thr, thr], color=:gray, linestyle=:dot)
         
        #=
        scatter!(ax, t[xm-7:xm+7], y[xm-7:xm+7], 
                marker=:cross, color=:red)
        vlines!(ax, [new_config["label"]["$loc"]["dt"]], 
                color=:red, linestyle=:dot)
        =#

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
    
    returns config with updated value
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

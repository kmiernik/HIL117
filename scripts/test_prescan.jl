using HIL117
using GLMakie

"""
    find_shifts(t_loc, config, prefix; dt_min=-1000, dt=1, dt_max=1000)
    
    Based on t_loc table (dt_min:dt:dt_max) calculates t_loc value for
    each channel, returns config Dict with updated values
"""
function find_shifts(t_loc, config, prefix; dt_min=-1000, dt=1, dt_max=1000,
                                    doplot=false)
    last_label = size(t_loc)[1]
    new_config = copy(config)
    for loc in 1:last_label
        if sum(t_loc[loc, :]) == 0 || !config["label"]["$loc"]["valid"]
            continue
        end
        y = vec(t_loc[loc, :])
        t = dt_min:dt:dt_max
        if config["label"]["$loc"]["type"] == "Ge"
            y = vec(rebin(y, 4))
            t = dt_min:dt*4:dt_max .+ dt .* 2
        end
        xm = argmax(y)
        new_config["label"]["$loc"]["dt"] = mean(t[xm-7:xm+7], 
                                            weights(y[xm-7:xm+7]))
        fig = Figure(size=(600, 400))
        ax = Axis(fig[1, 1]; title="$loc")
        scatter!(ax, t, y)
        scatter!(ax, t[xm-7:xm+7], y[xm-7:xm+7], 
                marker=:cross, color=:red)
        vlines!(ax, [new_config["label"]["$loc"]["dt"]], 
                color=:red, linestyle=:dot)
        if config["label"]["$loc"]["type"] == "Ge"
            xlims!(ax, t[xm-50], t[xm+50])
        elseif config["label"]["$loc"]["type"] == "BGO"
            xlims!(ax, t[xm-100], t[xm+100])
        elseif config["label"]["$loc"]["type"] == "DIAMANT"
            xlims!(ax, t[xm-100], t[xm+100])
        elseif config["label"]["$loc"]["type"] == "NEDA"
            xlims!(ax, t[xm-20], t[xm+20])
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

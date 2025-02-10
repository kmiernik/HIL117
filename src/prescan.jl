
"""
"""
function prescan(data_dir, configfile::String; args...)

    config = TOML.parsefile(configfile)
    prescan(data_dir, config; args...)
end


"""
    prescan(data_dir, config::Dict; 
                n_files=2, ref_label=34, E_min=20, E_max=8192, Estartmin=100, 
                                    dt_min=-1000.0, dt_max=1000.0, dt=1.0)
    
    Prescan run, take first n_files
    Returns Energy-loc and Time-loc tables
    Time-loc is calculated vs. reference detector (ref_label) 
        and is histogram dt_min:dt:dt_max

    Events in range E_min to E_max are taken to tables
    Reference events have energy larger than Estartmin
"""
function prescan(data_dir, config::Dict; 
                n_files=2, ref_label=34, E_min=20, E_max=8192, Estartmin=100, 
                                    dt_min=-1000.0, dt_max=1000.0, dt=1.0)

    last_label = 0
    for loc in keys(config["label"])
        if parse(Int64, loc) > last_label && config["label"]["$loc"]["valid"] 
            last_label = parse(Int64, loc)
        end
    end

    E_loc = zeros(Int64, last_label, length(1:E_max))
    t_loc = zeros(Int64, last_label, length(dt_min:dt:dt_max))
    if n_files == 0
        @warn "No files to be prescanned"
        return E_loc, t_loc
    end

    files_caen = readdir(data_dir, join=true)
    filter!(x->endswith(x, ".caendat"), files_caen)

    files_dia = readdir(data_dir, join=true)
    filter!(x->split(x, '.')[2] == "dat", files_dia)
    files_dia = [files_dia[1]; sort(files_dia[2:end], by=x->parse(Int64, split(x, '.')[end]))]

    if n_files > length(files_caen)
        @warn "There are no $n_files caendat file(s) to be prescanned, using $(length(files_caen)) file(s) instead"
        n_files = length(files_caen)
    end
    if n_files > length(files_dia)
        @warn "There are no $n_files diamant file(s) to be prescanned, using $(length(files_caen)) file(s) instead"
        n_files = length(files_caen)
    end

    cal_ge = zeros(7, 16)
    for loc in 1:2:32
        cal_ge[:, round(Int64, loc/2, RoundDown)+1] = config["label"]["$loc"]["cal"]
    end

    cfin = open(files_caen[1], "r")
    i_caen = 1
    if parse(Int64, split(files_caen[1], ['/', '_', '.'])[end-1]) == 0
        header = zeros(UInt32, 12)
        read!(cfin, header)
        agava_ts = header[4] % UInt64
        agava_ts = agava_ts << 32 + header[5] % UInt64
    else
        @warn "Could not read agava time stamp"
        return 0
    end
    csize = filesize(files_caen[1])

    i_file = 0
    t_ref = Float64[]
    block_caen = true
    block_dia = true
    n_print = 10_000_000
    time_start = Dates.Time(Dates.now())
    for filename in files_caen
        i_file += 1
        if i_file > n_files
            break
        end

        fin = open(filename, "r")

        if parse(Int64, split(filename, ['/', '_', '.'])[end-1]) == 0
            header = zeros(UInt32, 12)
            read!(fin, header)
        end

        n_hits = 0
        while !eof(fin) 
            hits = RawHit[]
            try
                hits = read_aggregate(fin, config)
            catch err
                break
            end
            for hit in hits
                n_hits += 1
                if n_hits % n_print == 0
                    if isopen(cfin)
                        cpos = position(cfin) / csize
                    else
                        cpos = 1.0
                    end
                    block_caen, block_dia = progress_dots(
                                                  "    \u25E6 ref. time  : ",
                                                    time_start, cpos, 0.0,
                                                    i_caen, block_caen, n_files,
                                                    0, block_dia, 0)
                end
                if hit.board * 16 + hit.ch + 1 != ref_label
                    continue
                end
                t = hit.ts * 4 + hit.tf / 256 
                if hit.E > Estartmin
                    push!(t_ref, t)
                end
            end
        end
        close(fin)
    end
    sort!(t_ref)
    i_ref = 1
    n_ref = length(t_ref)

    progress_dots("    \u25E6 ref. time  : ", time_start, 1.0, 0.0,
                    n_files, block_caen, n_files,
                    0, block_dia, 0)
    println()
    
    i_dia = 1
    dfin = open(files_dia[1], "r")
    dsize = filesize(files_dia[1])

    caen_good = true
    dia_good = true
    t_caen = 0
    t_dia = 0
    
    n_hits = 0

    while caen_good || dia_good

        hits = RawHit[]
        if (caen_good && t_dia >= t_caen) || !(dia_good)
            try
                hits = read_aggregate(cfin, config)
            catch err
                throw(err)
            end
            if size(hits)[1] > 0
                t_caen = hits[end].ts
            end
        elseif (dia_good && t_dia < t_caen) || !(caen_good)
            try
                hits = read_diahit(dfin, agava_ts, 100)
            catch err
                if !eof(dfin)
                    throw(err)
                end
            end
            if size(hits)[1] > 0
                t_dia = hits[end].ts
            end
        end

        for hit in hits
            n_hits += 1

            t = hit.ts * 4 + hit.tf / 256 
            E = hit.E
            loc = hit.board * 16 + hit.ch + 1
            if loc in 1:2:32
                E = (quadquad(hit.E, 
                              cal_ge[:, round(Int64, loc/2, RoundDown)+1])
                     + (rand() - 0.5) * 0.5)
            end
            if E_min < E < E_max
                E_loc[loc, round(Int64, E, RoundUp)] += 1

                if t > t_ref[i_ref]
                    while i_ref < n_ref
                        i_ref += 1
                        if t <= t_ref[i_ref]
                            break
                        end
                    end
                elseif t < t_ref[i_ref]
                    while i_ref > 1
                        i_ref -= 1
                        if t > t_ref[i_ref]
                            i_ref += 1
                            break
                        end
                    end
                end
                if i_ref == 1
                    continue
                end
                dti = ifelse(
                abs(t_ref[i_ref] - t) < abs(t_ref[i_ref-1] - t), 
                    t - t_ref[i_ref], t - t_ref[i_ref-1])
                if dt_min < dti < dt_max
                    it = round(Int, (dti - dt_min) / dt) + 1 
                    t_loc[loc, it] += 1
                end
            end

            if n_hits % n_print == 0
                if isopen(cfin)
                    cpos = position(cfin) / csize
                else
                    cpos = 1.0
                end
                if isopen(dfin)
                    dpos = position(dfin) / dsize
                else
                    dpos = 1.0
                end
                block_caen, block_dia = progress_dots("    \u25E6 prescan    : ",
                                                time_start, cpos, dpos,
                                                i_caen, block_caen, n_files,
                                                i_dia, block_dia, n_files)
            end

        end

        if eof(cfin)
            close(cfin)
            if i_caen < length(files_caen) && i_caen < n_files
                i_caen += 1
                cfin = open(files_caen[i_caen], "r")
                csize = filesize(files_caen[i_caen])
            else
                caen_good = false
            end
        end
        if eof(dfin)
            close(dfin)
            if i_dia < length(files_dia) && i_dia < n_files
                i_dia += 1
                dfin = open(files_dia[i_dia], "r")
                dsize = filesize(files_dia[i_dia])
            else
                dia_good = false
            end
        end
    end
    progress_dots("    \u25E6 prescan    : ", time_start, 1.0, 1.0,
                    i_caen, block_caen, n_files,
                    i_dia, block_dia, n_files)

    E_loc, t_loc
end


"""
    find_shifts(t_loc, config; dt_min=-1000, dt=1, dt_max=1000)
    
    Based on t_loc table (dt_min:dt:dt_max) calculates t_loc value for
    each channel, returns config Dict with updated values
"""
function find_shifts(t_loc, config; dt_min=-1000, dt=1, dt_max=1000,
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
        if doplot
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
            save(@sprintf("t_loc_%03d.png", loc), fig)
        end
    end
    new_config
end



"""
    find_period(t_loc, config; period_loc=35, dt_min=-1000, dt=1, 
                                        dt_max=1000)

    Find beam period based on NEDA detector timing

    512 ns samples in a "safe" range (wihtout main peak)
    detector next to ref (default 34, first NEDA) 
    shouldn't have scattering "echos"
    finding minima in the safe range will give off-beam period
    
    returns config with updated value
"""
function find_period(t_loc, config; period_loc=35, dt_min=-1000, dt=1, 
                                    t_safe=200,
                                    dt_max=1000, doplot=false)
    
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
    if doplot
        fig = Figure(size=(600, 400))
        ax = Axis(fig[1, 1]; xlabel="t (ns)", 
                title=@sprintf("%.3f ns", new_config["spectra"]["beam_period"]))
        scatter!(ax, tm, ym)
        vlines!(ax, localmins, color=:red, linestyle=:dash)
        save("period.png", fig)
    end

    new_config
end


function fine_cal(E_loc, configfile::String; args...)
    config = TOML.parsefile(configfile)
    fine_cal(E_loc, config; args...)
end


"""
    fine_cal(E_loc, configfile; doplot=false)

    Take E_loc table and perform fine calibration for the run
    Returns new config with fcal fields for Ge
"""
function fine_cal(E_loc, config::Dict; doplot=false)
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

            if doplot
                ix = round(Int, (ip-1)/4, RoundDown)+1
                iy = (ip-1)%4+1
                ax = Axis(fig[ix, iy]; title="$(lines[ip, 1])")
                stairs!(ax, E0-dE-10:E0+dE+10, 
                        vec(E_loc[loc, E0-dE-10:E0+dE+10]))
                lines!(ax, E0-dE:0.1:E0+dE, gausslin(E0-dE:0.1:E0+dE, 
                                                pf.param), color="red")
            end
        end

        lf = curve_fit(lin, Efit, lines[:, 1], [0.0, 1.0] )
        qf = curve_fit(quad, Efit, lines[:, 1], [0.0, 1.0, 0.0] )
        new_config["label"]["$loc"]["fcal"] = qf.param
        if doplot
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
            save(@sprintf("fcal_%02d.png", loc), fig)
        end
    end
    new_config
end

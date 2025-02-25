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
    Returns Energy-loc (raw ch numbers) and Time-loc tables
    Time-loc is calculated vs. reference detector (ref_label) 
        and is histogram dt_min:dt:dt_max

    Events in range E_min to E_max are taken to tables
    Reference events have energy larger than Estartmin
"""
function prescan(data_dir, config::Dict; 
                n_files=2, ref_label=34, E_min=20, E_max=32768, Estartmin=100, 
                                    dt_min=-1000.0, dt_max=1000.0, dt=1.0)

    last_label = 0
    for loc in keys(config["label"])
        if parse(Int64, loc) > last_label && config["label"]["$loc"]["valid"] 
            last_label = parse(Int64, loc)
        end
    end
    
    pidpars = PidPars(config)

    E_loc = zeros(Int64, last_label, length(1:E_max))
    t_loc = zeros(Int64, last_label, length(dt_min:dt:dt_max))
    if n_files == 0
        #@warn "No files to be prescanned"
        return E_loc, t_loc
    end

    files_caen = readdir(data_dir, join=true)
    filter!(x->endswith(x, ".caendat"), files_caen)

    if length(files_caen) == 0
        #@warn "No caendat files found"
        return E_loc, t_loc
    end

    files_dia = readdir(data_dir, join=true)
    filter!(x->split(basename(x), '.')[2] == "dat", files_dia)
    if length(files_dia) > 0
        files_dia = [files_dia[1]; sort(files_dia[2:end], by=x->parse(Int64, split(basename(x), '.')[end]))]
    end

    if n_files > length(files_caen)
        #@warn "There are no $n_files caendat file(s) to be prescanned, using $(length(files_caen)) file(s) instead"
        n_files = length(files_caen)
    end

    if n_files > length(files_caen) && n_files > length(files_dia)
        n_new = maximum(length(files_caen), length(files_dia))
        #@warn "There are no both $n_files caendat and diamant file(s) to be prescanned, using $(n_new) file(s) instead"
        n_files = n_new
    elseif n_files > length(files_dia)
        #@warn "There are no $n_files diamant file(s) to be prescanned"
    elseif n_files > length(files_caen)
        #@warn "There are no $n_files caendat file(s) to be prescanned"
    end

    cfin = open(files_caen[1], "r")
    i_caen = 1
    if parse(Int64, split(files_caen[1], ['/', '_', '.'])[end-1]) == 0
        header = zeros(UInt32, 12)
        read!(cfin, header)
        agava_ts = header[4] % UInt64
        agava_ts = agava_ts << 32 + header[5] % UInt64
    else
        #@warn "Could not read agava time stamp"
        return 1
    end
    csize = filesize(files_caen[1])

    totsize = 0
    donesize = 0
    for filename in files_caen[1:n_files]
        totsize += filesize(filename)
    end
    prog = Progress(totsize; dt=1.0, desc="\tRef. timing ",
                     barlen=30, color=:green)

    i_file = 0
    t_ref = Float64[]
    n_print = 10_000_000
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
                msg = @sprintf("File %s is possibly corrupted at position %d. Skipping rest of the file.", filename, position(fin)) 
                @warn msg
                seekend(fin)
            end
            for hit in hits
                n_hits += 1
                if n_hits % n_print == 0
                    current_pos = 0
                    if isopen(cfin)
                        current_pos = position(fin)
                    end
                    update!(prog, donesize + current_pos)
                end
                if hit.board * 16 + hit.ch + 1 != ref_label
                    continue
                end
                t = hit.ts * 4 + hit.tf / 256 
                if hit.E > Estartmin
                    pid = (Float64(hit.qshort) + randn()) / Float64(hit.E)
                    if (pidpars.g_low <= pid < pidpars.g_high)
                        push!(t_ref, t)
                    end
                end
            end
        end
        donesize += filesize(filename)
        close(fin)
    end
    update!(prog, totsize)
    sort!(t_ref)
    i_ref = 1
    n_ref = length(t_ref)
    
    caen_good = true
    dia_good = true
    i_dia = 1
    if length(files_dia) > 0
        dfin = open(files_dia[1], "r")
        dsize = filesize(files_dia[1])
    else
        #@warn "No diamant files found"
        dfin = IOBuffer()
        close(dfin)
        dsize = 0
        dia_good = false
        i_dia = n_files 
    end

    totsize = 0
    donesize = 0
    for filename in files_caen[1:n_files]
        totsize += filesize(filename)
    end
    n_dia = ifelse(length(files_dia) >= n_files, n_files, length(files_dia))
    for filename in files_dia[1:n_dia]
        totsize += filesize(filename)
    end
    prog = Progress(totsize; dt=1.0, desc="\tPrescan     ", 
                    barlen=30, color=:green)
    
    t_caen = 0
    t_dia = 0
    n_hits = 0

    while caen_good || dia_good

        hits = RawHit[]
        if (caen_good && t_dia >= t_caen) || !(dia_good)
            try
                hits = read_aggregate(cfin, config)
            catch err
                msg = @sprintf("File %s is possibly corrupted at position %d. Skipping rest of the file.", files_caen[i_caen], position(cfin)) 
                @warn msg
                seekend(cfin)
            end
            if size(hits)[1] > 0
                t_caen = hits[end].ts
            end
        elseif (dia_good && t_dia < t_caen) || !(caen_good)
            try
                hits = read_diahit(dfin, agava_ts, 100)
            catch err
                if !eof(dfin)
                    msg = @sprintf("File %s is possibly corrupted at position %d. Skipping rest of the file.", files_dia[i_dia], position(dfin)) 
                    @warn msg
                    seekend(dfin)
                end
            end
            if size(hits)[1] > 0
                t_dia = hits[end].ts
            end
        end

        for hit in hits
            n_hits += 1

            t = hit.ts * 4 + hit.tf / 256 
            loc = hit.board * 16 + hit.ch + 1
            if E_min < hit.E < E_max
                E_loc[loc, hit.E] += 1

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
                    it = round(Int, (dti - dt_min) / dt, RoundUp)
                    t_loc[loc, it] += 1
                end
            end

            if n_hits % n_print == 0
                current_pos = 0
                if isopen(cfin)
                    current_pos += position(cfin)
                end
                if isopen(dfin)
                    current_pos += position(dfin)
                end
                update!(prog, donesize + current_pos)
            end

        end

        if eof(cfin)
            close(cfin)
            if i_caen < length(files_caen) && i_caen < n_files
                donesize += filesize(files_caen[i_caen])
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
                donesize += filesize(files_dia[i_dia])
                i_dia += 1
                dfin = open(files_dia[i_dia], "r")
                dsize = filesize(files_dia[i_dia])
            else
                dia_good = false
            end
        end
    end
    update!(prog, totsize)

    E_loc, t_loc
end


"""
    find_shifts(t_loc, config; dt_min=-1000, dt=1, dt_max=1000)
    
    Based on t_loc table (dt_min:dt:dt_max) calculates t_loc value for
    each channel, returns config Dict with updated values
"""
function find_shifts(t_loc, config; dt_min=-1000, dt=1, dt_max=1000)
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
    end
    new_config
end


"""
    find_period(t_loc, config; period_loc=35, dt_min=-1000, dt=1, 
                                    t_safe=200, dt_max=1000)

    Find beam period based on NEDA detector timing

    512 ns samples in a "safe" range (wihtout main peak)
    detector next to ref (default 34, first NEDA) 
    shouldn't have scattering "echos"
    finding minima in the safe range will give off-beam period
    
    returns config with updated value
"""
function find_period(t_loc, ; period_loc=35, dt_min=-1000, dt=1, 
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
    period = mean(localmins[2:end] - localmins[1:end-1])
    return period
end


function find_cal(E_loc, configfile::String; args...)
    config = TOML.parsefile(configfile)
    find_cal(E_loc, config; args...)
end


"""
    find_cal(E_loc, configfile; doplot=false)

    Take E_loc table and perform calibration for the run
    Returns new config with cal fields for Ge (linear)
"""
function find_cal(E_loc, config::Dict)
    new_config = copy(config)
    
    E1 = config["calibration"]["E1"]
    E2 = config["calibration"]["E2"]
    last_label = size(E_loc)[1]
    x = collect(1:size(E_loc)[2]) .- 0.5
    for loc in 1:last_label
        if (!haskey(config["label"], "$loc")
            || !config["label"]["$loc"]["valid"] 
            || config["label"]["$loc"]["type"] != "Ge")
            continue
        end
        
        y = E_loc[loc, :]

        x1 = argmax(y[100:8000]) + 100
        ab = E1 / x1
        ch2 = round(Int, E2 / ab)

        x2 = argmax(y[ch2-100:ch2+100]) + ch2 - 100

        pf1 = curve_fit(gausslin, x[x1-50:x1+50], y[x1-50:x1+50], 
                       [5.0*y[x1], x1, 10.0, y[x1-50], 0.0])
        ch1 = pf1.param[2]

        guess = [5.0*y[x2], x2, 10.0, y[x2-50], 0.0]
        pf2 = curve_fit(gausslin, x[x2-50:x2+50], 
                           y[x2-50:x2+50], guess)
        ch2 = pf2.param[2]
        a2 = (E2 - E1) / (ch2 - ch1)
        a1 = E2 - a2 * ch2

        new_config["label"]["$loc"]["cal"] = [a1, a2]
    end
    new_config
end


function prescan_all(dirname, configfile; n_prescan=2, skipruns=Int64[])
    dirs = readdir(dirname, join=true)
    for dir in dirs
        if isfile(dir)
            continue
        end
        if !startswith(basename(dir), "run_0")
            continue
        end
        run_number = parse(Int64, split(basename(dir), ['_', '.'])[2])
        println("\u25CD Run $run_number ")
        if run_number in skipruns
            println("   skipping")
            continue
        end
        config = TOML.parsefile(configfile)
        E_loc, t_loc = prescan(dir, config; n_files=n_prescan)
        print("    \u25E6 Calculating shifts, calibration, and period (")
        config = find_shifts(t_loc, config)
        config = find_cal(E_loc, config)
        calc_period = find_period(t_loc)
        @printf("%.3f ns)\n", calc_period)
        if abs(config["event"]["beam_period"] - calc_period) > 1.0
            @warn "Calculated period is different than one in config!"
        end
        nicetoml(config, @sprintf("config/config_%03d.toml", run_number))
    end
end

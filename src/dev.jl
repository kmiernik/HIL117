#=
Notes
-----

1. test_co60() - data read test 
2. init_cal() - rough calibration with Co-60
2a. scan_raw_ge() - scan whole runs, only 16x16384 table for ge detectors
3. cal_energy() - calibration with Eu-152 source
4. test_walk() - walk correction per channel - not needed, no dependency
                of time shift on energy; different methods tested, e.g.
                Co60 coincidence, (one 1332, one 1173 or lower), Eu-152, etc.
5. find_shift_table() - calculates shift table for each channel relative to 
                    the selected NEDA detector
6. find_shift() - time shifts for each channel based on shift table
7. find_period() - beam period based on shift table

=#

function test_co60()
    data = read_raw_data("data/needleRU_i1329_0004_0000.caendat"; n_agg=100000)

    spectra = Dict{Int, Dict{Int, Vector{Int64}}}()
    for board in 0:5
        chans = Dict{Int, Vector{Int64}}()
        for ch in 0:15
            chans[ch] = zeros(Int64, 16384)
        end
        spectra[board] = chans
    end

    for hit in data
        if 1 < hit.E <= 16384
            spectra[hit.board][hit.ch][hit.E] += 1
        end
    end
    spectra
end


function scan_raw_ge(dirname, configfile; file_numbers=[], chmax=32768)
    config = TOML.parsefile(configfile)
    files = readdir(dirname, join=true)
    filter!(x->endswith(x, ".caendat"), files)

    if length(file_numbers) > 0
        filter!(x->parse(Int64, split(x, ['/', '_', '.'])[end-1]) in file_numbers, files)
    end

    spectra = zeros(Int64, 16, chmax)
    for filename in files
        fin = open(filename, "r")
        name = split(filename, ['/'])[end]
        println("\u25CB Processing file $name")

        if parse(Int64, split(filename, ['/', '_', '.'])[end-1]) == 0
            header = zeros(UInt32, 12)
            read!(fin, header)
        end

        i_hit = 0
        i_symbol = 0
        while !eof(fin) 
            hits = RawHit[]
            try
                hits = read_aggregate(fin, config)
            catch err
                println(err)
                break
            end
            for hit in hits
                i_hit += 1
                if i_hit % 1000000 == 0
                    i_symbol += 1
                    if i_symbol == 1
                        @printf("\r \u25D0 %6.2e ", i_hit)
                    elseif i_symbol == 2
                        @printf("\r \u25D3 %6.2e ", i_hit)
                    elseif i_symbol == 3
                        @printf("\r \u25D1 %6.2e ", i_hit)
                    elseif i_symbol == 4
                        @printf("\r \u25D2 %6.2e ", i_hit)
                        i_symbol = 0
                    end
                end
                if hit.board > 1 || hit.ch % 2 != 0
                    continue
                end
                ige = hit.board * 8 + round(Int, hit.ch / 2, RoundDown) + 1
                if 1 <= hit.E <= chmax
                    spectra[ige, hit.E] += 1
                end
            end
        end
        @printf("\r \u25CF %6.2e \n", i_hit)
    end
    spectra
end


"""

    Rough Ge calibration based on Co-60

    Returns vector of two elements vector
    Ge B0C0, B0C2, ..., B0C14, B1C0, ..., B1C12, B1C14
        1     2    ....  8      9          15    16
    ige = board * 8 + round(Int, ch / 2, RoundDown) + 1
    
"""
function init_cal(filename)
    data = read_raw_data(filename; n_agg=100000)

    spectra = Dict{Int, Dict{Int, Vector{Int64}}}()
    for board in 0:5
        chans = Dict{Int, Vector{Int64}}()
        for ch in 0:15
            chans[ch] = zeros(Int64, 16384)
        end
        spectra[board] = chans
    end
    for hit in data
        if 1 < hit.E <= 16384
            spectra[hit.board][hit.ch][hit.E] += 1
        end
    end

    cal_ge = Vector{Vector{Float64}}()
    for b in 0:1
        for c in 0:2:15
            push!(cal_ge, [0.0, 0.0])
        end
    end
    thr = 20
    for b in 0:1
        for c in 0:2:15
            fits = Float64[]
            if sum(spectra[b][c]) == 0
                continue
            end
            ds = spectra[b][c][2:end] .- spectra[b][c][1:end-1]
            ds = MTools.rebin(ds, 4)
            up = 0
            down = 0
            for i in 1:size(ds)[1]
                if ds[i] > thr 
                    up = i
                elseif ds[i] < -thr
                    down = i
                    if down - up < 10
                        x1 = (up - 10) * 4
                        x2 = (down + 10) *4
                        pf = curve_fit(gausslin, x1:x2, spectra[b][c][x1:x2],
                                       [sum(spectra[b][c][x1:x2]), (x1+x2)/2,
                                            10, spectra[b][c][x1], 0.0],
                                        lower=[0, x1, 1, -Inf, -Inf],
                                        upper=[Inf, x2, 20, Inf, Inf])
                        push!(fits, pf.param[2])
                    end
                    up = 0
                    down = 0
                end
            end
            lf = curve_fit(lin, fits, [1173.2, 1332.5], [0.0, 1.0])
            ige = b * 8 + round(Int, c / 2, RoundDown) + 1
            cal_ge[ige] = lf.param
        end
    end
    return cal_ge
end


function linquad(x::Vector, p)
    y = zeros(size(x))
    for (i, xi) in enumerate(x)
        y[i] = linquad(xi, p)
    end
    y
end


function linquad(x::Real, p)
    if x < p[6]
        return p[1] + p[2] * x
    else
        return p[3] + p[4] * x + p[5] * x^2
    end
end

function linlin(x::Vector, p)
    y = zeros(size(x))
    for (i, xi) in enumerate(x)
        y[i] = linlin(xi, p)
    end
    y
end


function linlin(x::Real, p)
    if x < p[5]
        return p[1] + p[2] * x
    else
        return p[3] + p[4] * x
    end
end

"""
    Calibrate using Eu-152 data (energy calibration)

    Data are loaded from output of scan_raw_ge (16x16384 array)
"""
function cal_energy(filename, configfile)

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
        
    end
    new_config
end


function test_walk(data_dir, configfile::String; args...)
    config = TOML.parsefile(configfile)
    test_walk(data_dir, config; args...)
end

"""
    ΔT vs. E of a given detector
    for DIAMANT ΔT vs. PID

    Testing time widths and walk
"""
function test_walk(data_dir, config::Dict; 
                n_files=2, ref_label=34, dE=8, Emax=16384, Estartmin=100, 
                                    dt_min=-99.0, dt_max=900.0, dt=1.0)

    last_label = 0
    for loc in keys(config["label"])
        if parse(Int64, loc) > last_label && config["label"]["$loc"]["valid"] 
            last_label = parse(Int64, loc)
        end
    end
    type_table, valid_table, cal_table, fcal_table, shift_table, pid_table = config_to_tables(config, last_label)

    walk = zeros(Int32, last_label, round(Int64, Emax / dE, RoundUp), 
                 length(dt_min:dt:dt_max))
    #walk = zeros(Int32, 64, round(Int64, Emax / dE, RoundUp), 
    #             length(dt_min:dt:dt_max), length(0:0.01:1.27))
    if n_files == 0
        @warn "No files to be prescanned"
        return walk
    end

    files_caen = readdir(data_dir, join=true)
    filter!(x->endswith(x, ".caendat"), files_caen)

    if length(files_caen) < 1
        @warn "No caendat files found"
        return 0
    end

    files_dia = readdir(data_dir, join=true)
    filter!(x->split(x, '.')[2] == "dat", files_dia)
    if length(files_dia) > 0
        files_dia = [files_dia[1]; sort(files_dia[2:end], by=x->parse(Int64, split(x, '.')[end]))]
    end

    if n_files > length(files_caen) && n_files > length(files_dia)
        n_new = maximum([length(files_caen), length(files_dia)])
        @warn "There are no both $n_files caendat and diamant file(s) to be prescanned, using $(n_new) file(s) instead"
        n_files = n_new
    elseif n_files > length(files_dia)
        @warn "There are no $n_files diamant file(s) to be prescanned"
    elseif n_files > length(files_caen)
        @warn "There are no $n_files caendat file(s) to be prescanned"
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
                    if isopen(fin)
                        cpos = position(fin) / csize
                    else
                        cpos = 1.0
                    end
                    block_caen, block_dia = progress_dots(
                                                  "    \u25E6 ref. time  : ",
                                                    time_start, cpos, 0.0,
                                                    i_file, block_caen, n_files,
                                                    0, block_dia, 0)
                end
                if hit.board * 16 + hit.ch + 1 != ref_label
                    continue
                end
                t = hit.ts * 4 + hit.tf / 256 
                if hit.E > Estartmin
                    pid = (Float64(hit.qshort) + randn()) / Float64(hit.E)
                    if hit.E > pid_table[1] && (0.75 <= pid < 0.95)
                        #Good gamma
                        push!(t_ref, t)
                    end
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
    
    caen_good = true
    dia_good = true
    t_caen = 0
    t_dia = 0
    i_dia = 1
    if length(files_dia) > 0
        dfin = open(files_dia[1], "r")
        dsize = filesize(files_dia[1])
    else
        dfin = IOBuffer()
        close(dfin)
        dsize = 0
        dia_good = false
        i_dia = n_files 
    end

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
            E = Float64(hit.E)
            loc = hit.board * 16 + hit.ch + 1
            pid = 0.0
            if type_table[loc] == 1
                ige = round(Int64, loc / 2, RoundDown) + 1
                E = quad(quadquad(
                    hit.E, cal_table[:, loc]) + (rand()-0.5), 
                            fcal_table[:, ige]) * dE
            elseif type_table[loc] == 3
                pid = (Float64(hit.qshort) + randn()) / E
                E = pid * 10000
            elseif type_table[loc] == 4
                T = Float64(hit.qshort)
                pid = (T / E 
                    - cal_table[1, loc] / E^2 
                    - cal_table[2, loc] + 0.5)
                E = pid * 10000
            else
                continue
            end
            if 1 < E < Emax
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
                # Only previous t_ref
                #dti = t - t_ref[i_ref-1]
                if dt_min < dti < dt_max 
                    it = round(Int, (dti - dt_min) / dt, RoundUp)
                    iE = round(Int, E / dE, RoundUp)
                    walk[loc, iE, it] += 1
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
                block_caen, block_dia = progress_dots(
                                                "    \u25E6 prescan    : ",
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
    walk
end


"""
    find_shift_table(config_file, dir_name; ref_label=34,
                                    Emin = 20, Emax = 8192, Estartmin = 100, 
                                    dt_min = -1000.0, dt_max = 1000.0, dt = 1.0)
    
    Calculate shift table (for each channel) relative to the ref_label
    detector (best first NEDA detector).

    Emin - minimal ch number of detector to be taken
    Emax - maximum ch number of detector to be taken
    Estartmin - minimal ch number in the ref detector to be taken 
                as a start signal
    dt_min / dt / dt_max - time table dt_min:dt_dt_max

    returns number of labels vs. time table 2D array

"""
function find_shift_table(config_file, dir_name; ref_label=34,
                                    Emin=20, Emax=8192, Estartmin=100, 
                                    dt_min=-1000.0, dt_max=1000.0, dt=1.0)

    config = TOML.parsefile(config_file)
    files = readdir(dir_name, join=true)
    filter!(x->endswith(x, ".caendat"), files)

    last_label = 0
    for loc in keys(config["label"])
        if parse(Int64, loc) > last_label && config["label"]["$loc"]["type"] == "NEDA"
            last_label = parse(Int64, loc)
        end
    end

    cal_table = zeros(2, last_label)
    for label in 1:last_label
        cal_table[:, label] = config["label"]["$label"]["rcal"]
    end

    shift = zeros(Int64, last_label, length(dt_min:dt:dt_max))

    i_file = 0
    for filename in files
        time_start = Dates.Time(Dates.now())
        i_file += 1
        fin = open(filename, "r")
        name = split(filename, ['/'])[end]
        println("\nProcessing file $name  $i_file / $(length(files))")

        if i_file == 1
            header = zeros(UInt32, 12)
            read!(fin, header)
        end

        t_start = Float64[]
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
                if hit.board * 16 + hit.ch + 1 != ref_label
                    continue
                end
                t = hit.ts * 4 + hit.tf / 256  # tf * 4 / 1024
                E = lin(hit.E, cal_table[:, ref_label])
                if hit.E > Estartmin
                    push!(t_start, t)
                end
            end
        end
        seek(fin, 0)
        if i_file == 1
            read!(fin, header)
        end
        sort!(t_start)
        n_start = size(t_start)[1]
        @show n_start
        @show n_hits

        if size(t_start)[1] == 0
            continue
        end

        i_start = 1
        i_hit = 0
        while !eof(fin) 
            hits = RawHit[]
            try
                hits = read_aggregate(fin, config)
            catch err
                break
            end
            for hit in hits
                i_hit += 1
                loc = hit.board * 16 + hit.ch + 1
                if loc > last_label
                    continue
                end
                t = hit.ts * 4 + hit.tf / 256 
                E = lin(hit.E, cal_table[:, ref_label])
                if Emin < hit.E < Emax
                    if t > t_start[i_start]
                        while i_start < n_start
                            i_start += 1
                            if t <= t_start[i_start]
                                break
                            end
                        end
                    elseif t < t_start[i_start]
                        while i_start > 1
                            i_start -= 1
                            if t > t_start[i_start]
                                i_start += 1
                                break
                            end
                        end
                    end
                    if i_start == 1
                        continue
                    end
                    dti = ifelse(
                    abs(t_start[i_start] - t) < abs(t_start[i_start-1] - t), 
                        t - t_start[i_start], t - t_start[i_start-1])
                    if dt_min < dti < dt_max
                        it = round(Int, (dti - dt_min) / dt) + 1 
                        shift[loc, it] += 1
                    end
                end
                if i_hit % 100000 == 0
                    dtime = (Dates.Time(Dates.now()) - time_start)
                    @printf("\r %8.2f%% %8.1f s ", i_hit / n_hits * 100,
                            dtime.value * 1e-9)
                end
            end
        end
    end
    println()
    shift
end


function test_neda(dir_name, configfile::String)
    config = TOML.parsefile(configfile)
    test_neda(dir_name, config)
end

"""

    Test NEDA PID n/g based on short/long ratio
"""
function test_neda(dir_name, config::Dict)
    files = readdir(dir_name, join=true)
    filter!(x->endswith(x, ".caendat"), files)

    last_label = 0
    neda_locs = Int64[]
    for loc in keys(config["label"])
        if config["label"]["$loc"]["type"] == "NEDA"
            push!(neda_locs, parse(Int64, loc))
        end

        if parse(Int64, loc) > last_label && config["label"]["$loc"]["type"] == "NEDA"
            last_label = parse(Int64, loc)
        end
    end

    Lmax = 2048
    Smax = 128
    neda = zeros(Int64, length(neda_locs), Lmax, Smax)

    totsize = 0
    for filename in files
        totsize += filesize(filename)
    end

    prog = Progress(totsize; dt=1.0, desc="Calculating NEDA PID spectra",
                    barglyphs=BarGlyphs("[=> ]"), barlen=25, color=:black)
    donesize = 0

    i_file = 0
    for filename in files
        time_start = Dates.Time(Dates.now())
        i_file += 1
        fin = open(filename, "r")
        name = split(filename, ['/'])[end]

        if parse(Int64, split(filename, ['/', '_', '.'])[end-1]) == 0
            header = zeros(UInt32, 12)
            read!(fin, header)
            agava_ts = header[4] % UInt64
            agava_ts = agava_ts << 32 + header[5] % UInt64
        end

        i_hit = 0
        while !eof(fin) 
            hits = RawHit[]
            try
                hits = read_aggregate(fin, config)
            catch err
                break
            end
            for hit in hits
                i_hit += 1
                loc = hit.board * 16 + hit.ch + 1
                if loc > last_label
                    continue
                end
                if i_hit % 10_000_000 == 0
                    update!(prog, donesize + position(fin))
                end
                if loc in neda_locs
                    ineda = loc - 32
                    L = hit.E
                    S = (hit.qshort + randn()) / L * 100
                    if 1 <= L < Lmax && 1 <= S < Smax
                        neda[ineda, round(Int64, L, RoundUp), 
                                    round(Int64, S, RoundUp)] += 1
                    end
                end
            end
        end
        donesize += filesize(filename)
    end
    println()
    neda
end


function dia_banana_fit(dia_dir_name, configfile::String;
                        args...)
    config = TOML.parsefile(configfile)
    dia_banana_fit(dia_dir_name, config; args...)
end

"""
"""
function dia_banana_fit(dia_dir_name, config::Dict; 
                        n_files=3, c_lim=3000, E_lim=3000, arr_mode=false,
                        pid_mode=false)

    time_start = Dates.Time(Dates.now())

    files_dia = readdir(dia_dir_name, join=true)
    filter!(x->split(x, '.')[2] == "dat", files_dia)
    files_dia = [files_dia[1]; sort(files_dia[2:end], by=x->parse(Int64, split(x, '.')[end]))]

    loc_dia_start = -1
    loc_dia_end = -1
    for loc in keys(config["label"])
        l = parse(Int64, loc)
        if config["label"]["$loc"]["type"] == "DIAMANT"
            if l > loc_dia_end || loc_dia_end < 0
                loc_dia_end = l
            end
            if l < loc_dia_start || loc_dia_start < 0
                loc_dia_start = l
            end
        end
    end

    cal_dia = zeros(2, loc_dia_end-loc_dia_start+1)
    if pid_mode
        for loc in keys(config["label"])
            l = parse(Int64, loc)
            if config["label"]["$loc"]["type"] == "DIAMANT"
                cal_dia[:, l-loc_dia_start+1] = config["label"]["$loc"]["cal"]
            end
        end
    end

    i_dia = 1
    dfin = open(files_dia[1], "r")

    dia_good = true
    
    dia_id = zeros(Int32, 64, 128, 128)

    block_dia = true
    i_hit = 0
    while dia_good

        try
            hit = read_diahit(dfin, 0)

            loc = hit.board * 16 + hit.ch + 1
            i_hit += 1
            if loc_dia_start <= loc <= loc_dia_end
                iE = round(Int64, hit.E / 100, RoundDown)
                iT = 0
                if pid_mode
                    if hit.E > 100
                        E = Float64(hit.E)
                        T = Float64(hit.qshort)
                        pid = (T / E 
                            - cal_dia[1, loc-loc_dia_start+1] / E^2 
                            - cal_dia[2, loc-loc_dia_start+1] + 0.5)
                        iT = round(Int64, pid * 100)
                    end
                else
                    iT = round(Int64, hit.qshort / hit.E * 100)
                end
                if 1 <= iE <= 128 && 1 <= iT <= 128
                    dia_id[loc-loc_dia_start+1, iE, iT] += 1
                end
            end
            if i_hit % 10_000_000 == 0
                print("\r")
                for i in 1:i_dia-1
                    print("\u25CD") 
                end
                if block_dia
                    print("\u25CE")
                    block_dia = false
                else
                    print("\u25E6")
                    block_dia = true
                end
                for i in 1:n_files-i_dia
                    print("\u25CC")
                end
                dtime = (Dates.Time(Dates.now()) - time_start)
                @printf(" %8.2f s ", dtime.value * 1e-9)
            end
        catch err
            if eof(dfin)
                close(dfin)
                if i_dia < n_files
                    i_dia += 1
                    dfin = open(files_dia[i_dia], "r")
                else
                    dia_good = false
                end
            else
                throw(err)
                break
            end
        end
    end

    if arr_mode
        return dia_id
    end

    new_config = copy(config)

    for loc in loc_dia_start:loc_dia_end
        di = loc - loc_dia_start + 1
        bp = fit_ban(dia_id[di, :, :], c_lim, E_lim, loc)
        new_config["label"]["$loc"]["cal"] = bp
    end
    new_config

end


function fit_ban(dia_spec, c_lim, E_lim, loc)        
    Ai = Float64[]
    Pi = Float64[]
    Ei = Float64[]
    fi = 1
    t = 0.00:0.01:1.27
    b(x, p) = @. p[1] / x^2 + p[2]
    for ie in 1:size(dia_spec)[1]
        if maximum(dia_spec[ie, :]) > c_lim && (100 * (ie - 1) < E_lim)
            tg = t[argmax(dia_spec[ie, :])]
            guess = [sum(dia_spec[ie, :]) * 0.015, tg, 0.04, 
                         sum(dia_spec[ie, :]) * 0.003, tg+0.12, 0.04]
            gf = curve_fit(ngauss, t, dia_spec[ie, :], guess)
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
    bf = curve_fit(b, Ei, Ai, [30000.0, 0.6])
    bf.param
end


"""
    Δt between "trigger" signals

"""
function dt_trig(data_dir, configfile; ref_label=93, n_files=1, t_max=1_000_000)

    dt = zeros(Int64, length(0:1:t_max))

    files_caen = readdir(data_dir, join=true)
    filter!(x->endswith(x, ".caendat"), files_caen)
    run_number = parse(Int64, split(files_caen[1], ['/', '_', '.'])[end-2])
    config = TOML.parsefile(configfile)

    i_file = 0
    block_caen = true
    block_dia = true
    n_print = 10_000_000
    time_start = Dates.Time(Dates.now())
    t_last = 0.0
    for filename in files_caen
        i_file += 1
        if i_file > n_files
            break
        end

        fin = open(filename, "r")
        cpos = 0.0
        csize = filesize(filename)

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
                    if isopen(fin)
                        cpos = position(fin) / csize
                    else
                        cpos = 1.0
                    end
                    block_caen, block_dia = progress_dots(
                                                  "    \u25E6 ref. time  : ",
                                                    time_start, cpos, 0.0,
                                                    i_file, block_caen, n_files,
                                                    0, block_dia, 0)
                end
                if hit.board * 16 + hit.ch + 1 != ref_label
                    continue
                end
                t = hit.ts * 4 + hit.tf / 256 
                it = round(Int64, t - t_last, RoundDown) + 1
                t_last = t
                if 1 <= it <= t_max
                    dt[it] += 1
                end
            end
        end
        close(fin)
    end
    dt
end


function neda_lims(neda, config; binsize=8, thr_g=80, thr_n=16)
    new_config = copy(config)
    for loc in 34:84
        E = collect(1:2048) .- 0.5
        pid = collect(0.0:0.01:1.27) .+ 0.005
    
        yn = rebin(vec(sum(neda[loc-32, :, 50:69], dims=2)), binsize)[:, 1]
        yg = rebin(vec(sum(neda[loc-32, :, 71:95], dims=2)), binsize)[:, 1]
        dyn = yn[2:end] .- yn[1:end-1]
        dyg = yg[2:end] .- yg[1:end-1]

        zcn = -1
        zcg = -1
        if maximum(yg[1:round(Int, 80/binsize)]) > thr_g
            # minimum mode
            for i in round(Int, 80/binsize):round(Int, 2000/binsize)
                if yg[i] > thr_g && dyg[i] < 0 && dyg[i+1] > 0
                    a = (dyg[i+1] - dyg[i]) / binsize
                    b = dyg[i] - ((i-1)*binsize + (binsize-1)/2) * a
                    zcg = -b/a
                    break
                end
            end
        else
            # thr_geshold mode
            for i in 1:round(Int, 2000/binsize)
                if zcg < 0 && yg[i] > thr_g && yg[i+1] > thr_g
                    zcg = (i-1)*binsize + (binsize-1)/2
                    break
                end
            end
        end
        if maximum(yn[1:round(Int, 80/binsize)]) > thr_n
            # minimum mode
            for i in round(Int, 80/binsize):round(Int, 2000/binsize)
                if yn[i] > thr_n && dyn[i] < 0 && dyn[i+1] > 0
                    a = (dyn[i+1] - dyn[i]) / binsize
                    b = dyn[i] - ((i-1)*binsize + (binsize-1)/2) * a
                    zcn = -b/a
                    break
                end
            end
        else
            # thr_neshold mode
            for i in 1:round(Int, 2000/binsize)
                if yn[i] > thr_n && yn[i+1] > thr_n
                    zcn = (i-1)*binsize + (binsize-1)/2
                    break
                end
            end
        end

        # Manual intervention!
        if loc == 82
            zcg = 80.0
        end
        new_config["label"]["$loc"]["cal"] = [zcg, zcn]
    end  
    new_config
end



"""
    Δt between group of detectors (Ge, Neda, Diamant)

"""
function dt_groups(data_dir, configfile; n_files=2, dt_max=1023.0, dt=1.0, Emin=100, chunk_size=100_000)
    config = TOML.parsefile(configfile)

    last_label = 0
    for loc in keys(config["label"])
        if parse(Int64, loc) > last_label && config["label"]["$loc"]["valid"] 
            last_label = parse(Int64, loc)
        end
    end

    type_table, valid_table, cal_table, shift_table = config_to_tables(config, last_label)
    
    dtgroup = zeros(Int64, 4, length(0:dt:dt_max))
    if n_files == 0
        @warn "No files to be prescanned"
        return dtgroup
    end

    files_caen = readdir(data_dir, join=true)
    filter!(x->endswith(x, ".caendat"), files_caen)

    if length(files_caen) < 1
        @warn "No caendat files found"
        return 0
    end

    files_dia = readdir(data_dir, join=true)
    filter!(x->split(x, '.')[2] == "dat", files_dia)
    files_dia = [files_dia[1]; sort(files_dia[2:end], by=x->parse(Int64, split(x, '.')[end]))]

    if n_files > length(files_caen)
        @warn "There are no $n_files caendat file(s) to be prescanned, using $(length(files_caen)) file(s) instead"
        n_files = length(files_caen)
    end

    if n_files > length(files_caen) && n_files > length(files_dia)
        n_new = maximum(length(files_caen), length(files_dia))
        @warn "There are no both $n_files caendat and diamant file(s) to be prescanned, using $(n_new) file(s) instead"
        n_files = n_new
    elseif n_files > length(files_dia)
        @warn "There are no $n_files diamant file(s) to be prescanned"
    elseif n_files > length(files_caen)
        @warn "There are no $n_files caendat file(s) to be prescanned"
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

    caen_good = true
    dia_good = true
    i_dia = 1
    if length(files_dia) > 0
        dfin = open(files_dia[1], "r")
        dsize = filesize(files_dia[1])
    else
        @warn "No diamant files found"
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
    
    prog = Progress(totsize; dt=1.0, desc="\tΔt groups ",
                     barlen=30, color=:green)
    t_caen = 0
    t_dia = 0

    last_ge = 0
    last_bgo = 0
    last_neda = 0
    last_dia = 0

    i_chunk = 0
    chunk = zeros(RawHit, chunk_size)
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
            i_chunk += 1

            if i_chunk > chunk_size
                current_pos = 0
                if isopen(cfin)
                    current_pos += position(cfin)
                end
                if isopen(dfin)
                    current_pos += position(dfin)
                end
                update!(prog, donesize + current_pos)
                sort!(chunk, by=x->(x.ts * 4 + x.tf / 256 - shift_table[x.board * 16 + x.ch + 1]))

                for chit in chunk
                    loc = chit.board * 16 + chit.ch + 1
                    if !valid_table[loc]
                        continue
                    end
                    if chit.E < Emin
                        continue
                    end
                    t = chit.ts * 4 + chit.tf / 256 - shift_table[loc]
                    if type_table[loc] == GE
                        tl = t - last_ge
                        if 0 < tl <= dt_max
                            it = round(Int64, tl / dt, RoundUp) 
                            dtgroup[1, it] += 1
                        end
                        last_ge = t
                    elseif type_table[loc] == BGO
                        tl = t - last_bgo
                        if 0 < tl <= dt_max
                            it = round(Int64, tl / dt, RoundUp) 
                            dtgroup[2, it] += 1
                        end
                        last_bgo = t
                    elseif type_table[loc] == NEDA
                        tl = t - last_neda
                        if 0 < tl <= dt_max
                            it = round(Int64, tl / dt, RoundUp) 
                            dtgroup[3, it] += 1
                        end
                        last_neda = t
                    elseif type_table[loc] == DIAMANT 
                        tl = t - last_dia
                        if 0 < tl <= dt_max
                            it = round(Int64, tl / dt, RoundUp) 
                            dtgroup[4, it] += 1
                        end
                        last_dia = t
                    end
                end
                i_chunk = 0
            else
                chunk[i_chunk] = hit
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
    dtgroup
end
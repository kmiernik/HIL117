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


"""
    test_walk(dir_name, outname)

    
    No walk correction observed and needed. Plain shift is sufficient
"""
function test_walk(dir_name, outname)

    config = TOML.parsefile("config/base.toml")
    files = readdir(dir_name, join=true)
    filter!(x->endswith(x, ".caendat"), files)

    #=
    cal_ge = zeros(3, 16)
    for label in 0:2:30
        cal_ge[:, Int(label/2) + 1] = config["label"]["$label"]["cal"]
    end
    =#

    dE = 8.0
    Estartmin = 100
    Emin = 20
    Emax = 32768
    dt_min = -1000.0
    dt_max = 1000.0
    dt = 1
    ref_board = 2
    ref_ch = 1
    walk = zeros(Int32, 16*6, round(Int64, Emax / dE, RoundUp), 
                 length(dt_min:dt:dt_max))

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
                if hit.board != ref_board || hit.ch != ref_ch
                    continue
                end
                t = hit.ts * 4 + hit.tf / 256  # tf * 4 / 1024
                if hit.E > Estartmin
                    push!(t_start, t)
                end
                #=
                ige = hit.board * 8 + round(Int, hit.ch / 2, RoundDown) + 1
                if ige == ref_ge
                    E = lin(hit.E, cal_ge[:, ige])
                    t = hit.ts * 4 + hit.tf / 256 
                    #if 500 < E
                    if abs(E - 1332.5) < dE
                        push!(t_start, t)
                    end
                end
                =#
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
                t = hit.ts * 4 + hit.tf / 256 
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
                        iE = round(Int, hit.E / dE)
                        ich = hit.board * 16 + hit.ch + 1
                        walk[ich, iE, it] += 1
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
    fout = open(outname, "w")
    write(fout, walk)
    close(fout)
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

"""
    pid - particle id
    0 - gamma
    1 - neutron
    2 - proton
    3 - alpha
"""
struct CalHit
    loc::UInt8
    E::Float64
    t::Float64
    pid::UInt8
end


function Base.zero(::Type{CalHit})
    return CalHit(zero(UInt8), 0.0, 0.0, zero(UInt8))
end

function test_scan(dir_name, config_file::String)
    config = TOML.parsefile(config_file)
    test_scan(dir_name, config)
end

"""

Test scan (caen only): E_loc and gamma_gamma 
"""
function test_scan(dir_name, config::Dict)
    files = readdir(dir_name, join=true)
    filter!(x->endswith(x, ".caendat"), files)

    last_label = 0
    for loc in keys(config["label"])
        if parse(Int64, loc) > last_label && config["label"]["$loc"]["type"] == "NEDA"
            last_label = parse(Int64, loc)
        end
    end

    ge_locs = [x for x in 1:2:32]
    cal_table = zeros(7, 16)
    fcal_table = zeros(3, 16)
    shift_table = zeros(last_label)
    for label in 1:last_label
        if label in ge_locs
            ige = round(Int64, label / 2, RoundDown) + 1
            cal_table[:, ige] = config["label"]["$label"]["cal"]
            if haskey(config["label"]["$label"], "fcal")
                fcal_table[:, ige] = config["label"]["$label"]["fcal"]
            else
                fcal_table[:, ige] = [0.0, 1.0, 0.0]
            end
        end
        shift_table[label] = config["label"]["$label"]["dt"]
    end

    period = config["spectra"]["beam_period"]
    t_min = -1000.0
    dt = 1.0
    t_max = 1000.0
    #t_loc = zeros(Int64, last_label, length(t_min:dt:t_max))

    E_min = 1
    E_max = 4096
    E_loc = zeros(Int64, last_label, length(E_min:E_max))

    gg = zeros(Int64, E_max, E_max)

    chunk_size = 10_000
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

        i_hit = 0
        i_chunk = 0
        chunk = zeros(CalHit, chunk_size)
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
                i_chunk += 1
                t = hit.ts * 4 + hit.tf / 256 - shift_table[loc]
                E = hit.E
                if loc in ge_locs
                    ige = round(Int64, loc / 2, RoundDown) + 1
                    E = quad(quadquad(
                        hit.E, cal_table[:, ige]) + (rand()-0.5), 
                             fcal_table[:, ige])
                else
                    E = E / 10
                end
                if E_min < E <= E_max
                    iE = round(Int64, E)
                    E_loc[loc, iE] += 1
                end
                if i_chunk <= chunk_size
                    chunk[i_chunk] = CalHit(UInt8(loc), E, t, zero(UInt8))
                else
                    sort!(chunk, by=x->x.t)

                    event = [1]
                    for i in 2:chunk_size
                        if chunk[i].t - chunk[event[1]].t < period
                            if chunk[i].loc in ge_locs
                                push!(event, i)
                            end
                        else
                            ev_size = length(event)
                            for j in 1:ev_size
                                if !(chunk[event[j]].loc in ge_locs)
                                    continue
                                end
                                    if E_min < chunk[event[j]].E <= E_max
                                        for k in j+1:ev_size
                                            if !(chunk[event[k]].loc in ge_locs)
                                                continue
                                            end
                                            if E_min < chunk[event[k]].E <= E_max
                                                iE1 = round(Int64, chunk[event[j]].E)
                                                iE2 = round(Int64, chunk[event[k]].E)
                                                gg[iE1, iE2] += 1
                                                gg[iE2, iE1] += 1
                                            end
                                        end
                                    end
                            end
                            event = [i]
                        end
                    end

                    ## Testing distance between two consecutive hits
                    #=
                    for j in 2:chunk_size-1
                        if E_min < chunk[j].E <= E_max
                            iE = round(Int64, chunk[j].E)
                            E_loc[chunk[j].loc, iE] += 1
                        end

                        t_near = ifelse(
                        chunk[j].t - chunk[j-1].t < chunk[j+1].t - chunk[j].t,
                        chunk[j].t - chunk[j-1].t, 
                        chunk[j].t - chunk[j+1].t)
                        if t_min < t_near < t_max
                            it = round(Int64, (t_near - t_min) / dt, RoundUp)
                            t_loc[chunk[j].loc, it] += 1
                        end

                        last = chunk[end]
                        chunk[1] = last
                        chunk[2] = CalHit(UInt8(loc), E, t)
                        i_chunk = 2
                    end
                    =#
                    last = chunk[end]
                    chunk[1] = last
                    chunk[2] = CalHit(UInt8(loc), E, t, zero(UInt8))
                    i_chunk = 2
                    dtime = (Dates.Time(Dates.now()) - time_start)
                    @printf("\r %8.2e %8.1f s ", i_hit, dtime.value * 1e-9)
                end
            end
        end
        #=
        for j in 2:i_chunk-1
            if E_min < chunk[j].E <= E_max
                iE = round(Int64, chunk[j].E)
                E_loc[chunk[j].loc, iE] += 1
            end

            t_near = ifelse(
            chunk[j].t - chunk[j-1].t < chunk[j+1].t - chunk[j].t,
            chunk[j].t - chunk[j-1].t, 
            chunk[j].t - chunk[j+1].t)
            if t_min < t_near < t_max
                it = round(Int64, (t_near - t_min) / dt, RoundUp)
                t_loc[chunk[j].loc, it] += 1
            end
        end
        =#
        event = [1]
        for i in 2:i_chunk
            if chunk[i].t - chunk[event[1]].t < period
                if chunk[i].loc in ge_locs
                    push!(event, i)
                end
            else
                ev_size = length(event)
                for j in 1:ev_size
                    if !(chunk[event[j]].loc in ge_locs)
                        continue
                    end
                    if E_min < chunk[event[j]].E <= E_max
                        for k in j+1:ev_size
                            if !(chunk[event[k]].loc in ge_locs)
                                continue
                            end
                            if E_min < chunk[event[k]].E <= E_max
                                iE1 = round(Int64, chunk[event[j]].E)
                                iE2 = round(Int64, chunk[event[k]].E)
                                gg[iE1, iE2] += 1
                                gg[iE2, iE1] += 1
                            end
                        end
                    end
                end
                event = [i]
            end
        end
    end
    println()
    E_loc, gg
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

    i_file = 0
    for filename in files
        time_start = Dates.Time(Dates.now())
        i_file += 1
        fin = open(filename, "r")
        name = split(filename, ['/'])[end]
        println("\nProcessing file $name  $i_file / $(length(files))")

        if parse(Int64, split(filename, ['/', '_', '.'])[end-1]) == 0
            header = zeros(UInt32, 12)
            read!(fin, header)
            println(header)
            agava_ts = header[4] % UInt64
            agava_ts = agava_ts << 32 + header[5] % UInt64
            @printf("AGAVA: %x\n", agava_ts)
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
    Read both caen and dia files
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

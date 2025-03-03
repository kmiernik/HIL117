"""
    distance(config, valid_table, last_label)

    Create distance_table between detectors, additional index (last)
    is for distance to the target

    returns table in distances (meters)
    
"""
function distance(config, valid_table, last_label; def_value=3.0)
    distance_table = ones(Float64, last_label+1, last_label+1) .* def_value
    distance_table[last_label+1, last_label+1] = 0.0
    for i in 1:last_label
        if !valid_table[i]
            continue
        end
        distance_table[i, i] = 0.0
        for j in i+1:last_label+1
            if j == last_label+1
                if (config["label"]["$i"]["type"] == "Ge" 
                    || config["label"]["$i"]["type"] == "BGO")
                    distance_table[i, j] = 1.0
                    distance_table[j, i] = 1.0
                elseif config["label"]["$i"]["type"] == "NEDA" 
                    distance_table[i, j] = 1.0
                    distance_table[j, i] = 1.0
                elseif config["label"]["$i"]["type"] == "DIAMANT" 
                    distance_table[i, j] = 0.0
                    distance_table[j, i] = 0.0
                end
                continue
            end
            if !valid_table[j]
                continue
            end
            if (config["label"]["$i"]["ring"] != "0"
              && config["label"]["$i"]["ring"] == config["label"]["$j"]["ring"] 
              && config["label"]["$i"]["phi"] == config["label"]["$j"]["phi"]) 
                distance_table[i, j] = 0.0
                distance_table[j, i] = 0.0
            elseif (config["label"]["$i"]["type"] == "DIAMANT"
                    && config["label"]["$j"]["type"] == "DIAMANT")
                distance_table[i, j] = 0.0
                distance_table[j, i] = 0.0
            end
        end
    end
    distance_table
end


"""
    prepare_spectra_file(config, spectranameout, pars)

    Prepare HDF5 spectra file:
    * Ge singles vs. multi 
    * Ge singles vs. particle type 
    * netron energy
    * E-loc
    * t-loc
    * Ge-Ge all/protons/alphas/no_alphas
    * dia-pid
    * neda-pid
    * multiplicity
"""
function prepare_spectra_file(config, spectranameout, pars)
    fout = h5open(spectranameout, "w")

    eventpars = ["ge_low", "bgo_low", "neda_low", "dia_low", "beam_period", 
                "ge_dt", "bgo_dt", "neda_g_dt", "neda_n_dt", "dia_dt",
                "t_delay"]
    for par in eventpars
        HDF5.attributes(fout)[par] = config["event"][par]
    end

    pidpars = ["g_low", "g_high", "n_low", "n_high",
               "a_low", "a_high", "p_low", "p_high"]
    for par in pidpars
        HDF5.attributes(fout)[par] = config["pid"][par]
    end

    data = zeros(UInt64, 
                 size(0:pars.dE:pars.Emax)[1]-1, 
                 size(1:pars.Mmax)[1]-1)

    write(fout, "gM", data)
    dset = fout["gM"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.Emax
    HDF5.attributes(dset)["ymin"] = 1
    HDF5.attributes(dset)["dy"] = 1
    HDF5.attributes(dset)["ymax"] = pars.Mmax

    # Gamma-particle type
    write(fout, "gP", data)
    dset = fout["gP"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.Emax
    HDF5.attributes(dset)["ymin"] = 1
    HDF5.attributes(dset)["dy"] = 1
    HDF5.attributes(dset)["ymax"] = pars.Mmax

    data = zeros(UInt64, size(0:pars.dE:pars.Emax)[1]-1, 4)
    write(fout, "partE", data)
    dset = fout["partE"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.Emax
    HDF5.attributes(dset)["ymin"] = 0
    HDF5.attributes(dset)["dy"] = 1
    HDF5.attributes(dset)["ymax"] = 4

    data = zeros(UInt64, 
                 size(0:pars.dE:pars.Emax)[1]-1,
                 size(1:pars.last_label)[1])
    write(fout, "Eloc", data)
    dset = fout["Eloc"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.Emax
    HDF5.attributes(dset)["ymin"] = 1
    HDF5.attributes(dset)["dy"] = 1
    HDF5.attributes(dset)["ymax"] = pars.last_label+1

    data = zeros(UInt64,
                 size(0:pars.dt:pars.tmax)[1]-1, 
                 size(1:pars.last_label)[1])
    write(fout, "tloc", data)
    dset = fout["tloc"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dt
    HDF5.attributes(dset)["xmax"] = pars.tmax
    HDF5.attributes(dset)["ymin"] = 1
    HDF5.attributes(dset)["dy"] = 1
    HDF5.attributes(dset)["ymax"] = pars.last_label+1

    data = zeros(UInt32,
                 size(0:pars.dE:pars.E2max)[1]-1,
                 size(0:pars.dE:pars.E2max)[1]-1)
    write(fout, "gg", data)
    dset = fout["gg"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.E2max
    HDF5.attributes(dset)["ymin"] = 0.0
    HDF5.attributes(dset)["dy"] = pars.dE
    HDF5.attributes(dset)["ymax"] = pars.E2max

    write(fout, "gg_alpha", data)
    dset = fout["gg_alpha"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.E2max
    HDF5.attributes(dset)["ymin"] = 0.0
    HDF5.attributes(dset)["dy"] = pars.dE
    HDF5.attributes(dset)["ymax"] = pars.E2max

    write(fout, "gg_no_alpha", data)
    dset = fout["gg_no_alpha"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.E2max
    HDF5.attributes(dset)["ymin"] = 0.0
    HDF5.attributes(dset)["dy"] = pars.dE
    HDF5.attributes(dset)["ymax"] = pars.E2max

    write(fout, "gg_proton", data)
    dset = fout["gg_proton"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.E2max
    HDF5.attributes(dset)["ymin"] = 0.0
    HDF5.attributes(dset)["dy"] = pars.dE
    HDF5.attributes(dset)["ymax"] = pars.E2max

    data = zeros(UInt64,
                 size(0:pars.dE:pars.Emax)[1]-1,
                 size(0:pars.dpid:pars.pidmax)[1]-1)

    write(fout, "dia_pid", data)
    dset = fout["dia_pid"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.Emax
    HDF5.attributes(dset)["ymin"] = 0.0
    HDF5.attributes(dset)["dy"] = pars.dpid
    HDF5.attributes(dset)["ymax"] = pars.pidmax

    write(fout, "neda_pid", data)
    dset = fout["neda_pid"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.Emax
    HDF5.attributes(dset)["ymin"] = 0.0
    HDF5.attributes(dset)["dy"] = pars.dpid
    HDF5.attributes(dset)["ymax"] = pars.pidmax

    data = zeros(UInt64,
                 size(0:pars.dt:pars.tmax)[1]-1,
                 size(0:pars.dpid:pars.pidmax)[1]-1)

    write(fout, "neda_tof", data)
    dset = fout["neda_tof"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dt
    HDF5.attributes(dset)["xmax"] = pars.tmax
    HDF5.attributes(dset)["ymin"] = 0.0
    HDF5.attributes(dset)["dy"] = pars.dpid
    HDF5.attributes(dset)["ymax"] = pars.pidmax

    data = zeros(UInt64, 
                 size(0:pars.Mmax)[1]-1,
                 size(0:pars.Mmax)[1]-1)
    write(fout, "MP", data)
    dset = fout["MP"]
    HDF5.attributes(dset)["xmin"] = 0
    HDF5.attributes(dset)["dx"] = 1
    HDF5.attributes(dset)["xmax"] = pars.Mmax
    HDF5.attributes(dset)["ymin"] = 0
    HDF5.attributes(dset)["dy"] = 1
    HDF5.attributes(dset)["ymax"] = pars.Mmax

    k = round(Int64, pars.E3max / pars.dE3)
    n = round(Int64, k * (k + 1) * (k + 2) / 6)
    data = zeros(UInt16, n, 1)

    write(fout, "ggg_prompt", data)
    dset = fout["ggg_prompt"]
    attributes(dset)["xmin"] = 0.0
    attributes(dset)["dx"] = 1.0
    attributes(dset)["xmax"] = n
    attributes(dset)["ymin"] = 0.0
    attributes(dset)["dy"] = 1.0
    attributes(dset)["ymax"] = 1.0

    close(fout)
    println("    \u25E6 spectra file ", spectranameout, " created")
end


"""
    prepare_spectra(spectranameout)

    Prepare spectra dictionary from HDF5 spectra file
"""
function prepare_spectra(spectranameout)
    fout = h5open(spectranameout, "r+")

    spectra = NamedTuple()
    for key in keys(fout)
        xmin = HDF5.attributes(fout[key])["xmin"][]
        dx = HDF5.attributes(fout[key])["dx"][]
        xmax = HDF5.attributes(fout[key])["xmax"][]
        ymin = HDF5.attributes(fout[key])["ymin"][]
        dy = HDF5.attributes(fout[key])["dy"][]
        ymax = HDF5.attributes(fout[key])["ymax"][]
        dtype = Int64
        if datatype(fout[key]) == HDF5.Datatype(HDF5.H5T_STD_U16LE)
            dtype = UInt16
        elseif datatype(fout[key]) == HDF5.Datatype(HDF5.H5T_STD_I16LE)
            dtype = Int16
        elseif datatype(fout[key]) == HDF5.Datatype(HDF5.H5T_STD_U32LE)
            dtype = UInt32
        elseif datatype(fout[key]) == HDF5.Datatype(HDF5.H5T_STD_I32LE)
            dtype = Int32
        elseif datatype(fout[key]) == HDF5.Datatype(HDF5.H5T_STD_U64LE)
            dtype = UInt64
        elseif datatype(fout[key]) == HDF5.Datatype(HDF5.H5T_IEEE_F32LE)
            dtype = Float32
        elseif datatype(fout[key]) == HDF5.Datatype(HDF5.H5T_IEEE_F64LE)
            dtype = Float64
        end
        spectra = merge(spectra, (Symbol(key) => zeros(dtype, 
                           size(xmin:dx:xmax)[1]-1, size(ymin:dy:ymax)[1]-1), ))
    end

    spectra
end

"""
    config_to_tables(config)

    Convert config to location based fast access tables 

    returns:
    type_table, valid_table, cal_table, shift_table, pid_table
"""
function config_to_tables(config, last_label)
    type_table = zeros(DetectorType, last_label)
    valid_table = zeros(Bool, last_label)
    cal_table = zeros(Float64, 2, last_label)
    shift_table = zeros(Float64, last_label)

    for label in keys(config["label"])
        loc = parse(Int64, label)
        if !config["label"]["$loc"]["valid"]
            continue
        else
            valid_table[loc] = true
            shift_table[loc] = config["label"]["$loc"]["dt"]
        end
        cal_table[:, loc] = config["label"]["$loc"]["cal"]
        if config["label"]["$loc"]["type"] == "Ge"
            type_table[loc] = GE
        elseif config["label"]["$loc"]["type"] == "BGO"
            type_table[loc] = BGO
        elseif config["label"]["$loc"]["type"] == "NEDA"
            type_table[loc] = NEDA
        elseif config["label"]["$loc"]["type"] == "DIAMANT"
            type_table[loc] = DIAMANT
        end
    end

    type_table, valid_table, cal_table, shift_table
end


function scan_run(data_dir, configfile::String, prefix; args...)
    config = TOML.parsefile(configfile)
    scan_run(data_dir, config, prefix; args...)
end


"""
    scan_run(data_dir, config::Dict, prefix; 
             chunk_size=10_000, dia_buf=100, n_prescan=2)
    scan_run(data_dir, config::Dict, prefix;
                    chunk_size=100_000, dia_buf=100, n_prescan=2, edfmode=false,
                    saveconfig=false)

    Scan run pointed by data_dir (notice that both caendat and diamant files
    must be present).
    First a prescan is performed on n_prescan.
    Both caen and dia files are read simultanuesly (keeping hits time about
    the same). Dia is read in dia_buf bunches.
    Most of the parameters are configured in the config TOML file.
    chunk_size is number of hits to be broken into events (of length equal to
    beam period).

    Returns 0 if no errors, 1 in case of error #(along with warn message)
"""
function scan_run(data_dir, config::Dict, prefix;
                    chunk_size=100_000, dia_buf=100, n_prescan=2, edfmode=false,
                    saveconfig=false)

    files_caen = readdir(data_dir, join=true)
    filter!(x->endswith(x, ".caendat"), files_caen)
    if length(files_caen) == 0
        #@warn "Could not find caendat files, aborting scan"
        return 1
    end

    run_number = parse(Int64, split(basename(files_caen[1]), 
                                    ['/', '_', '.'])[end-2])

    files_dia = readdir(data_dir, join=true)
    filter!(x->split(basename(x), '.')[2] == "dat", files_dia)
    if length(files_dia) > 0
        files_dia = [files_dia[1]; sort(files_dia[2:end], 
                            by=x->parse(Int64, split(basename(x), '.')[end]))]
    end

    println("\u25CD Run $run_number ")
    println("\u25CD Prescanning $n_prescan file(s)")
    E_loc, t_loc = prescan(data_dir, config; n_files=n_prescan)
    if n_prescan > 0
        print("    \u25E6 Calculating shifts, calibration, and period (")
        config = find_shifts(t_loc, config)
        config = find_cal(E_loc, config)
        calc_period = find_period(t_loc)
        @printf("%.3f ns)\n", calc_period)
        if abs(config["event"]["beam_period"] - calc_period) > 1.0
            @warn "Calculated period is different than the one in config!"
        end
    end
    if saveconfig
        nicetoml(config, @sprintf("config_%03d.toml", run_number))
    end

    cfin = open(files_caen[1], "r")
    csize = filesize(files_caen[1])
    i_caen = 1
    if parse(Int64, split(files_caen[1], ['/', '_', '.'])[end-1]) == 0
        header = zeros(UInt32, 12)
        read!(cfin, header)
        agava_ts = header[4] % UInt64
        agava_ts = agava_ts << 32 + header[5] % UInt64
    else
        @warn "Could not read agava time stamp"
        return 1
    end

    println("\u25CD Scanning run")
    specpars = SpectraPars(config)
    eventpars = EventPars(config)
    if edfmode
        edffile = open(prefix * ".edf", "w")
        println("    \u25E6 edf file ", prefix * ".edf", " created")
    else
        spectranameout = prefix * "_s.h5"
        if !isfile(spectranameout)
            prepare_spectra_file(config, spectranameout, specpars)
        end
        spectra = prepare_spectra(spectranameout)
    end

    i_dia = 1
    if length(files_dia) > 0
        dfin = open(files_dia[1], "r")
        dsize = filesize(files_dia[1])
    else
        dfin = IOBuffer()
        close(dfin)
        dsize = 0
        dia_good = false
    end

    totsize = 0
    donesize = 0
    for filename in files_caen
        totsize += filesize(filename)
    end
    for filename in files_dia
        totsize += filesize(filename)
    end
    prog = Progress(totsize; dt=1.0, desc="\tMain scan   ", barlen=30, 
                    color=:red)

    caen_good = true
    dia_good = true
    t_caen = 0
    t_dia = 0

    last_label = 0
    for label in keys(config["label"])
        if (parse(Int64, label) > last_label 
            && config["label"][label]["valid"])
            last_label = parse(Int64, label)
        end
    end
    
    pidpars = PidPars(config)

    type_table, valid_table, cal_table, shift_table = config_to_tables(config, last_label)

    distance_table = distance(config, valid_table, last_label)

    chunk = zeros(RawHit, chunk_size)
    i_chunk = 0
    chunk_number = 1
    last_event = Hit[]

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
                hits = read_diahit(dfin, agava_ts, dia_buf)
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
            i_chunk += 1

            if i_chunk > chunk_size
                if edfmode
                    last_event = event_builder!(chunk, last_event, 
                                                eventpars, specpars, edffile;
                                                valid_table=valid_table, 
                                                cal_table=cal_table,
                                                shift_table=shift_table,
                                                type_table=type_table,
                                                distance_table=distance_table,
                                                pidpars=pidpars)
                else
                    last_event = event_builder!(chunk, last_event, 
                                                eventpars, specpars, spectra;
                                                valid_table=valid_table, 
                                                cal_table=cal_table,
                                                shift_table=shift_table,
                                                type_table=type_table,
                                                distance_table=distance_table,
                                                pidpars=pidpars)
                end
                i_chunk = 0
                    
                chunk_number += 1
                if chunk_number % 100 == 0
                    current_pos = 0
                    if isopen(cfin)
                        current_pos += position(cfin)
                    end
                    if isopen(dfin)
                        current_pos += position(dfin)
                    end
                    update!(prog, donesize + current_pos)
                end
            else
                chunk[i_chunk] = hit
            end
        end

        if eof(cfin)
            close(cfin)
            if i_caen < length(files_caen)
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
            if i_dia < length(files_dia)
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

    if edfmode
        close(edffile)
    else
        println("\u25CD Updating spectra")
        specfile = h5open(spectranameout, "r+")
        for field in keys(spectra)
            key = String(field)
            @printf("\r    \u25E6 %20s ", key)
            dset = specfile[key]
            data = read(dset)
            data[:, :] .+= spectra[field][:, :]
            dset[:, :] = data[:, :]
            @printf("+%20d", sum(spectra[field]))
            GC.gc()
        end
        println()
        close(specfile)
    end
    println("\u25CF Done run $run_number ")

    return 0
end


"""
    scan_all(dirname, configfile; use_configs=false, n_prescan=2,
                                      skipruns=Int64[])

    dirname - base directory with data (in subdirectories, per each run)
    configfile - base config file, see below
    use_configs - if true a config directory will be search for config
                  per run e.g. config_012.toml (no prescan will be used).
                  If not found base config + prescan will be run
    n_prescan - see `scan_run`
    skipruns  - run numbers to be ommitted, default skipruns is based on 
                elog entries are removes bad data, short runs, test runs, 
                decay runs etc.
"""
function scan_all(dirname, configfile; use_configs=false, n_prescan=2, 
        edfmode=false,
        skipruns=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 21, 22, 23, 24, 26, 27, 33, 50, 72, 83, 85, 161, 162, 190, 194, 197, 198, 199, 200, 201, 202, 203, 300, 301])
    dirs = readdir(dirname, join=true)
    scan_dirs = String[]
    scan_confs = String[]
    scan_prefixes = String[]
    scan_prescans = Int64[]
    for dir in dirs
        if isfile(dir)
            continue
        end
        if !startswith(basename(dir), "run_0")
            continue
        end
        run_number = parse(Int64, split(basename(dir), ['_', '.'])[2])
        if run_number in skipruns
            continue
        end
        push!(scan_dirs, dir)
        push!(scan_prefixes, @sprintf("run_%03d", run_number))
        if use_configs
            if isfile(@sprintf("config/config_%03d.toml", run_number))
                push!(scan_confs, 
                      @sprintf("config/config_%03d.toml", run_number))
                push!(scan_prescans, 0)
            else
                push!(scan_confs, configfile)
                push!(scan_prescans, n_prescan)
            end
        else
            push!(scan_confs, configfile)
            push!(scan_prescans, n_prescan)
        end
    end
    pmap((d, c, p, r)->scan_run(d, c, p; n_prescan=r, edfmode=edfmode), 
         [d for d in scan_dirs],
         [c for c in scan_confs],
         [p for p in scan_prefixes], 
         [r for r in scan_prescans])
end

"""
Parametrs of spectra for scan 

    dE: energy step (eg. 1 keV)
    Emax: max. energy (e.g. 4000 keV)
    E2max: max. energy (e.g. 3000 keV) for gamma-gamma
    dt: time step for time spectra (1 ns)
    dpid: PID step for pid spectra (e.g. 0.01)
    pidmax: max PID (0-max), e.g. (0-1.27 -> 128 bins)
    tmax: max. time (e.g 100 ns) for time spectra
    Mmax: max. multiplicity
    ge_low: lowest energy in keV for Ge detector to accept
    bgo_low: lowest energy in channels for BGO detector to use in compton supp.
    neda_low: lowest energy in channels for NEDA detector to use
    dia_low: lowest energy in channels for DIAMANT detector to use
    last_label: last valid label
"""
struct SpectraPars
    dE::Float64
    Emax::Float64
    E2max::Float64
    dt::Float64
    tmax::Float64
    dpid::Float64
    pidmax::Float64
    Mmax::Int64
    ge_low::Float64
    bgo_low::Float64
    neda_low::Float64
    dia_low::Float64
    last_label::Int64
end


function SpectraPars(config::Dict{String, Any})
    SpectraPars(config["spectra"]["dE"],
                config["spectra"]["Emax"],
                config["spectra"]["E2max"],
                config["spectra"]["dt"],
                config["spectra"]["tmax"],
                config["spectra"]["dpid"],
                config["spectra"]["pidmax"],
                config["spectra"]["Mmax"],
                config["spectra"]["ge_low"],
                config["spectra"]["bgo_low"],
                config["spectra"]["neda_low"],
                config["spectra"]["dia_low"],
                maximum(parse.(Int64, keys(config["label"])))
               )
end


"""
    scattering!(event, hits, type_table, distance_table, specpars)

    Perform Compton suppresion and scattering suppresion on event_hits
    * type_table - table of types of detectors (on given location)
    * distance_table - table of distances between detectors 
    * specpars - SpecPars struct

"""
function scattering!(event, hits, type_table, distance_table, specpars)
    M = length(event)

    del_list = Int64[]

    for i in 1:M
        iloc = hits[event[i]].loc
        if type_table[iloc] == 2
            #BGO
            push!(del_list, i)
            for j in i+1:M
                # Search for Ge
                jloc = hits[event[j]].loc
                if (!(j in del_list) 
                    && type_table[jloc] == 1
                    && distance_table[iloc, jloc] == 0)
                    push!(del_list, j)
                end
            end
        elseif type_table[hits[event[i]].loc] == 3
            #NEDA
            for j in i+1:M
                jloc = hits[event[j]].loc
                #NEDA-NEDA gamma-scattering
                if ( !(j in del_list) 
                     && type_table[jloc] == 3
                     && hits[event[i]].pid == 1
                     && hits[event[j]].pid == 1
                     && abs(abs(hits[event[i]].t - hits[event[j]].t) 
                            - distance_table[iloc, jloc] / 0.29979 ) < 1.0)
                        push!(del_list, j)
                end
                #NEDA-NEDA neutron-scattering
                if ( !(j in del_list) 
                     && type_table[jloc] == 3
                     && hits[event[i]].pid == 2
                     && hits[event[j]].pid == 2)
                    vn = sqrt(2 * hits[event[i]].E / 939565.0) * 0.29979
                    if ( abs(abs(hits[event[i]].tof - hits[event[j]].tof) 
                             - distance_table[iloc, jloc] / vn)  < 1.0)
                        push!(del_list, j)
                    end
                end
            end
        end
    end
    deleteat!(event, sort(del_list))
end



"""
    function distance(config, valid_table, last_label)

    Create distance_table between detectors,

    returns table in distances (meters)
    
"""
function distance(config, valid_table, last_label; def_value=3.0)
    distance_table = ones(Float64, last_label, last_label) .* def_value
    for i in 1:last_label
        if !valid_table[i]
            continue
        end
        distance_table[i, i] = 0.0
        for j in i+1:last_label
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

"""
function prepare_spectra_file(config, spectranameout, pars)
    fout = h5open(spectranameout, "w")

    data = zeros(UInt32, 
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

    data = zeros(UInt32, size(0:pars.dE:pars.Emax)[1]-1, 4)
    write(fout, "partE", data)
    dset = fout["partE"]
    HDF5.attributes(dset)["xmin"] = 0.0
    HDF5.attributes(dset)["dx"] = pars.dE
    HDF5.attributes(dset)["xmax"] = pars.Emax
    HDF5.attributes(dset)["ymin"] = 0
    HDF5.attributes(dset)["dy"] = 1
    HDF5.attributes(dset)["ymax"] = 4

    data = zeros(UInt32, 
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

    data = zeros(UInt32,
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

    data = zeros(UInt32,
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

    data = zeros(UInt32,
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
    pid - particle id
    0 - unknown
    1 - gamma
    2 - neutron
    3 - proton
    4 - alpha
"""
struct Hit
    loc::UInt8
    E::Float64
    t::Float64
    tof::Float64
    pid::UInt8
end


function Hit()
    return Hit(zero(UInt8), 0.0, 0.0, 0.0, zero(UInt8))
end


function Base.zero(::Type{Hit})
    return Hit(zero(UInt8), 0.0, 0.0, 0.0, zero(UInt8))
end


"""
    event_builder!(chunk, last_event, spectra, specpars, t_event;
                       valid_table, cal_table, fcal_table, shift_table,
                       type_table, distance_table, pid_table)

    Keeps last event in chunk, swaps it to the front,
    returns i_chunk (position of last unused hit in chunk)
"""
function event_builder!(chunk, last_event, spectra, specpars, t_event;
                       valid_table, cal_table, fcal_table, shift_table,
                       type_table, distance_table, pid_table)
    d_target_neda = 1.0
    hits = [last_event;]
    sort!(chunk, by=x->(x.ts * 4 + x.tf / 256 - shift_table[x.board * 16 + x.ch + 1]))
    t_last_neda_g = 0.0
    for rawhit in chunk
        loc = rawhit.board * 16 + rawhit.ch + 1
        if loc > specpars.last_label
            continue
        elseif !valid_table[loc]
            continue
        end
        t = rawhit.ts * 4 + rawhit.tf / 256 - shift_table[loc]
        if type_table[loc] == 1
            #Ge
            ige = round(Int64, loc / 2, RoundDown) + 1
            E = quad(quadquad(
                rawhit.E, cal_table[:, loc]) + (rand()-0.5), 
                        fcal_table[:, ige])
            if E < specpars.ge_low
                continue
            end
            push!(hits, Hit(UInt8(loc), E, t, 0.0, UInt8(1)))
        elseif type_table[loc] == 2
            #BGO
            if rawhit.E < specpars.bgo_low
                continue
            end
            iE = round(Int64, rawhit.E / specpars.dE, RoundDown)
            if (1 <= iE < specpars.Emax)
                spectra.Eloc[iE, loc] += 1
            end
            push!(hits, Hit(UInt8(loc), rawhit.E, t, 0.0, UInt8(1)))
        elseif type_table[loc] == 3
            #NEDA
            E = Float64(rawhit.E)
            if E < specpars.neda_low
                continue
            end
            pid = (Float64(rawhit.qshort) + randn()) / E
            iE = round(Int64, E / specpars.dE, RoundDown)
            ip = round(Int64, pid / specpars.dpid, RoundUp)
            if (1 <= iE < specpars.Emax 
                && 1 <= ip < specpars.pidmax / specpars.dpid)
                spectra.neda_pid[iE, ip] += 1
            end
            tof = t - t_last_neda_g
            itof = round(Int64, tof / specpars.dt, RoundDown)
            if (1 <= itof < specpars.tmax 
                && 1 <= ip < specpars.pidmax / specpars.dpid)
                spectra.neda_tof[itof, ip] += 1
            end
            type = UInt8(0)
            if E > pid_table[1] && pid_table[2] <= pid < pid_table[3]
                #Neutron
                E = (5227.121 * (d_target_neda / tof)^2) * 1000 + randn()
                if E > 0.0 && !isinf(E) 
                    push!(hits, Hit(UInt8(loc), E, t_last_neda_g, 
                                    tof, UInt8(2)))
                end
            else 
                # Gamma
                push!(hits, Hit(UInt8(loc), E, t, 0.0, UInt8(1)))
                t_last_neda_g = t
            end
        elseif type_table[loc] == 4
            #DIAMANT
            E = Float64(rawhit.E)
            T = Float64(rawhit.qshort)
            pid = (T / E 
                - cal_table[1, loc] / E^2 
                - cal_table[2, loc] + 0.5)
            iE = round(Int64, E / specpars.dE, RoundDown)
            ip = round(Int64, pid / specpars.dpid, RoundUp)
            if (1 <= iE < specpars.Emax 
                && 1 <= ip < specpars.pidmax / specpars.dpid)
                spectra.dia_pid[iE, ip] += 1
            end
            type = UInt8(0)
            if pid_table[4] < pid <= pid_table[5]
                # Alpha
                push!(hits, Hit(UInt8(loc), E, t, 0.0, UInt8(4)))
            elseif pid_table[6] < pid <= pid_table[7]
                # Proton
                push!(hits, Hit(UInt8(loc), E, t, 0.0, UInt8(3)))
            end
        end
    end
    sort!(hits, by=x->x.t)
    n_hits = size(hits)[1]

    event = [1]
    for i in 2:n_hits
        if hits[i].t - hits[event[1]].t < t_event
            push!(event, i)
        else
            Mraw = length(event)
            scattering!(event, hits, type_table, distance_table, specpars)
            M = length(event)
            ges = 0
            gammas = 0
            neutrons = 0
            protons = 0
            alphas = 0
            for j in 1:M
                jloc = hits[event[j]].loc
                jE = round(Int64, hits[event[j]].E / specpars.dE, RoundDown)
                jt = round(Int64, 
                           (hits[event[j]].t - hits[event[1]].t) / specpars.dt,
                           RoundDown) + 1
                if 1 <= jt <= specpars.tmax
                    spectra.tloc[jt, jloc] += 1
                end
                if type_table[jloc] == 1
                    ges += 1
                    if 1 <= jE <= specpars.Emax
                        spectra.gP[jE, 1] += 1
                    end
                elseif type_table[jloc] == 3
                    if hits[event[j]].pid == 1
                        gammas += 1
                    elseif hits[event[j]].pid == 2
                        neutrons +=1
                        if 1 <= jE <= specpars.Emax
                            spectra.partE[jE, 2] += 1
                        end
                    end
                elseif type_table[jloc] == 4
                    if hits[event[j]].pid == 3
                        protons +=1
                        if 1 <= jE <= specpars.Emax
                            spectra.partE[jE, 3] += 1
                        end
                    elseif hits[event[j]].pid == 4
                        alphas += 1
                        if 1 <= jE <= specpars.Emax
                            spectra.partE[jE, 4] += 1
                        end
                    end
                end
            end
            for j in 1:M
                jloc = hits[event[j]].loc
                if type_table[jloc] != 1
                    continue
                end
                jE = round(Int64, hits[event[j]].E / specpars.dE, RoundDown)
                if jE < 1 || jE > specpars.Emax
                    continue
                end
                spectra.Eloc[jE, jloc] += 1
                if neutrons > 0
                    spectra.gP[jE, 2] += 1
                end
                if protons > 0
                    spectra.gP[jE, 3] += 1
                end
                if alphas > 0
                    spectra.gP[jE, 4] += 1
                end
                if 1 <= M < specpars.Mmax
                    spectra.gM[jE, M] += 1
                else
                    spectra.gM[jE, specpars.Mmax-1] += 1
                end
                if 1 <= jE <= specpars.E2max 
                    for k in j+1:M
                        kloc = hits[event[k]].loc
                        if type_table[kloc] != 1
                            continue
                        end
                        kE = round(Int64, 
                                   hits[event[k]].E / specpars.dE, RoundDown)
                        if 1 <= kE <= specpars.E2max
                            spectra.gg[jE, kE] += 1
                            spectra.gg[kE, jE] += 1
                            if protons > 0
                                spectra.gg_proton[jE, kE] += 1
                                spectra.gg_proton[kE, jE] += 1
                            end
                            if alphas > 0
                                spectra.gg_alpha[jE, kE] += 1
                                spectra.gg_alpha[kE, jE] += 1
                            end
                            if alphas == 0
                                spectra.gg_no_alpha[jE, kE] += 1
                                spectra.gg_no_alpha[kE, jE] += 1
                            end
                        end
                    end
                end
            end
            event = [i]
        end
    end
    
    left_hits = Hit[]
    for i in event
        push!(left_hits, hits[i])
    end
    
    return left_hits
end


function progress_dots(prefix, time_start, cpos, dpos,
                       i_caen, block_caen, n_caen,
                       i_dia, block_dia, n_dia)
    print("\r", prefix)
    for i in 1:i_caen-1
        print("\u25C9") 
    end
    if 0 <= cpos < 0.25
        if block_caen
            print("\u25D4")
            block_caen = false
        else
            print("\u25CB")
            block_caen = true
        end
    elseif 0.25 <= cpos < 0.5
        if block_caen
            print("\u25D1")
            block_caen = false
        else
            print("\u25D4")
            block_caen = true
        end
    elseif 0.5 <= cpos < 0.75
        if block_caen
            print("\u25D5")
            block_caen = false
        else
            print("\u25D1")
            block_caen = true
        end
    elseif cpos < 1.0
        if block_caen
            print("\u25CF")
            block_caen = false
        else
            print("\u25D5")
            block_caen = true
        end
    else
        print("\u25C9")
        block_caen = false
    end
    for i in 1:n_caen-i_caen
        print("\u25CB")
    end
    print(" ")
    for i in 1:i_dia-1
        print("\u25C9") 
    end
    if 0 <= dpos < 0.25
        if block_dia
            print("\u25D4")
            block_dia = false
        else
            print("\u25CB")
            block_dia = true
        end
    elseif 0.25 <= dpos < 0.5
        if block_dia
            print("\u25D1")
            block_dia = false
        else
            print("\u25D4")
            block_dia = true
        end
    elseif 0.5 <= dpos < 0.75
        if block_dia
            print("\u25D5")
            block_dia = false
        else
            print("\u25D1")
            block_dia = true
        end
    elseif dpos < 1.0
        if block_dia
            print("\u25CF")
            block_dia = false
        else
            print("\u25D5")
            block_dia = true
        end
    else
        print("\u25C9")
        block_dia = false
    end
    for i in 1:n_dia-i_dia
        print("\u25CB")
    end
    dtime = (Dates.Time(Dates.now()) - time_start)
    @printf(" %8.2f s ", dtime.value * 1e-9)

    return block_caen, block_dia
end


function scan_run(data_dir, configfile::String, prefix; args...)
    config = TOML.parsefile(configfile)
    scan_run(data_dir, config, prefix; args...)
end


"""
    scan_run(data_dir, config::Dict, prefix; chunk_size=10_000, n_prescan)

    Scan run pointed by data_dir (notice that both caendat and diamant files
    must be present).
    First a prescan is performed on n_prescan.
    Both caen and dia files are read simultanuesly (keeping hits time about
    the same).
    Most of the parameters are configured in the config TOML file.
    chunk_size is number of hits to be broken into events (of length equal to
    beam period).

    Returns 0 if no errors, 1 in case of error (along with warn message)
"""
function scan_run(data_dir, config::Dict, prefix; 
                    chunk_size=10_000, n_prescan=2)

    time_start = Dates.Time(Dates.now())
    files_caen = readdir(data_dir, join=true)
    filter!(x->endswith(x, ".caendat"), files_caen)
    run_number = parse(Int64, split(files_caen[1], ['/', '_', '.'])[end-2])

    files_dia = readdir(data_dir, join=true)
    filter!(x->split(x, '.')[2] == "dat", files_dia)
    files_dia = [files_dia[1]; sort(files_dia[2:end], by=x->parse(Int64, split(x, '.')[end]))]

    println("\u25CD Run $run_number ")
    println("\u25CD Prescanning $n_prescan file(s)")
    E_loc, t_loc = prescan(data_dir, config; n_files=n_prescan)
    println()
    if n_prescan > 0
        print("    \u25E6 Calculating shifts, fine calibration, and period (")
        config = find_shifts(t_loc, config)
        config = fine_cal(E_loc, config)
        config = find_period(t_loc, config)
        @printf("%.3f ns)\n", config["spectra"]["beam_period"])
    end
    nicetoml(config, @sprintf("config_%03d.toml", run_number))

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
    beam_period = config["spectra"]["beam_period"]
    specpars = SpectraPars(config)
    spectranameout = prefix * "_s.h5"
    if !isfile(spectranameout)
        prepare_spectra_file(config, spectranameout, specpars)
    end
    spectra = prepare_spectra(spectranameout)

    i_dia = 1
    dfin = open(files_dia[1], "r")
    dsize = filesize(files_dia[1])

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
    
    type_table = zeros(UInt8, last_label)
    valid_table = zeros(Bool, last_label)
    cal_table = zeros(Float64, 7, last_label)
    fcal_table = zeros(Float64, 3, 16)
    shift_table = zeros(Float64, last_label)
    pid_table = zeros(Float64, 7)
    pid_table[1] = config["pid"]["En_low"]
    pid_table[2] = config["pid"]["n_low"]
    pid_table[3] = config["pid"]["n_high"]
    pid_table[4] = config["pid"]["a_low"]
    pid_table[5] = config["pid"]["a_high"]
    pid_table[6] = config["pid"]["p_low"]
    pid_table[7] = config["pid"]["p_high"]

    for label in keys(config["label"])
        loc = parse(Int64, label)
        if !config["label"]["$loc"]["valid"]
            continue
        else
            valid_table[loc] = true
            shift_table[loc] = config["label"]["$loc"]["dt"]
        end
        if config["label"]["$loc"]["type"] == "Ge"
            type_table[loc] = 1
            cal_table[:, loc] = config["label"]["$loc"]["cal"]
            ige = round(Int64, loc / 2, RoundDown) + 1
            if haskey(config["label"]["$loc"], "fcal")
                fcal_table[:, ige] = config["label"]["$loc"]["fcal"]
            else
                @warn "No fine calibration for Ge $ige"
                fcal_table[:, ige] = [0.0, 1.0, 0.0]
            end
        elseif config["label"]["$loc"]["type"] == "BGO"
            type_table[loc] = 2
            cal_table[1:2, loc] = config["label"]["$loc"]["cal"]
        elseif config["label"]["$loc"]["type"] == "NEDA"
            type_table[loc] = 3
            cal_table[1:2, loc] = config["label"]["$loc"]["cal"]
        elseif config["label"]["$loc"]["type"] == "DIAMANT"
            type_table[loc] = 4
            cal_table[1:2, loc] = config["label"]["$loc"]["cal"]
        end
    end

    distance_table = distance(config, valid_table, last_label)

    chunk = zeros(RawHit, chunk_size)
    i_chunk = 0
    chunk_number = 1
    empty = RawHit()
    last_event = Hit[]

    block_caen = true
    block_dia = true
    while caen_good || dia_good

        hits = RawHit[]
        if (caen_good && t_dia >= t_caen) || !(dia_good)
            try
                hits = read_aggregate(cfin, config)
            catch err
                println(err)
            end
            t_caen = hits[end].ts
        elseif (dia_good && t_dia < t_caen) || !(caen_good)
            while true
                try
                    hit = read_diahit(dfin, agava_ts)
                    if hit !== empty
                        push!(hits, hit)
                        t_dia = hit.ts
                        break
                    end
                catch err
                    if eof(dfin)
                        break
                    else
                        println(err)
                    end
                end
            end
        end

        for hit in hits
            i_chunk += 1

            if i_chunk > chunk_size
                last_event = event_builder!(chunk, last_event, spectra, 
                                            specpars, beam_period;
                                        valid_table=valid_table, 
                                        cal_table=cal_table,
                                        fcal_table=fcal_table, 
                                        shift_table=shift_table,
                                        type_table=type_table,
                                        distance_table=distance_table,
                                        pid_table=pid_table)
                i_chunk = 0
                    
                chunk_number += 1
                if chunk_number % 100 == 0
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
                                        "    \u25E6 scanning : ", 
                                        time_start, cpos, dpos,
                                        i_caen, block_caen, length(files_caen),
                                        i_dia, block_dia, length(files_dia))
                end
            else
                chunk[i_chunk] = hit
            end
        end

        if eof(cfin)
            close(cfin)
            if i_caen < length(files_caen)
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
                i_dia += 1
                dfin = open(files_dia[i_dia], "r")
                dsize = filesize(files_dia[i_dia])
            else
                dia_good = false
            end
        end
    end

    progress_dots("    \u25E6 done     : ", time_start, 1.0, 1.0,
                    i_caen, block_caen, length(files_caen),
                    i_dia, block_dia, length(files_dia))
    println()
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

    return 0
end

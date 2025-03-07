"""
    All fields UInt8

    Mraw - raw multiplicity
    M - clean multiplicity (after addback, Compton suppresion and removal 
        of bad hits)
    ge - clean Ge multiplicity (addback, compton, low threshold)
    gammas - multiplicity of clean gammas in ge and neda
    neutrons - same
    protons - same
    alphas - same
"""
mutable struct EDFHeader
    Mraw::UInt8
    M::UInt8
    ge::UInt8
    gammas::UInt8
    neutrons::UInt8
    protons::UInt8
    alphas::UInt8
end


function EDFHeader()
    return EDFHeader(UInt8(0), UInt8(0), UInt8(0), UInt8(0), 
                  UInt8(0), UInt8(0), UInt8(0))
end


function Base.show(io::IO, head::EDFHeader)
    print(io, "Mraw     : ", Int(head.Mraw), "\n",
              "M        : ", Int(head.M), "\n",
              "ge       : ", Int(head.ge), "\n",
              "gammas   : ", Int(head.gammas), "\n",
              "neutrons : ", Int(head.neutrons), "\n",
              "protons  : ", Int(head.protons), "\n",
              "alphas   : ", Int(head.alphas), "\n"
              ) 
end


function Base.write(io::IO, head::EDFHeader)
    write(io, head.Mraw)
    write(io, head.M)
    write(io, head.ge)
    write(io, head.gammas)
    write(io, head.neutrons)
    write(io, head.protons)
    write(io, head.alphas)
end


function Base.read(io::IO, ::Type{EDFHeader})
    EDFHeader(read(io, UInt8), 
           read(io, UInt8),
           read(io, UInt8),
           read(io, UInt8),
           read(io, UInt8),
           read(io, UInt8),
           read(io, UInt8)
           )
end


"""
    loc: as in config
    pid: 0 for Ge, for NEDA/DIAMANT is PID alphas are negative
    E: calibrated energy in 100 ev (0-65535 -> 0-6535.5 keV), 0.1 keV resolution
    t: time in 10 ps (0-65535 -> 0-655.35 ns), 0.01 ns resolution
"""
struct EDFHit
    loc::UInt8
    pid::Int8
    E::UInt16
    t::UInt16
end


function Base.show(io::IO, hit::EDFHit)
    print(io, "loc  : ", Int(hit.loc), "\n",
              "pid  : ", Int(hit.pid), "\n",
              "E    : ", Int(hit.E), "\n",
              "t    : ", Int(hit.t), "\n")
              
end


function Base.write(io::IO, hit::EDFHit)
    write(io, hit.loc)
    write(io, hit.pid)
    write(io, hit.E)
    write(io, hit.t)
end


function Base.read(io::IO, ::Type{EDFHit})
    EDFHit(read(io, UInt8), 
        read(io, Int8),
        read(io, UInt16),
        read(io, UInt16))
end


"""
    event_builder!(chunk, last_event, eventpars, specpars, 
                        edffile::IO;
                       valid_table, cal_table, shift_table,
                       type_table, distance_table, pidpars)

    Keeps last event in chunk, swaps it to the front,
    returns i_chunk (position of last unused hit in chunk)
"""
function event_builder!(chunk, last_event, eventpars, specpars, 
                        edffile::IO;
                       valid_table, cal_table, shift_table,
                       type_table, distance_table, pidpars)
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
        if type_table[loc] !== NEDA
            t += eventpars.t_delay
        end
        if type_table[loc] == GE
            E = lin(rawhit.E, cal_table[:, loc]) + 0.5 * rand() - 0.25
        else
            E = Float64(rawhit.E)
        end
        iE = round(Int64, E / specpars.dE, RoundUp)
        if type_table[loc] == GE
            if E < eventpars.ge_low
                continue
            end
            push!(hits, Hit(UInt8(loc), E, t, zero(Int8)))
        elseif type_table[loc] == BGO
            if rawhit.E < eventpars.bgo_low
                continue
            end
            push!(hits, Hit(UInt8(loc), rawhit.E, t, zero(Int8)))
        elseif type_table[loc] == NEDA
            if E < eventpars.neda_low
                continue
            end
            pid = (Float64(rawhit.qshort) + randn()) / E
            if E > cal_table[2, loc] && pidpars.n_low <= pid < pidpars.n_high
                push!(hits, Hit(UInt8(loc), E, t, round(Int8, -pid * 100)))
            elseif (E > cal_table[1, loc] 
                    && pidpars.g_low <= pid < pidpars.g_high)
                push!(hits, Hit(UInt8(loc), E, t, round(Int8, pid * 100)))
            end
        elseif type_table[loc] == DIAMANT
            if E < eventpars.dia_low
                continue
            end
            t += (rand() - 0.5) * 5.0 
            T = Float64(rawhit.qshort)
            pid = (T / E 
                - cal_table[1, loc] / E^2 
                - cal_table[2, loc] + 0.5)
            if pidpars.a_low < pid <= pidpars.a_high
                push!(hits, Hit(UInt8(loc), E, t, round(Int8, -pid * 100)))
            elseif 0 <= pid <= 1.27
                push!(hits, Hit(UInt8(loc), E, t, round(Int8, pid * 100)))
            end
        end
    end
    sort!(hits, by=x->x.t)
    n_hits = size(hits)[1]

    event = [1]
    for i in 2:n_hits
        closeevent = false
        if type_table[hits[i].loc] == NEDA && hits[i].pid > 0
            if hits[i].t - hits[event[1]].t < eventpars.neda_g_dt
                push!(event, i)
            else
                closeevent = true
            end
        elseif type_table[hits[i].loc] == NEDA hits[i].pid < 0
            if hits[i].t - hits[event[1]].t < eventpars.neda_n_dt
                push!(event, i)
            else
                closeevent = true
            end
        elseif type_table[hits[i].loc] == DIAMANT 
            if hits[i].t - hits[event[1]].t < eventpars.dia_dt
                push!(event, i)
            else
                closeevent = true
            end
        elseif type_table[hits[i].loc] == GE
            if hits[i].t - hits[event[1]].t < eventpars.ge_dt
                push!(event, i)
            else
                closeevent = true
            end
        elseif type_table[hits[i].loc] ==  BGO
            if hits[i].t - hits[event[1]].t < eventpars.bgo_dt
                push!(event, i)
            else
                closeevent = true
            end
        end
        if closeevent
            update_spectra!(event, hits, edffile, pidpars,
                            type_table, distance_table)
            event = [i]
        end
    end
    
    left_hits = Hit[]
    for i in event
        push!(left_hits, hits[i])
    end
    
    return left_hits
end



function update_spectra!(event, hits, edffile::IO, pidpars, 
                         type_table, distance_table)
    Mraw = length(event)
    scattering!(event, hits, type_table, distance_table)
    M = length(event)
    ges = 0
    gammas = 0
    neutrons = 0
    protons = 0
    alphas = 0
    for j in 1:M
        jloc = hits[event[j]].loc
        if type_table[jloc] == GE
            ges += 1
            gammas += 1
        elseif type_table[jloc] == NEDA
            if hits[event[j]].pid < 0
                neutrons +=1
            else
                gammas += 1
            end
        elseif type_table[jloc] == DIAMANT
            if pidpars.a_low < abs(hits[event[j]].pid / 100.0) <= pidpars.a_high
                alphas += 1
            elseif pidpars.p_low < hits[event[j]].pid / 100.0 <= pidpars.p_high
                protons +=1
            end
        end
    end
    if ges > 0
        write(edffile, EDFHeader(UInt8(Mraw), UInt8(M), UInt8(ges),
                            UInt8(gammas), UInt8(neutrons),
                            UInt8(protons), UInt8(alphas)))
        for j in 1:M
            Eedf = zero(UInt16)
            if hits[event[j]].E > 6553.5
                Eedf = typemax(UInt16)
            else
                # 0.1 keV unit 100 keV -> 1000, max is 6553.5 keV
                Eedf = round(UInt16, hits[event[j]].E * 10)
            end
            tedf = zero(UInt16)
            if hits[event[j]].t - hits[event[1]].t > 655.35
                tedf = typemax(UInt16)
            else
                # 10 ps unit 10.15 ns -> 1015, max is 655.35 ns
                tedf = round(UInt16, 
                             (hits[event[j]].t - hits[event[1]].t) * 100)
            end
            write(edffile, EDFHit(hits[event[j]].loc, hits[event[j]].pid,
                                  Eedf, tedf))
        end
    end
    return 0
end


"""
    edf_scan_to_spectra(filename, configfile, fileout)

    Example of EDF scanning to spectra:
    gamma and gamma-neutron gate 

"""
function edf_scan_to_spectra(filename, configfile, fileout)
    fin = open(filename, "r")
    totsz = filesize(filename)

    config = TOML.parsefile(configfile)
    
    last_label = 0
    for label in keys(config["label"])
        if (parse(Int64, label) > last_label 
            && config["label"][label]["valid"])
            last_label = parse(Int64, label)
        end
    end

    type_table = zeros(DetectorType, last_label)
    valid_table = zeros(Bool, last_label)

    for label in keys(config["label"])
        loc = parse(Int64, label)
        if !config["label"]["$loc"]["valid"]
            continue
        else
            valid_table[loc] = true
        end
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
    
    Emax = 3096
    spectra = (g=zeros(Int64, Emax, 1),
               gg=zeros(Int64, Emax, Emax),
               ga=zeros(Int64, Emax, 1),
               gga=zeros(Int64, Emax, Emax),
               g2a=zeros(Int64, Emax, 1),
               g0a=zeros(Int64, Emax, 1),
               gg0a=zeros(Int64, Emax, Emax),
               gp0a=zeros(Int64, Emax, 1),
               ggp0a=zeros(Int64, Emax, Emax),
               type_multi=zeros(Int64, 16, 8),
               neda_pid=zeros(Int64, 127, 1024),
               dia_pid=zeros(Int64, 127, 1024))

    szhit = sizeof(EDFHit)
    t0 = Dates.Time(Dates.now())
    lastperc = 0.0
    dperc = 0.001
    prog = Progress(filesize(filename); 
                    dt=1.0, desc="\tScan to HDF5 ", barlen=30, 
                    color=:orange)

    while !eof(fin)
        head = read(fin, EDFHeader)

        if head.Mraw < 16
            spectra.type_multi[head.Mraw+1, 1] += 1
            spectra.type_multi[head.M+1, 2] += 1
            spectra.type_multi[head.ge+1, 3] += 1
            spectra.type_multi[head.gammas+1, 4] += 1
            spectra.type_multi[head.neutrons+1, 5] += 1
            spectra.type_multi[head.protons+1, 6] += 1
            spectra.type_multi[head.alphas+1, 7] += 1
        end

        gehits = EDFHit[]
        for i in 1:head.M
            hit = read(fin, EDFHit)
            if type_table[hit.loc] == GE
                push!(gehits, hit)
            elseif type_table[hit.loc] == NEDA
                E = round(Int64, hit.E / 32 + randn())
                pid = round(Int64, abs(hit.pid))
                if 1 <= E <= 1024 && 1 <= pid <= 127
                    spectra.neda_pid[pid, E] += 1
                end
            elseif type_table[hit.loc] == DIAMANT
                E = round(Int64, hit.E / 32 + randn())
                pid = round(Int64, abs(hit.pid))
                if 1 <= E <= 1024 && 1 <= pid <= 127
                    spectra.dia_pid[pid, E] += 1
                end
            end
        end
        Nge = length(gehits)
        for i in 1:Nge
            iE = round(Int64, gehits[i].E / 10, RoundDown)
            if 1 <= iE <= Emax
                spectra.g[iE, 1] += 1
                if head.alphas > 0
                    spectra.ga[iE, 1] += 1
                    if head.alphas > 1
                        spectra.g2a[iE, 1] += 1
                    end
                else
                    spectra.g0a[iE, 1] += 1
                    if head.protons > 1
                        spectra.gp0a[iE, 1] += 1
                    end
                end
                for j in i+1:Nge
                    jE = round(Int64, gehits[j].E / 10, RoundDown)
                    if 1 <= jE <= Emax
                        spectra.gg[iE, jE] += 1
                        spectra.gg[jE, iE] += 1
                        if head.alphas > 0
                            spectra.gga[iE, jE] += 1
                            spectra.gga[jE, iE] += 1
                        else
                            spectra.gg0a[iE, jE] += 1
                            spectra.gg0a[jE, iE] += 1
                            if head.protons > 1
                                spectra.ggp0a[iE, jE] += 1
                                spectra.ggp0a[jE, iE] += 1
                            end
                        end
                    end
                end
            end
        end

        pos = position(fin) / totsz
        if pos > lastperc + dperc
            t1 = Dates.Time(Dates.now())
            dt = (t1 - t0)
            update!(prog, position(fin); 
                    showvalues=[("R", 
                        round(position(fin) / dt.value * 1e3, digits=1))])
            lastperc += dperc
        end
    end
    update!(prog, filesize(filename))

    if !isfile(fileout)
        fout = h5open(fileout, "w")
        for key in keys(spectra)
            fout[String(key)] = spectra[key]
        end
        close(fout)
    else
        fout = h5open(fileout, "r+")
        for key in keys(spectra)
            @printf("\r  \u25E6 %20s ", key)
            dset = fout[String(key)]
            data = read(dset)
            data[:, :] .+= spectra[key][:, :]
            dset[:, :] = data[:, :]
            @printf("+%20d", sum(spectra[key]))
            GC.gc()
        end
        close(fout)
    end
end

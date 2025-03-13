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

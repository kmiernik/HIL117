"""
    event_builder!(chunk, last_event, spectra, specpars, t_event;
                       valid_table, cal_table, shift_table,
                       type_table, distance_table, pidpars)

    Keeps last event in chunk, swaps it to the front,
    returns i_chunk (position of last unused hit in chunk)
"""
function event_builder!(chunk, last_event, eventpars, specpars, 
                        spectra::NamedTuple;
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
        if (1 <= iE < specpars.Emax)
            spectra.Eloc[iE, loc] += 1
        end
        if type_table[loc] == GE
            if E < eventpars.ge_low
                continue
            end
            push!(hits, Hit(UInt8(loc), E, t, Int8(GAMMA)))
        elseif type_table[loc] == BGO
            if rawhit.E < eventpars.bgo_low
                continue
            end
            push!(hits, Hit(UInt8(loc), rawhit.E, t, Int8(GAMMA)))
        elseif type_table[loc] == NEDA
            if E < eventpars.neda_low
                continue
            end
            pid = (Float64(rawhit.qshort) + randn()) / E
            ip = round(Int64, pid / specpars.dpid, RoundUp)
            if (1 <= iE < specpars.Emax 
                && 1 <= ip < specpars.pidmax / specpars.dpid)
                spectra.neda_pid[iE, ip] += 1
            end
            tof = t - t_last_neda_g
            if (1 <= tof < specpars.tmax 
                && 1 <= ip < specpars.pidmax / specpars.dpid)
                itof = round(Int64, tof / specpars.dt, RoundUp)
                spectra.neda_tof[itof, ip] += 1
            end
            if E > cal_table[2, loc] && pidpars.n_low <= pid < pidpars.n_high
                push!(hits, Hit(UInt8(loc), E, t, Int8(NEUTRON)))
            elseif (E > cal_table[1, loc] 
                    && pidpars.g_low <= pid < pidpars.g_high)
                push!(hits, Hit(UInt8(loc), E, t, Int8(GAMMA)))
                t_last_neda_g = t
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
            ip = round(Int64, pid / specpars.dpid, RoundUp)
            if (1 <= iE < specpars.Emax 
                && 1 <= ip < specpars.pidmax / specpars.dpid)
                spectra.dia_pid[iE, ip] += 1
            end
            if pidpars.a_low < pid <= pidpars.a_high
                push!(hits, Hit(UInt8(loc), E, t, Int8(ALPHA)))
            elseif pidpars.p_low < pid <= pidpars.p_high
                push!(hits, Hit(UInt8(loc), E, t, Int8(PROTON)))
            end
        end
    end
    sort!(hits, by=x->x.t)
    n_hits = size(hits)[1]

    event = [1]
    for i in 2:n_hits
        closeevent = false
        if (type_table[hits[i].loc] == NEDA 
            && ParticleType(hits[i].pid) == GAMMA)
            if hits[i].t - hits[event[1]].t < eventpars.neda_g_dt
                push!(event, i)
            else
                closeevent = true
            end
        elseif (type_table[hits[i].loc] == NEDA 
                && ParticleType(hits[i].pid) == NEUTRON)
            if hits[i].t - hits[event[1]].t < eventpars.neda_n_dt
                push!(event, i)
            else
                closeevent = true
            end
        elseif (ParticleType(hits[i].pid) == PROTON
                || ParticleType(hits[i].pid) == ALPHA)
            if hits[i].t - hits[event[1]].t < eventpars.dia_dt
                push!(event, i)
            else
                closeevent = true
            end
        elseif (type_table[hits[i].loc] == GE 
                && ParticleType(hits[i].pid) == GAMMA)
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
            update_spectra!(event, hits, spectra, type_table, distance_table, specpars)
            event = [i]
        end
    end
    
    left_hits = Hit[]
    for i in event
        push!(left_hits, hits[i])
    end
    
    return left_hits
end


"""
    scattering!(event, hits, type_table, distance_table)

    Perform Compton suppresion and scattering suppresion on event_hits
    * type_table - table of types of detectors (on given location)
    * distance_table - table of distances between detectors 
    * specpars - SpecPars struct

"""
function scattering!(event, hits, type_table, distance_table)
    M = length(event)

    del_list = Int64[]

    for i in 1:M
        iloc = hits[event[i]].loc
        if type_table[iloc] == BGO
            push!(del_list, i)
            for j in i+1:M
                jloc = hits[event[j]].loc
                if (!(j in del_list) 
                    && type_table[jloc] == GE
                    && distance_table[iloc, jloc] == 0)
                    push!(del_list, j)
                end
            end
        elseif type_table[hits[event[i]].loc] == NEDA
            for j in i+1:M
                jloc = hits[event[j]].loc
                if ( !(j in del_list) 
                     && type_table[jloc] == NEDA
                     && ((hits[event[i]].pid < 0 && hits[event[j]].pid < 0)
                         || (hits[event[i]].pid == Int8(NEUTRON)
                             && hits[event[j]].pid == Int8(NEUTRON))))
                    # Neighboring detectors are removed if dt > 2.0 ns
                    if (0.14 < distance_table[iloc, jloc] < 0.15
                        &&  abs(hits[event[j]].t - hits[event[i]].t) > 2.0)
                        push!(del_list, j)
                    end
                end
            end
        end
    end
    deleteat!(event, sort(del_list))
end


function update_spectra!(event, hits, spectra::NamedTuple, type_table, distance_table, specpars)
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
        jE = round(Int64, hits[event[j]].E / specpars.dE, RoundUp)
        jt = round(Int64, 
                    (hits[event[j]].t - hits[event[1]].t) / specpars.dt,
                    RoundUp)
        if 1 <= jt <= specpars.tmax
            spectra.tloc[jt, jloc] += 1
        end
        if type_table[jloc] == GE
            ges += 1
            if 1 <= jE <= specpars.Emax
                spectra.gP[jE, 1] += 1
            end
        elseif type_table[jloc] == NEDA
            if ParticleType(hits[event[j]].pid) == GAMMA
                gammas += 1
            elseif ParticleType(hits[event[j]].pid) == NEUTRON
                neutrons +=1
                if 1 <= jE <= specpars.Emax
                    spectra.partE[jE, 2] += 1
                end
            end
        elseif type_table[jloc] == DIAMANT
            if ParticleType(hits[event[j]].pid) == PROTON
                protons +=1
                if 1 <= jE <= specpars.Emax
                    spectra.partE[jE, 3] += 1
                end
            elseif ParticleType(hits[event[j]].pid) == ALPHA
                alphas += 1
                if 1 <= jE <= specpars.Emax
                    spectra.partE[jE, 4] += 1
                end
            end
        end
    end
    for (itype, m) in enumerate([Mraw, M, ges, gammas, neutrons,
                                    protons, alphas])
        iM = ifelse(m < specpars.Mmax-2, m+1, specpars.Mmax-1)
        spectra.MP[iM, itype] += 1
    end
    if ges < 1
        return 0
    end

    for j in 1:M
        jloc = hits[event[j]].loc
        Ej = hits[event[j]].E
        if type_table[jloc] != GE
            continue
        end
        if Ej < 1 || Ej > specpars.Emax
            continue
        end
        jE = round(Int64, Ej / specpars.dE, RoundUp)
        if neutrons == 0 && protons == 0 && alphas == 0
            spectra.gP[jE, 1] += 1
        end
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
        if 1 <= Ej <= specpars.E2max 
            for k in j+1:M
                kloc = hits[event[k]].loc
                if type_table[kloc] != GE
                    continue
                end
                Ek = hits[event[k]].E
                if 1 <= Ek <= specpars.E2max
                    kE = round(Int64, 
                                hits[event[k]].E / specpars.dE, RoundUp)
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
                for m in k+1:M
                    mloc = hits[event[m]].loc
                    if type_table[mloc] != GE
                        continue
                    end
                    Em = hits[event[m]].E
                    if (   Ej <= specpars.E3max && Ek <= specpars.E3max
                        && 1 <= Em <= specpars.E3max)
                        j3E = round(Int64, Ej / specpars.dE3, RoundUp)
                        k3E = round(Int64, Ek / specpars.dE3, RoundUp)
                        m3E = round(Int64, Em / specpars.dE3, RoundUp)
                        ic = icube(sort([j3E, k3E, m3E])...)
                        spectra.ggg_prompt[ic] += 1
                    end
                end
            end
        end
    end
    return 0
end

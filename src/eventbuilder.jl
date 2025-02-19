"""
    event_builder!(chunk, last_event, spectra, specpars, t_event;
                       valid_table, cal_table, shift_table,
                       type_table, distance_table, pidpars)

    Keeps last event in chunk, swaps it to the front,
    returns i_chunk (position of last unused hit in chunk)
"""
function event_builder!(chunk, last_event, spectra, specpars, t_event;
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
        if type_table[loc] == GE
            E = lin(rawhit.E, cal_table[:, loc]) + 0.5 * rand() - 0.25
        else
            E = Float64(rawhit.E)
        end
        iE = round(Int64, E / specpars.dE, RoundDown)
        if (1 <= iE < specpars.Emax)
            spectra.Eloc[iE, loc] += 1
        end
        if type_table[loc] == GE
            if E < specpars.ge_low
                continue
            end
            push!(hits, Hit(UInt8(loc), E, t, 0.0, GAMMA))
        elseif type_table[loc] == BGO
            if rawhit.E < specpars.bgo_low
                continue
            end
            push!(hits, Hit(UInt8(loc), rawhit.E, t, 0.0, GAMMA))
        elseif type_table[loc] == NEDA
            if E < specpars.neda_low
                continue
            end
            pid = (Float64(rawhit.qshort) + randn()) / E
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
            if E > cal_table[2, loc] && pidpars.n_low <= pid < pidpars.n_high
                E = (5227.121 * (d_target_neda / tof)^2) * 1000 + randn()
                if E > 0.0 && !isinf(E) 
                    push!(hits, Hit(UInt8(loc), E, t_last_neda_g, 
                                    tof, NEUTRON))
                end
            elseif (E > cal_table[1, loc] 
                    && pidpars.g_low <= pid < pidpars.g_high)
                push!(hits, Hit(UInt8(loc), E, t, 0.0, GAMMA))
                t_last_neda_g = t
            end
        elseif type_table[loc] == DIAMANT
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
                push!(hits, Hit(UInt8(loc), E, t, 0.0, ALPHA))
            elseif pidpars.p_low < pid <= pidpars.p_high
                push!(hits, Hit(UInt8(loc), E, t, 0.0, PROTON))
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
                if type_table[jloc] == GE
                    ges += 1
                    if 1 <= jE <= specpars.Emax
                        spectra.gP[jE, 1] += 1
                    end
                elseif type_table[jloc] == NEDA
                    if hits[event[j]].pid == GAMMA
                        gammas += 1
                    elseif hits[event[j]].pid == NEUTRON
                        neutrons +=1
                        if 1 <= jE <= specpars.Emax
                            spectra.partE[jE, 2] += 1
                        end
                    end
                elseif type_table[jloc] == DIAMANT
                    if hits[event[j]].pid == PROTON
                        protons +=1
                        if 1 <= jE <= specpars.Emax
                            spectra.partE[jE, 3] += 1
                        end
                    elseif hits[event[j]].pid == ALPHA
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
                event = [i]
                continue
            end

            for j in 1:M
                jloc = hits[event[j]].loc
                if type_table[jloc] != GE
                    continue
                end
                jE = round(Int64, hits[event[j]].E / specpars.dE, RoundDown)
                if jE < 1 || jE > specpars.Emax
                    continue
                end
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
                if 1 <= jE <= specpars.E2max 
                    for k in j+1:M
                        kloc = hits[event[k]].loc
                        if type_table[kloc] != GE
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
                     && hits[event[i]].pid == GAMMA
                     && hits[event[j]].pid == GAMMA
                     && abs(abs(hits[event[i]].t - hits[event[j]].t) 
                            - distance_table[iloc, jloc] / 0.29979 ) < 1.0)
                        push!(del_list, j)
                end
                if ( !(j in del_list) 
                     && type_table[jloc] == NEDA
                     && hits[event[i]].pid == NEUTRON
                     && hits[event[j]].pid == NEUTRON)
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

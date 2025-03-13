using HIL117
using HDF5
using Printf
using ProgressMeter
using TOML


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


function edf_gamma_spectra(dirname, configfile, fileout)
    files = readdir(dirname, join=true)
    filter!(x->endswith(x, ".edf"), files)
    config = TOML.parsefile(configfile)
    last_label = 0
    for label in keys(config["label"])
        if (parse(Int64, label) > last_label 
            && config["label"][label]["valid"])
            last_label = parse(Int64, label)
        end
    end

    type_table = zeros(HIL117.DetectorType, last_label)
    valid_table = zeros(Bool, last_label)

    for label in keys(config["label"])
        loc = parse(Int64, label)
        if !config["label"]["$loc"]["valid"]
            continue
        else
            valid_table[loc] = true
        end
        if config["label"]["$loc"]["type"] == "Ge"
            type_table[loc] = HIL117.GE
        elseif config["label"]["$loc"]["type"] == "BGO"
            type_table[loc] = HIL117.BGO
        elseif config["label"]["$loc"]["type"] == "NEDA"
            type_table[loc] = HIL117.NEDA
        elseif config["label"]["$loc"]["type"] == "DIAMANT"
            type_table[loc] = HIL117.DIAMANT
        end
    end


    Emax = 3096
    spectra = (g_npa=zeros(Int64, Emax, 1000),)
    
    totsz = 0
    for filename in files
        totsz += filesize(filename)
    end
    donesz = 0
    prog = Progress(totsz; 
                    dt=1.0, desc="\tScan to HDF5 ", barlen=30, 
                    color=:orange)
    szhit = sizeof(HIL117.EDFHit)
    dperc = 0.001

    for filename in files
        fin = open(filename, "r")

        lastperc = 0.0

        while !eof(fin)
            head = read(fin, HIL117.EDFHeader)


            gehits = HIL117.EDFHit[]
            for i in 1:head.M
                hit = read(fin, HIL117.EDFHit)
                if type_table[hit.loc] == HIL117.GE
                    push!(gehits, hit)
                end
            end

            pid = 1
            if head.neutrons >= 9
                pid += 9
            else 
                pid += head.neutrons
            end
            if head.protons >= 9
                pid += 90
            else 
                pid += 10 * head.protons
            end
            if head.alphas >= 9
                pid += 900
            else
                pid += 100 * head.alphas
            end
            Nge = length(gehits)
            for i in 1:Nge
                iE = round(Int64, gehits[i].E / 10, RoundDown)
                if 1 <= iE <= Emax
                    spectra.g_npa[iE, pid] += 1
                end
            end

            pos = position(fin) / totsz
            if pos > lastperc + dperc
                update!(prog, position(fin) + donesz)
                lastperc += dperc
            end
        end
        donesz += filesize(filename)
        update!(prog, donesz)
    end

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

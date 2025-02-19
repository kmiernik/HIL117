"""
    bareconfig(tomlfile)
    Create bare TOML config file
    
"""
function zeroconfig(tomlfile)
    fout = open(tomlfile, "w")
    header = """
title = "Zero config for experiment HIL117

[boards]
    PHA = [0, 1]
    PSD = [2, 3, 4, 5]
    DIA = [6, 7, 8, 9]

[spectra]
	dE = 1.0
	Emax = 4096.0
    E2max = 2560.0
	dt = 1.0
	tmax = 100.0
    dpid = 0.01
    pidmax = 1.27
	Mmax = 10
	ge_low = 30.0
	bgo_low = 30.0
	neda_low = 10.0
	dia_low = 100.0
	beam_period = 73.89

"""
    print(fout, header)
    print(fout, "[label]\n")
    #TODO: Ring (theta) / Phi
    R101 = [[2, 10], [0, 12]]
    R117 = [[0, 6, 12], [2, 10, 14]]
    R143 = [[4, 8], [4, 6, 8]]
    for board in [0, 1]
        for ch in 0:15
            if board == 0 && ch >= 14
                continue
            end
            print(fout, "\t[label.", board*16 + ch + 1, "]\n")
            print(fout, "\t\ttype = \"", ifelse(ch % 2 == 0, "Ge", "BGO"), "\"\n")
            ring = ""
            if 2 * round(Int64, ch/2, RoundDown) in R101[board+1]
                ring = "101"
            elseif 2 * round(Int64, ch/2, RoundDown) in R117[board+1]
                ring = "117"
            elseif 2 * round(Int64, ch/2, RoundDown) in R143[board+1]
                ring = "143"
            end
            print(fout, "\t\tring = \"", ring, "\"\n")
            print(fout, "\t\tphi = \"", 0, "\"\n")
            print(fout, "\t\tdt = 0.0\n")
            print(fout, "\t\tcal = [0.0, 1.0, 0.0]\n")
            print(fout, "\t\tvalid = true\n")
            print(fout, "\t\tcomment = \" \"\n")
            print(fout, "\n")
        end
    end
    for board in [2, 3, 4, 5]
        for ch in 0:15
            type = "NEDA"
            if board == 5 
                if ch == 8
                    type = "Sum"
                elseif ch == 12
                    type = "Trig"
                elseif ch > 3
                    continue
                end
            end
            print(fout, "\t[label.", board*16 + ch + 1, "]\n")
            print(fout, "\t\ttype = \"", type, "\"\n")
            ring = ""
            print(fout, "\t\tring = \"", ring, "\"\n")
            print(fout, "\t\tphi = \"", 0, "\"\n")
            print(fout, "\t\tdt = 0.0\n")
            print(fout, "\t\tcal = [0.0, 1.0, 0.0]\n")
            print(fout, "\t\tvalid = true\n")
            print(fout, "\t\tcomment = \" \"\n")
            print(fout, "\n")
        end
    end
    for board in [6, 7, 8, 9]
        for ch in 0:15
            type = "DIAMANT"
            print(fout, "\t[label.", board*16 + ch + 1, "]\n")
            print(fout, "\t\ttype = \"", type, "\"\n")
            ring = ""
            print(fout, "\t\tring = \"", ring, "\"\n")
            print(fout, "\t\tphi = \"", 0, "\"\n")
            print(fout, "\t\tdt = 0.0\n")
            print(fout, "\t\tcal = [0.0, 1.0, 0.0]\n")
            print(fout, "\t\tvalid = true\n")
            print(fout, "\t\tcomment = \" \"\n")
            print(fout, "\n")
        end
    end
    close(fout)
end


"""
    nicetoml(config, tomlfile, additional_comment="")

    Create nice TOML config file (sorted by location)

"""
function nicetoml(config, tomlfile, additional_comment="")
    labels = sort(collect(config["label"]), by=x->parse(Int64, x.first))
    fout = open(tomlfile, "w")

    print(fout, "title = \"", config["title"], "\"\n")
    comment = config["comment"] * " " * additional_comment
    print(fout, "comment = \"$comment\"\n")
    print(fout, "\n")

    print(fout, "[boards]\n")
    print(fout, "    PHA = ", config["boards"]["PHA"], "\n")
    print(fout, "    PSD = ", config["boards"]["PSD"], "\n")
    print(fout, "    DIA = ", config["boards"]["DIA"], "\n")
    print(fout, "\n")

    print(fout, "[spectra]\n")
    parnames = ["dE", "Emax", "E2max", "dt", "tmax", "dpid", "pidmax",
                "Mmax", "ge_low", "bgo_low", "neda_low", "dia_low",
                "beam_period"]
    for par in parnames
        print(fout, "\t$par = ", config["spectra"][par], "\n")
    end
    print(fout, "\n")

    print(fout, "[pid]\n")
    pidnames = ["g_low", "g_high", 
                "n_low", "n_high", 
                "a_low", "a_high",
                "p_low", "p_high"]
    for par in pidnames
        print(fout, "\t$par = ", config["pid"][par], "\n")
    end
    print(fout, "\n")

    for entry in labels
        print(fout, "[label.", entry.first, "]\n")
        print(fout, "\ttype = \"", entry.second["type"], "\"\n")
        print(fout, "\tring = \"", entry.second["ring"], "\"\n")
        print(fout, "\tphi = \"", entry.second["phi"], "\"\n")
        print(fout, "\tdt = ", round(entry.second["dt"], digits=3), "\n")
        if haskey(entry.second, "rcal")
            print(fout, "\trcal = ", entry.second["rcal"], "\n")
        elseif entry.first in 1:2:32
            print(fout, "\trcal = [0.0, 1.0]\n")
        end
        if haskey(entry.second, "cal")
            print(fout, "\tcal = ", entry.second["cal"], "\n")
        elseif entry.first in 1:2:32
            print(fout, "\tcal = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0]\n")
        else
            print(fout, "\tcal = [0.0, 1.0]\n")
        end
        if haskey(entry.second, "fcal")
            print(fout, "\tfcal = ", entry.second["fcal"], "\n")
        end
        print(fout, "\tvalid = ", entry.second["valid"], "\n")
        if haskey(entry.second, "comment")
            print(fout, "\tcomment = \"", entry.second["comment"], "\"\n")
        else
            print(fout, "\tcomment = \" \"\n")
        end
        print(fout, "\n")
    end
    close(fout)
end


function progress_dots(prefix, time_start, cpos, dpos,
                       i_caen, block_caen, n_caen,
                       i_dia, block_dia, n_dia)
    print("\r", prefix)
    for i in 1:i_caen-1
        print("\u25C9") 
    end
    if 0 < cpos < 0.25
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
    elseif 0.75 <= cpos < 1.0
        if block_caen
            print("\u25CF")
            block_caen = false
        else
            print("\u25D5")
            block_caen = true
        end
    elseif cpos == 1.0
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
    if 0 < dpos < 0.25
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
    elseif 0.75 <= dpos < 1.0
        if block_dia
            print("\u25CF")
            block_dia = false
        else
            print("\u25D5")
            block_dia = true
        end
    elseif dpos == 1.0
        print("\u25C9")
        block_dia = false
    end
    for i in 1:n_dia-i_dia
        print("\u25CB")
    end
    dtime = (Dates.Time(Dates.now()) - time_start)
    @printf(" %-8.1f s ", dtime.value * 1e-9)

    return block_caen, block_dia
end



"""
    bareconfig(tomlfile)
    Create bare TOML config file

"""
function bareconfig(tomlfile)
    fout = open(tomlfile, "w")
    header = """
title = "Zero config for experiment HIL117

[boards]
    PHA = [0, 1]
    PSD = [2, 3, 4, 5]

[spectra]
    dE = 1.0
    Emax = 4096.0
    dt = 1.0
    tmax = 100.0
    Mmax = 10
    coin = [40.0, 0.0, 0.0]
    bgo_low = 30.0
    ge_low = 10.0
    neda_low = 10.0
    dia_low = 10.0

"""
    print(fout, header)
    print(fout, "[label]\n")
    #TODO: Ring (theta) / Phi
    R101 = [[2, 10], [0, 8, 12]]
    R117 = [[0, 6, 12], [2, 10]]
    R143 = [[4, 8], [4, 6, 10]]
    for board in [0, 1]
        for ch in 0:15
            if board == 0 && ch >= 14
                continue
            end
            print(fout, "\t[label.", board*16 + ch, "]\n")
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
            print(fout, "\t\twalk = [0.0, 1.0, 0.0]\n")
            print(fout, "\t\tvalid = true\n")
            print(fout, "\t\tcomment = \" \"\n")
            print(fout, "\n")
        end
    end
    for board in [2, 3, 4, 5]
        for ch in 0:15
            print(fout, "\t[label.", board*16 + ch, "]\n")
            type = "NEDA"
            if board == 5 
                if ch == 8
                    type = "Sum"
                elseif ch == 12
                    type = "Trig"
                end
            end
            print(fout, "\t\ttype = \"", type, "\"\n")
            ring = ""
            print(fout, "\t\tring = \"", ring, "\"\n")
            print(fout, "\t\tphi = \"", 0, "\"\n")
            print(fout, "\t\tdt = 0.0\n")
            print(fout, "\t\tcal = [0.0, 1.0, 0.0]\n")
            print(fout, "\t\twalk = [0.0, 1.0, 0.0]\n")
            print(fout, "\t\tvalid = true\n")
            print(fout, "\t\tcomment = \" \"\n")
            print(fout, "\n")
        end
    end
    close(fout)
end

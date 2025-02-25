"""
    lin(x, p)

    Linear function
        p[1] + p[2] * x
"""
function lin(x, p)
    @. p[1] + p[2] * x
end


"""
    quad(x, p)

    Quadratic function
        p[1] + p[2] * x + p[3] * x^2
"""
function quad(x, p)
    @. p[1] + p[2] * x + p[3] * x^2
end


"""
    gauss(x, p)

    Gaussian shape function with linear background
    (p[1] / sqrt(2 * pi * p[3]^2) * exp(-(x - p[2])^2 / (2 * p[3]^2)) 

    * p[1]: Area (A)
    * p[2]: Mean (μ)
    * p[3]: Deviation (σ)
"""
function gauss(x, p)
    @. p[1] / sqrt(2 * pi * p[3]^2) * exp(-(x - p[2])^2 / (2 * p[3]^2)) 
end


"""
    gausslin(x, p)

    Gaussian shape function with linear background
"""
function gausslin(x, p)
    gauss(x, p[1:3]) .+ p[4] .+ p[5] * x
end


"""
    rebin(data, bins)
    
    Rebin 1-D data using bin size `bins`
"""
function rebin(data, bins)
    sy = round(Int64, size(data)[1] / bins, RoundUp)
    datac = copy(data)
    if sy * bins > size(data)[1]
        append!(datac, zeros(eltype(data), sy * bins - size(data)[1]))
    end
    collect(mean(reshape(datac, (bins, round(Int64, sy))), dims=1)')
end


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

[event]
	ge_low = 30.0
	bgo_low = 30.0
	neda_low = 10.0
	dia_low = 100.0
	beam_period = 74.0
	ge_dt = 148.0
	bgo_dt = 74.0
	neda_g_dt = 10.0
	neda_n_dt = 37.0
	dia_dt = 37.0
    t_delay = 20.0

[pid]
	g_low = 0.71
	g_high = 0.95
	n_low = 0.5
	n_high = 0.69
	a_low = 0.42
	a_high = 0.56
	p_low = 0.56
	p_high = 0.67

[calibration]
	E1 = 511.0
	E2 = 1256.69

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
    specnames = ["dE", "Emax", "E2max", "dt", "tmax", "dpid", "pidmax",
                "Mmax"]
    for par in specnames
        print(fout, "\t$par = ", config["spectra"][par], "\n")
    end
    print(fout, "\n")

    print(fout, "[event]\n")
    evnames = ["ge_low", "bgo_low", "neda_low", "dia_low", "beam_period", 
                "ge_dt", "bgo_dt", "neda_g_dt", "neda_n_dt", "dia_dt",
                "t_delay"]
    for par in evnames
        print(fout, "\t$par = ", config["event"][par], "\n")
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

    print(fout, "[calibration]\n")
    calnames = ["E1", "E2"]
    for par in calnames
        print(fout, "\t$par = ", config["calibration"][par], "\n")
    end
    print(fout, "\n")

    for entry in labels
        print(fout, "[label.", entry.first, "]\n")
        print(fout, "\ttype = \"", entry.second["type"], "\"\n")
        print(fout, "\tring = \"", entry.second["ring"], "\"\n")
        print(fout, "\tphi = \"", entry.second["phi"], "\"\n")
        print(fout, "\tdt = ", round(entry.second["dt"], digits=3), "\n")
        if haskey(entry.second, "cal")
            print(fout, "\tcal = ", entry.second["cal"], "\n")
        else
            print(fout, "\tcal = [0.0, 1.0]\n")
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


"""
    Modify prescanned config files for each run.
    Change function according to needs.
"""
function modify_config_bulk()
	files = readdir("config/", join=true)
	filter!(x->startswith(basename(x), "config_"), files)
	for file in files
	    config = TOML.parsefile(file)
	    config["event"]["dia_dt"] = 128.0
	    #config["calibration"] = Dict{String, Any}()
	    #config["calibration"]["E1"] = 511.0
	    #config["calibration"]["E2"] = 1256.69
	    nicetoml(config, file)
	end
end


"""
    sum_spectra(data_dir, sumfilename; prefix="", runs=Int64[])

    Add all spectra in 'data_dir' into new file 'sumfilename'. The new file
    keeps the structure of spectra files.
    * prefix - default="", files name prefix , e.g. "run_19_s.h5" - "run_" 
    * runs - default=[], if empty list it will use all files; if some values
        are given (array of Int64) it will scan only the given files
"""
function sum_spectra(data_dir, sumfilename; prefix="", skip=Int64[])

    data_path = readdir(data_dir, join=true)
    filter!(p->endswith(p, "_s.h5"), data_path)
    filter!(p->startswith(split(p, '/')[end], prefix), data_path)
    if length(skip) > 0
        filter!(p -> (isnothing(tryparse(Int64, replace(split(data_path[1], '/')[end], "_s.h5"=>"", ""=>"")))) || !(parse(Int64, replace(split(data_path[1], '/')[end], "_s.h5"=>"", ""=>"")) in skip), data_path)
    end
    fn = length(data_path)

    if fn == 0
        return fn
    end

    prog = Progress(length(data_path); dt=1.0, desc="Summing files ", barlen=30)

    if !isfile(sumfilename)
        fout = h5open(sumfilename, "w")
        fin = h5open(data_path[1], "r")
        for key in keys(fin)
            inattr = attrs(fin[key])
            indset = fin[key]
            
            create_dataset(fout, key, datatype(fin[key]), size(indset))
            outdset = fout[key]
            for at_key in inattr
                attributes(outdset)[at_key.first] = at_key.second
            end
        end
        close(fin)
        close(fout)
    end
    fout = h5open(sumfilename, "r+")

    i_file = 1
    for path in data_path
        fin = h5open(path, "r")
        for key in keys(fin)
            indset = fin[key]
            outdset = fout[key]

            indata = read(indset)
            outdata = read(outdset)
            outdata[:, :] .+= indata[:, :]
            outdset[:, :] = outdata[:, :]
            GC.gc()
        end
        close(fin)
        update!(prog, i_file)
        i_file += 1
    end
    close(fout)
    return fn
end



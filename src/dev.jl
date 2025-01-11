
function test_co60()
    data = read_raw_data("data/needleRU_i1329_0004_0000.caendat"; n_agg=100000)

    spectra = Dict{Int, Dict{Int, Vector{Int64}}}()
    for board in 0:5
        chans = Dict{Int, Vector{Int64}}()
        for ch in 0:15
            chans[ch] = zeros(Int64, 16384)
        end
        spectra[board] = chans
    end

    for hit in data
        if 1 < hit.E <= 16384
            spectra[hit.board][hit.ch][hit.E] += 1
        end
    end
    spectra
end

function test_times()
    # 0-5 boards, 0-15 channels
    # x = board * 16 + channel + 1
    tmax = 10_000_000
    dt = 10.0
    spectra = zeros(UInt16, 6*16, Int(tmax / dt))
    #times = Float64[]

    fin = open("data/needleRU_i1363_0036_0000.caendat", "r")
    config = TOML.parsefile("config/base.toml")
    header = zeros(UInt32, 12)
    read!(fin, header)

    t0 = -1.0
    while !eof(fin) 
        hits = read_aggregate(fin, config)
        for hit in hits
            t = hit.ts * 4 + hit.tf / 1024 
            if t0 == -1.0
                t0 = t
            else
                t -= t0
            end
            if 1 < t <= tmax
                spectra[hit.board * 16 + hit.ch + 1, 
                        round(Int64, t / dt, RoundDown)] += 1
            end
            #=
            if hit.board == 2
                push!(times, t)
            end
            =#
        end
    end
    spectra
    #times
end


"""

    Rough Ge calibration based on Co-60

    Returns vector of two elements vector
    Ge B0C0, B0C2, ..., B0C14, B1C0, ..., B1C12, B1C14
        1     2    ....  8      9          15    16
    ige = board * 8 + round(Int, ch / 2, RoundDown) + 1
    
"""
function init_cal()
    data = read_raw_data("data/needleRU_i1329_0004_0000.caendat"; n_agg=100000)

    spectra = Dict{Int, Dict{Int, Vector{Int64}}}()
    for board in 0:5
        chans = Dict{Int, Vector{Int64}}()
        for ch in 0:15
            chans[ch] = zeros(Int64, 16384)
        end
        spectra[board] = chans
    end
    for hit in data
        if 1 < hit.E <= 16384
            spectra[hit.board][hit.ch][hit.E] += 1
        end
    end

    cal_ge = Vector{Vector{Float64}}()
    for b in 0:1
        for c in 0:2:15
            push!(cal_ge, [0.0, 0.0])
        end
    end
    thr = 20
    for b in 0:1
        for c in 0:2:15
            fits = Float64[]
            if sum(spectra[b][c]) == 0
                continue
            end
            ds = spectra[b][c][2:end] .- spectra[b][c][1:end-1]
            ds = MTools.rebin(ds, 4)
            up = 0
            down = 0
            for i in 1:size(ds)[1]
                if ds[i] > thr 
                    up = i
                elseif ds[i] < -thr
                    down = i
                    if down - up < 10
                        x1 = (up - 10) * 4
                        x2 = (down + 10) *4
                        pf = curve_fit(gausslin, x1:x2, spectra[b][c][x1:x2],
                                       [sum(spectra[b][c][x1:x2]), (x1+x2)/2,
                                            10, spectra[b][c][x1], 0.0],
                                        lower=[0, x1, 1, -Inf, -Inf],
                                        upper=[Inf, x2, 20, Inf, Inf])
                        push!(fits, pf.param[2])
                    end
                    up = 0
                    down = 0
                end
            end
            lf = curve_fit(lin, fits, [1173.2, 1332.5], [0.0, 1.0])
            ige = b * 8 + round(Int, c / 2, RoundDown) + 1
            cal_ge[ige] = lf.param
        end
    end
    return cal_ge
end

"""
    Test timing diff between Ge detectors, based on Co-60 lines
"""
function test_ge_time()

    cal_ge = init_cal()

    fin = open("data/needleRU_i1329_0004_0000.caendat", "r")

    config = TOML.parsefile("config/base.toml")
    header = zeros(UInt32, 12)
    read!(fin, header)

    # A - 1173, B - 1332
    tr_A = Float64[]
    tr_B = Float64[]
    dE = 2.0
    Emin = 100.0
    dt_max = 1000.0
    ref_ge = 1

    while !eof(fin) 
        hits = read_aggregate(fin, config)
        for hit in hits
            ige = hit.board * 8 + round(Int, hit.ch / 2, RoundDown) + 1
            if ige == ref_ge
                E = lin(hit.E, cal_ge[ige])
                t = hit.ts * 4 + hit.tf / 1024 
                if abs(E - 1173.2) < dE
                    push!(tr_A, t)
                elseif abs(E - 1332.5) < dE
                    push!(tr_B, t)
                end
            end
        end
    end
    seek(fin, 0)
    read!(fin, header)
    sort!(tr_A)
    sort!(tr_B)

    dts  = Dict{Int, Vector{Float64}}()
    for i in 1:16
        dts[i] = Float64[]
    end

    while !eof(fin) 
        hits = read_aggregate(fin, config)
        for hit in hits
            if hit.board > 1 || hit.ch % 2 != 0
                continue
            end
            ige = hit.board * 8 + round(Int, hit.ch / 2, RoundDown) + 1
            E = lin(hit.E, cal_ge[ige])
            t = hit.ts * 4 + hit.tf / 1024 
            if ige == ref_ge
                continue
            else
                if abs(E - 1173.2) < dE
                    ip = findfirst(x -> x > t, tr_B)
                    if isnothing(ip) || ip == 1
                        continue
                    end
                    dti = ifelse(abs(tr_B[ip] - t) < abs(tr_B[ip-1] - t), 
                                t - tr_B[ip], t - tr_B[ip-1])
                    if abs(dti) < dt_max
                        push!(dts[ige], dti)
                    end
                elseif abs(E - 1332.5) < dE
                    ip = findfirst(x -> x > t, tr_A)
                    if isnothing(ip) || ip == 1
                        continue
                    end
                    dti = ifelse(abs(tr_A[ip] - t) < abs(tr_A[ip-1] - t), 
                                t - tr_A[ip], t - tr_A[ip-1])
                    if abs(dti) < dt_max
                        push!(dts[ige], dti)
                    end
                end
            end
        end
    end
    dts
end

struct ChFormatInfo
    samples::Int64
    dp::UInt8
    ap2::UInt8
    ap1::UInt8
    ex::UInt8
    es::Bool
    e2::Bool
    et::Bool
    ee::Bool
    dt::Bool
end

ChFormatInfo() = ChFormatInfo(0, UInt8(0), UInt8(0), UInt8(0), UInt8(0),
                              false, false, false, false, false)
"""
    As read from boards (caen) 

    In diamant top is written to qshort, tf is 0
"""
struct RawHit
    board::UInt8
    ch::UInt8
    E::UInt16
    ts::UInt64
    tf::Int16
    qshort::UInt16
end


function RawHit()
    return RawHit(zero(UInt8), zero(UInt8), zero(UInt16), 
                  zero(UInt64), zero(Int16), zero(UInt16))
end


function Base.zero(::Type{RawHit})
    return RawHit(zero(UInt8), zero(UInt8), zero(UInt16), 
                  zero(UInt64), zero(Int16), zero(UInt16))
end


function Base.write(io::IO, hit::RawHit)
    write(io, hit.board)
    write(io, hit.ch)
    write(io, hit.E)
    write(io, hit.ts)
    write(io, hit.tf)
    write(io, hit.qshort)
end


function Base.read(io::IO, ::Type{RawHit})
    board = read(io, UInt8)
    ch = read(io, UInt8)
    E = read(io, UInt16)
    ts = read(io, UInt64)
    tf = read(io, Int16)
    qshort = read(io, UInt16)
    return RawHit(board, ch, E, ts, tf, qshort)
end


function Base.show(io::IO, hit::RawHit)
    print(io, Int64(hit.board), " ",
              Int64(hit.ch), " ",
              Int64(hit.E), " ",
              Int64(hit.ts), " ",
              Int64(hit.tf), " ",
              Int64(hit.qshort))
end

@enum ParticleType begin
    UNKNOWN = 0x00
    GAMMA = 0x01
    NEUTRON = 0x02
    PROTON = 0x03
    ALPHA = 0x04
end

function Base.zero(::Type{ParticleType})
    return UNKNOWN
end

"""
    Calibrated/identified Hit
    pid - particle id see `ParticleType`
"""
struct Hit
    loc::UInt8
    E::Float64
    t::Float64
    tof::Float64
    pid::ParticleType
end


function Hit()
    return Hit(zero(UInt8), 0.0, 0.0, 0.0, zero(UInt8))
end


function Base.zero(::Type{Hit})
    return Hit(zero(UInt8), 0.0, 0.0, 0.0, zero(UInt8))
end
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


struct PidPars
    g_low::Float64
    g_high::Float64
    n_low::Float64
    n_high::Float64
    a_low::Float64
    a_high::Float64
    p_low::Float64
    p_high::Float64
end


function PidPars(config::Dict{String, Any})
    PidPars(config["pid"]["g_low"],
            config["pid"]["g_high"],
            config["pid"]["n_low"],
            config["pid"]["n_high"],
            config["pid"]["a_low"],
            config["pid"]["a_high"],
            config["pid"]["p_low"],
            config["pid"]["p_high"])
end


@enum DetectorType begin
    UNDEFINED = 0x00
    GE = 0x01
    BGO = 0x02
    NEDA = 0x03
    DIAMANT = 0x04
end


function Base.zero(::Type{DetectorType})
    return UNDEFINED
end

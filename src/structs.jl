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
    last_label: last valid label
"""
struct SpectraPars
    dE::Float64
    Emax::Float64
    E2max::Float64
    dE3::Float64
    E3max::Float64
    dt::Float64
    tmax::Float64
    dpid::Float64
    pidmax::Float64
    Mmax::Int64
    last_label::Int64
end


function SpectraPars(config::Dict{String, Any})
    SpectraPars(config["spectra"]["dE"],
                config["spectra"]["Emax"],
                config["spectra"]["E2max"],
                config["spectra"]["dE3"],
                config["spectra"]["E3max"],
                config["spectra"]["dt"],
                config["spectra"]["tmax"],
                config["spectra"]["dpid"],
                config["spectra"]["pidmax"],
                config["spectra"]["Mmax"],
                maximum(parse.(Int64, keys(config["label"])))
               )
end


"""
Parametrs of event building for scan 
    * ge_low - lowest energy in Ge to consider further
    * bgo_low - lowest ch number in BGO to consider further
    * neda_low - same
    * dia_low - same
	* beam_period - beam period, used only for testing if calculated result
                    matches this value
    * ge_dt::Float64 - correlation window for Ge hits (for event building) 
    * bgo_dt::Float64 - same
    * neda_g_dt::Float64 - same for NEDA-gamma
    * neda_n_dt::Float64 - same for NEDA-neutron
    * dia_dt::Float64 - same
    * t_delay::Float64 - shift of non-NEDA hits (Ge, BGO, DIAMANT), so they can 
                        correlate with the beam pulse (NEDA signal) including 
                        the natural plus-minus spread

"""
struct EventPars
    ge_low::Float64
    bgo_low::Float64
    neda_low::Float64
    dia_low::Float64
	beam_period::Float64
    ge_dt::Float64
    bgo_dt::Float64
    neda_g_dt::Float64
    neda_n_dt::Float64
    dia_dt::Float64
    t_delay::Float64
end


function EventPars(config::Dict{String, Any})
    EventPars(config["event"]["ge_low"],
                config["event"]["bgo_low"],
                config["event"]["neda_low"],
                config["event"]["dia_low"],
                config["event"]["beam_period"],
                config["event"]["ge_dt"],
                config["event"]["bgo_dt"],
                config["event"]["neda_g_dt"],
                config["event"]["neda_n_dt"],
                config["event"]["dia_dt"],
                config["event"]["t_delay"]
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

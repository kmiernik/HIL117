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
mutable struct Header
    Mraw::UInt8
    M::UInt8
    ge::UInt8
    gammas::UInt8
    neutrons::UInt8
    protons::UInt8
    alphas::UInt8
end


function Header()
    return Header(UInt8(0), UInt8(0), UInt8(0), UInt8(0), 
                  UInt8(0), UInt8(0), UInt8(0))
end


function Base.show(io::IO, head::Header)
    print(io, "Mraw     : ", Int(head.Mraw), "\n",
              "M        : ", Int(head.M), "\n",
              "ge       : ", Int(head.ge), "\n",
              "gammas   : ", Int(head.gammas), "\n",
              "neutrons : ", Int(head.neutrons), "\n",
              "protons  : ", Int(head.protons), "\n",
              "alphas   : ", Int(head.alphas), "\n"
              ) 
end


function Base.write(io::IO, head::Header)
    write(io, head.Mraw)
    write(io, head.M)
    write(io, head.ge)
    write(io, head.gammas)
    write(io, head.neutrons)
    write(io, head.protons)
    write(io, head.alphas)
end


function Base.read(io::IO, ::Type{Header})
    Header(read(io, UInt8), 
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
    E: calibrated energy in 100 ev (0-65535 -> 0-6535.5 keV), 0.1 keV resolution
    t: time in 10 ps (0-65535 -> 0-655.35 ns), 0.01 ns resolution
    type: 0-unknown, 1-gamma, 2-neutron, 3-proton, 4-alpha
"""
struct EDFHit
    loc::UInt8
    E::UInt16
    t::UInt16
    type::UInt8
end


function Base.show(io::IO, hit::EDFHit)
    type = "?"
    if hit.type == UInt8(1)
        type = "γ"
    elseif hit.type == UInt8(2)
        type = "n"
    elseif hit.type == UInt8(3)
        type = "p"
    elseif hit.type == UInt8(4)
        type = "α"
    end
    print(io, "loc  : ", Int(hit.loc), "\n",
              "E    : ", Int(hit.E), "\n",
              "t    : ", Int(hit.t), "\n",
              "type : ", type, "\n")
              
end


function Base.write(io::IO, hit::EDFHit)
    write(io, hit.loc)
    write(io, hit.E)
    write(io, hit.t)
    write(io, hit.type)
end


function Base.read(io::IO, ::Type{EDFHit})
    Hit(read(io, UInt8), 
        read(io, UInt16),
        read(io, UInt16),
        read(io, UInt8))
end

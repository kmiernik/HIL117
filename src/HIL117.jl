module HIL117

using TOML
using MTools
using GLMakie
using LsqFit
using Infiltrator
using Printf
using Dates
using StatsBase

include("caendat.jl")
export RawHit, read_aggregate

include("diadat.jl")

include("tools.jl")
export read_raw_data

include("dev.jl")

end

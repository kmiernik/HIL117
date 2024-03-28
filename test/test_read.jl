using Plots
using TOML


function test_read(n_agg; raw=true)
    fin = open("data/needleRU_i1363_0036_0000.caendat", "r")
    config = TOML.parsefile("config/base.toml")
    header = zeros(UInt32, 12)
    read!(fin, header)
    data = Dict{Int, Dict{Int, Vector{Int64}}}()
    for board in 0:5
        chans = Dict{Int, Vector{Int64}}()
        for ch in 0:15
            chans[ch] = zeros(Int64, 16384)
        end
        data[board] = chans
    end
    if raw
        hits = Hit[]
        for i in 1:n_agg
            append!(hits, read_aggregate(fin, config))
        end
        return hits
    else
        for i in 1:n_agg
            hits = read_aggregate(fin, config)
            for hit in hits
                if 1 < hit.E <= 16384
                    data[hit.board][hit.ch][hit.E] += 1
                end
            end
        end
        return data
    end
end

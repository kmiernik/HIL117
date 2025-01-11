using StatsBase
using GLMakie


struct Hit
    board::UInt8
    ch::UInt8
    E::UInt16
    ts::UInt64
    tf::Int16
    qshort::UInt16
end


function Base.write(io::IO, hit::Hit)
    write(io, hit.board)
    write(io, hit.ch)
    write(io, hit.E)
    write(io, hit.ts)
    write(io, hit.tf)
    write(io, hit.qshort)
end

function Base.read(io::IO, ::Type{Hit})
    board = read(io, UInt8)
    ch = read(io, UInt8)
    E = read(io, UInt16)
    ts = read(io, UInt64)
    tf = read(io, Int16)
    qshort = read(io, UInt16)
    return Hit(board, ch, E, ts, tf, qshort)
end

function Base.show(io::IO, hit::Hit)
    print(io, Int64(hit.board), " ",
              Int64(hit.ch), " ",
              Int64(hit.E), " ",
              Float64(hit.ts), " ",
              Int64(hit.tf), " ",
              Int64(hit.qshort))
end

function Hit()
    Hit(UInt8(0), UInt8(0), UInt16(0), UInt64(0), Int16(0), UInt16(0))
end

function Hit(r::Bool)
    Hit(rand(UInt8), rand(UInt8), rand(UInt16), rand(UInt64), 
        rand(Int16), rand(UInt16))
end


function test(n)
    fout = open("test.bin", "w")
    for i in 1:n
        hit = Hit(true)
        write(fout, hit)
        println(hit)
    end
    close(fout)
    println()
    fin = open("test.bin", "r")
    for i in 1:n
        hit = read(fin, Hit)
        println(hit)
    end
end

function test_co60()
    fin = open("co60_test.dat", "r")
    d = Dict(i => Int64[] for i in 1:16)

    n = 100
    k = 0
    while true
        hit = read(fin, Hit)
        i = hit.board * 8 + round(Int, hit.ch / 2, RoundDown) + 1
        if eof(fin)
            break
        end
        if 0 < i <= 16 && hit.ch % 2 == 0
            push!(d[i], hit.E)
        end
    end
    close(fin)

    fig = Figure(size=(800, 800))
    for i in 1:16
        h = fit(Histogram, d[i], 0:16384)
        ax = Axis(fig[(i-1)%4+1, round(Int, (i-1)/4+1, RoundDown)])
        stairs!(ax, h.edges[1][1:end-1], h.weights, label="$i")
        axislegend(ax)
    end
    fig
end

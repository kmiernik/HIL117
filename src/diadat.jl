function cast_u16(data, bigendian=false)
    if bigendian
        return ((data[1] % UInt16) << 8) | data[2]
    else
        return reinterpret(UInt16, data)[1]
    end
end


function cast_u24(data, bigendian=false)
    if bigendian
        return (  ((data[1] % UInt32) << 16) 
                | ((data[2] % UInt32) <<  8) 
                |   data[3])
    else
        return (  ((data[3] % UInt32) << 16) 
                | ((data[2] % UInt32) <<  8) 
                |   data[1])
    end
end


function cast_u32(data, bigendian=false)
    if bigendian
        return (  ((data[1] % UInt32) << 24) 
                | ((data[2] % UInt32) << 16) 
                | ((data[3] % UInt32) << 8) 
                |   data[4] )
    else
        return reinterpret(UInt32, data)[1]
    end
end


function cast_u48(data, bigendian=false)
    if bigendian
        return (  ((data[1] % UInt64) << 40) 
                | ((data[2] % UInt64) << 32) 
                | ((data[3] % UInt64) << 24) 
                | ((data[4] % UInt64) << 16) 
                | ((data[5] % UInt64) << 8) 
                |   data[6] )
    else
        return (  ((data[6] % UInt64) << 40) 
                | ((data[5] % UInt64) << 32) 
                | ((data[4] % UInt64) << 24) 
                | ((data[3] % UInt64) << 16) 
                | ((data[2] % UInt64) << 8) 
                |   data[1] )
    end
end


"""
Full information diamanat hit

"""
struct DiaHit
    metatype::UInt8
    framesize::UInt32
    frametype::UInt16
    datasource::UInt8
    revision::UInt8
    datasize::UInt32
    eventidx::UInt32
    eventtime::UInt64
    boardid::UInt16
    crystalid::UInt16
    status1::UInt8
    status2::UInt8
    energy::UInt32
    top::UInt32
    checksum::UInt16
end


function read_diahit_full(fin)
    header = read(fin, 8)

    metatype = header[1]
    bigendian = ifelse((metatype & 0x80) == 0x00, true, false)

    framesize = cast_u24(header[2:4], bigendian)
    frametype = cast_u16(header[6:7], bigendian)
    datasource = header[5]
    revision = header[8]
    datasize = (1 << (metatype & 0x0F)) * framesize

    if frametype == 0x0016
        data = read(fin, datasize-8)

        eventidx = cast_u32(data[1:4], bigendian)
        eventtime = cast_u48(data[5:10], bigendian)

        crystal_board_id = cast_u16(data[11:12], bigendian)
        boardid = ((crystal_board_id >> 5) & 0x07ff) % UInt16
        crystalid = (crystal_board_id & 0x001f) % UInt16
        status1 = data[13]
        status2 = data[14]
        energy = cast_u32(data[15:18], bigendian)
        top = cast_u32(data[19:22], bigendian)
        checksum = cast_u16(data[23:24], bigendian)

        return DiaHit(metatype, framesize, frametype, datasource, 
                      revision, datasize,
                      eventidx, eventtime,
                      boardid, crystalid, status1, status2,
                      energy, top, checksum)
    else
        skip(fin, datasize - 8)
        return DiaHit(metatype, framesize, frametype, datasource, 
                      revision, datasize,
                      zero(UInt32), zero(UInt64),
                      zero(UInt16), zero(UInt16), zero(UInt8), zero(UInt8),
                      zero(UInt32), zero(UInt32), zero(UInt16))
    end
end



board_translation = Dict{UInt16, UInt8}(
                        UInt16(124) => UInt8(6),
                        UInt16(148) => UInt8(7),
                        UInt16(156) => UInt8(8),
                        UInt16(204) => UInt8(9)
                       )

"""
    read_diahit(fin)

    Read next hit in fin file. Only necessary diamont fields are accessed, 
    for full processing see `read_diahit_full`.
    Energy and Top fields are multiplied by 0.5 compared to the original.

    Times are shifted by the agava_ts.
    Hits before are returned as empty (zeros)
    timescale = 2.5 - Times are converted to 4 ns time stamp 
                    (from 10 ns time stamp) (t * 10 / 4 -> t * 2.5)

    A `RawHit` struct is returned (with zeros, if not a diamont data).
    with a following reinterpretation of fields
    (512 taken from diatoroot, 64 adjusted so top/energy about 0.5)
    board -> boardid
    ch -> crystalid
    E -> energy/energy_div
    t -> eventtime
    tf -> 0
    qshort -> top/top_div

    boarid translation (hardcoded)
    124 -> 6
    148 -> 7
    156 -> 8
    204 -> 9
"""
function read_diahit(fin, agava_ts; time_scale=2.5, energy_div=512, top_div=64)
    header = read(fin, 8)

    metatype = header[1]
    bigendian = ifelse((metatype & 0x80) == 0x00, true, false)

    framesize = cast_u24(header[2:4], bigendian)
    frametype = cast_u16(header[6:7], bigendian)
    datasize = (1 << (metatype & 0x0F)) * framesize

    if frametype == 0x0016
        data = read(fin, datasize-8)

        eventtime = cast_u48(data[5:10], bigendian)
        if eventtime < agava_ts
            return RawHit()
        else
            eventtime -= agava_ts
        end

        crystal_board_id = cast_u16(data[11:12], bigendian)
        boardid = ((crystal_board_id >> 5) & 0x07ff) % UInt16
        crystalid = (crystal_board_id & 0x001f) % UInt16
        energy = cast_u32(data[15:18], bigendian)
        top = cast_u32(data[19:22], bigendian)

        if (energy / energy_div > typemax(UInt16) 
            || top / top_div > typemax(UInt16))
            return RawHit()
        end

        return RawHit(board_translation[boardid], UInt8(crystalid),
                      round(UInt16, energy / energy_div), 
                      round(UInt64, eventtime * time_scale), 
                      zero(Int16),
                      round(UInt16, top / top_div))
    else
        skip(fin, datasize - 8)
        return RawHit()
    end
end


"""
    read_diahit(fin, agava_ts, buf_size; 
                    time_scale=2.5, energy_div=512, top_div=64)

    Buffered version reading buf_size hits
"""
function read_diahit(fin, agava_ts, buf_size; 
                    time_scale=2.5, energy_div=512, top_div=64)
    hits = RawHit[]
    for i in 1:buf_size
        header = read(fin, 8)

        metatype = header[1]
        bigendian = ifelse((metatype & 0x80) == 0x00, true, false)

        framesize = cast_u24(header[2:4], bigendian)
        frametype = cast_u16(header[6:7], bigendian)
        datasize = (1 << (metatype & 0x0F)) * framesize

        if frametype == 0x0016
            data = read(fin, datasize-8)

            eventtime = cast_u48(data[5:10], bigendian)
            if eventtime < agava_ts
                continue
            else
                eventtime -= agava_ts
            end

            crystal_board_id = cast_u16(data[11:12], bigendian)
            boardid = ((crystal_board_id >> 5) & 0x07ff) % UInt16
            crystalid = (crystal_board_id & 0x001f) % UInt16
            energy = cast_u32(data[15:18], bigendian)
            top = cast_u32(data[19:22], bigendian)

            if (energy / energy_div > typemax(UInt16) 
                || top / top_div > typemax(UInt16))
                continue
            end

            push!(hits, RawHit(board_translation[boardid], UInt8(crystalid),
                        round(UInt16, energy / energy_div), 
                        round(UInt64, eventtime * time_scale), 
                        zero(Int16),
                        round(UInt16, top / top_div)))
        else
            skip(fin, datasize - 8)
            continue
        end
    end
    hits
end

"""
    read_aggregate(fin)

    Read one Board Aggregate frame from input file
    Header must be read before, so BA starts right on

"""
function read_aggregate(fin, config)
    pha_boards = config["boards"]["PHA"]
    psd_boards = config["boards"]["PSD"]

    board_header = zeros(UInt32, 4)
    read!(fin, board_header)

    signature = ((board_header[1] & 0xA0000000) >> 28) % UInt8
    if signature !== UInt8(10)
        current_pos = position(fin)
        @warn "Wrong board aggregate signature $signature at $current_pos"
        # Move back file pointer, skip just one byte (header is 16 long)
        # perhaps good board header will be found in the next call
        skip(fin, -15)
        return RawHit[]
    end
    bsize = (board_header[1] & (0x0FFFFFFF)) % UInt32

    dualchmask = (board_header[2] & 0x000000FF) % UInt8
    pattern = ((board_header[2] & 0x007FFFFF00) >> 8) % UInt16
    bf = (board_header[2] >> 26) & 1 != 0 
    boardid = ((board_header[2] & 0xF8000000) >> 27) % UInt8

    bacounter = (board_header[3] & (0x007FFFFF))
    btimetag = board_header[4]

    hits = RawHit[]
    channels = UInt8[]
    for i in 0:7
        if (dualchmask >> i) & 1 != 0
            push!(channels, Int8(i))
        end
    end

    chdata = zeros(UInt32, Int64(bsize) - 4)
    read!(fin, chdata)

    pos = 1
    for channel in channels
        chastart = pos
        chasize = (chdata[pos] & (0x7FFFFFFF)) % UInt32
        fl = (chdata[pos] >> 31) & 1 != 0 

        if fl
            nsamples = (chdata[pos+1] & 0x0000FFFF) % UInt16
            dp = ((chdata[pos+1] &  0x000F0000) >> 16) % UInt8
            ap2 = ((chdata[pos+1] & 0x00300000) >> 20) % UInt8
            ap1 = ((chdata[pos+1] & 0x00C00000) >> 22) % UInt8
            ex = ((chdata[pos+1] &  0x07000000) >> 24) % UInt8
            es = (chdata[pos+1] >> 27) & 1 != 0 
            e2 = (chdata[pos+1] >> 28) & 1 != 0 
            et = (chdata[pos+1] >> 29) & 1 != 0 
            ee = (chdata[pos+1] >> 30) & 1 != 0 
            dt = (chdata[pos+1] >> 31) & 1 != 0 
            info = ChFormatInfo(Int64(nsamples) * 8, dp, ap2, ap1, ex, 
                                es, e2, et, ee, dt)
            pos += 2
        else
            pos += 1
            info = ChFormatInfo()
        end

        while pos - chastart < chasize
            trigtime = (chdata[pos] & (0x7FFFFFFF))
            ch = (chdata[pos] >> 31) & 1 != 0 

            swords =  Int64(info.samples / 2)
            samples = zeros(UInt16, info.samples)
            for i in 1:swords
                samples[i] =    (chdata[pos+i] & 0x00003FFF) % UInt16
                samples[i+1] = ((chdata[pos+i] & 0x3FFF0000) >> 16) % UInt16
            end
            pos += swords + 1

            reject = false
            if Int64(boardid) in pha_boards
                extras2 = chdata[pos]
                energy = (chdata[pos+1] & (0x00007FFF)) % UInt16
                pu = (chdata[pos+1] >> 15) & 1 != 0 
                
                if (chdata[pos+1] >> 17) & 1 != 0 
                    # roll-over time-stamp
                    # could be useful?
                end
                if (chdata[pos+1] >> 19) & 1 != 0 
                    # fake event
                    reject = true
                end
                if (chdata[pos+1] >> 20) & 1 != 0 
                    # saturated
                    reject = true
                end
                if (chdata[pos+1] >> 25) & 1 != 0 
                    # pileup
                    reject = true
                end

                tstamp = UInt64(0)
                tfine = UInt16(0)
                if info.ex == 0x0
                    tstamp = (((extras2 & 0x00000000FFFF0000) << 15) | trigtime)
                    #baseline = (extras2 & 0x0000ffff)
                elseif info.ex == 0x2
                    tstamp = (((extras2 & 0x00000000FFFF0000) << 15) | trigtime)
                    tfine = (extras2 & 0x000003FF) % UInt16
                elseif info.ex == 0x4
                    tstamp = UInt64(trigtime)
                    #losttrigcounter = ((extras2 & 0xffff0000) >> 16) % UInt16
                    #totaltrigcounter = (extras2 & 0x0000ffff) % UInt16
                elseif info.ex == 0x5
                    tstamp = (((extras2 & 0x00000000FFFF0000) << 15) | trigtime)
                    evbefore = Float64(extras & 0x0000ffff)
                    evafter = Float64(extras & 0xffff0000)
                    tfine = round(UInt16, 
                                  (8192 - evbefore) / (evafter - evbefore) * 4)
                end
                qshort = UInt16(0)
            elseif Int64(boardid) in psd_boards
                extras = chdata[pos]
                qshort = (chdata[pos+1] & (0x00007FFF)) % UInt16
                pur = (chdata[pos+1] >> 15) & 1 != 0 
                qlong = ((chdata[pos+1] & 0xFFFF0000) >> 16) % UInt16
                if pur
                    # pile-up or saturated
                    reject = true
                end

                tstamp = UInt64(0)
                tfine = UInt16(0)
                if info.ex == 0x0
                    tstamp = (((extras & 0x00000000FFFF0000) << 15) | trigtime)
                    #baseline = (extras & 0x0000ffff)
                elseif info.ex == 0x1
                    tstamp = (((extras & 0x00000000FFFF0000) << 15) | trigtime)
                elseif info.ex == 0x2
                    tstamp = (((extras & 0x00000000FFFF0000) << 15) | trigtime)
                    tfine = (extras & 0x000003FF) % UInt16
                elseif info.ex == 0x4
                    tstamp = UInt64(trigtime)
                    #losttrigcounter = ((extras & 0xffff0000) >> 16) % UInt16
                    #totaltrigcounter = (extras & 0x0000ffff) % UInt16
                elseif info.ex == 0x5
                    tstamp = UInt64(trigtime)
                    evbefore = Float64(extras & 0x0000ffff)
                    evafter = Float64(extras & 0xffff0000)
                    tfine = round(UInt16, 
                                  (8192 - evbefore) / (evafter - evbefore) * 4)
                end

                energy = qlong
            else
                current_pos = position(fin)
                @warn "Unknown board $boardid at $current_pos"
                return hits
            end
            pos += 2

            chnumber = 2 * channel + ifelse(ch, Int8(1), Int8(0))
            if !reject
                push!(hits, 
                  RawHit(boardid, Int8(chnumber), 
                         energy, tstamp, tfine, qshort))
            end
        end
    end
    return hits
end

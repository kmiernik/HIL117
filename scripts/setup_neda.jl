

"""
    neda_lims(neda; binsize=8, thr_g=80, thr_n=16)

    Find NEDA low limits for gammas/neutrons

    neda array is an output of src/dev.jl:test_walk function
"""
function neda_lims(neda; binsize=8, thr_g=80, thr_n=16)
    for loc in 34:84
        fig = Figure(size=(800, 800))
        set_theme!(pub_theme)
        ax1 = Axis(fig[1, 1:2]; ylabel="PID")
        ax2a = Axis(fig[2, 1]; xlabel="E", ylabel="PID", yscale=log10)
        ax2b = Axis(fig[2, 2]; xlabel="PID", ylabel="E")

        E = collect(1:2048) .- 0.5
        pid = collect(0.0:0.01:1.27) .+ 0.005
    
        heatmap!(ax1, E, pid, log.(neda[loc-32, :, :]), colormap=:jet1)
        yn = rebin(vec(sum(neda[loc-32, :, 50:69], dims=2)), binsize)[:, 1]
        yg = rebin(vec(sum(neda[loc-32, :, 71:95], dims=2)), binsize)[:, 1]
        dyn = yn[2:end] .- yn[1:end-1]
        dyg = yg[2:end] .- yg[1:end-1]

        zcn = -1
        zcg = -1
        if maximum(yg[1:round(Int, 80/binsize)]) > thr_g
            # minimum mode
            for i in round(Int, 80/binsize):round(Int, 2000/binsize)
                if yg[i] > thr_g && dyg[i] < 0 && dyg[i+1] > 0
                    a = (dyg[i+1] - dyg[i]) / binsize
                    b = dyg[i] - ((i-1)*binsize + (binsize-1)/2) * a
                    zcg = -b/a
                    break
                end
            end
        else
            # thr_geshold mode
            for i in 1:round(Int, 2000/binsize)
                if zcg < 0 && yg[i] > thr_g && yg[i+1] > thr_g
                    zcg = (i-1)*binsize + (binsize-1)/2
                    break
                end
            end
        end
        if maximum(yn[1:round(Int, 80/binsize)]) > thr_n
            # minimum mode
            for i in round(Int, 80/binsize):round(Int, 2000/binsize)
                if yn[i] > thr_n && dyn[i] < 0 && dyn[i+1] > 0
                    a = (dyn[i+1] - dyn[i]) / binsize
                    b = dyn[i] - ((i-1)*binsize + (binsize-1)/2) * a
                    zcn = -b/a
                    break
                end
            end
        else
            # thr_neshold mode
            for i in 1:round(Int, 2000/binsize)
                if yn[i] > thr_n && yn[i+1] > thr_n
                    zcn = (i-1)*binsize + (binsize-1)/2
                    break
                end
            end
        end
        @printf("%03d %.1f %.1f\n", loc, zcg, zcn)
        
        stairs!(ax2a, E, vec(sum(neda[loc-32, :, 50:69], dims=2)), 
                color=:steelblue)
        stairs!(ax2a, E, vec(sum(neda[loc-32, :, 71:95], dims=2)), 
                color=:orange)

        #=
        scatterlines!(ax2a, 0:binsize:2048-2*binsize, 
                      yn[2:end] .- yn[1:end-1], color=:steelblue)
        scatterlines!(ax2a, 0:binsize:2048-2*binsize,
                      yg[2:end] .- yg[1:end-1], color=:orange)
        ylims!(ax2a, -500, 500)
        hlines!(ax2a, [0], color=:black)
        =#
        scatterlines!(ax2a, 0:binsize:2048-binsize, 
                      yn, color=:steelblue)
        vlines!(ax2a, [zcn], color=:steelblue, linestyle=:dash)
        vlines!(ax2a, [zcg], color=:orange, linestyle=:dash)
        xlims!(ax2a, 0, 300)
        ylims!(ax2a, 1e1, 1e5)
        text!(ax2a, [zcg-5, zcn+15], [1e4, 1e4], 
              text=[@sprintf("%.1f", zcg), @sprintf("%.1f", zcn)], 
              color=[:orange, :steelblue],
              rotation=pi/2)

        stairs!(ax2b, pid, vec(neda[loc-32, round(Int, zcn), :]), 
                color=:steelblue)
        vlines!(ax2b, [0.5, 0.69], color=:steelblue, linestyle=:dash)
        vlines!(ax2b, [0.71, 0.95], color=:orange, linestyle=:dash)
        xlims!(ax2b, 0, 1.1)

        vlines!(ax1, [zcn], color=:steelblue, linestyle=:dash)
        hlines!(ax1, [0.5, 0.69], color=:steelblue, linestyle=:dash)
        vlines!(ax1, [zcg], color=:orange, linestyle=:dash)
        hlines!(ax1, [0.71, 0.95], color=:orange, linestyle=:dash)
        
        save("neda_lims_$(@sprintf("%03d", loc)).png", fig)
    end  
end

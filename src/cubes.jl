"""
    icube(i, j, k)

    Position of i, j, k cube in a 3D symmetrized ggg matrix
"""
function icube(i, j, k)
    if 1 <= i <= j <= k
        return Int64(k * (k + 1) * (k + 2) / 6 
                    - k * (k + 1) / 2 + 1 
                    + (i - 1) * k 
                    - (i - 1) * (i - 2) / 2 
                    + j - i)
    else
        throw(DomainError((i, j, k),
                          "Parameters must fulfill 1 <= i <= j <= k"))
    end
end


"""
    cubei(c)

    Given index `c` in a 3D symmetrized ggg matrix returns position i, j, k

"""
function cubei(c)
    a = sqrt(729 * c^2 - 3) + 27 * c
    t = round(Int64, cbrt(a) / 3^(2/3) + 1 / cbrt(3 * a) - 1, RoundUp)
    d = round(Int64, c - (t - 1) * t * (t + 1) / 6)
    r = round(Int64, - 1/2 * sqrt((2 * t + 1)^2 - 8 * d) + t - 1 /2, RoundUp) + 1
    p = round(Int64, (r - 1) * t - (r - 1) * (r - 2) / 2 + 1)
    k = d - p + r
    return r, k, t
end


"""
    ggcube(E1, E2, dE, Emax, cubes)

    Return vector of gammas in coincidence with E1 and E2
"""
function ggcube(E1, E2, dE, Emax, cubes)
    E1, E2 = sort([E1, E2])
    i = round(Int64, E1 / dE, RoundUp)
    j = round(Int64, E2 / dE, RoundUp)
    n = round(Int64, Emax / dE, RoundUp)
    r = zeros(Int64, n)
    for k in 1:i-1
        r[k] += cubes[icube(k, i, j)]
    end
    for k in i:j-1
        r[k] += cubes[icube(i, k, j)]
    end
    for k in j:n
        r[k] += cubes[icube(i, j, k)]
    end
    return r
end



"""
    gcube(E1, dE, Emax, cubes)

    Return matrix of gamma-gamma in coincidence with E1
"""
function gcube(E1, dE, Emax, cubes)
    i = round(Int64, E1 / dE, RoundUp)
    n = round(Int64, Emax / dE, RoundUp)
    r = zeros(Int64, n, n)
    for j in 1:n
        for k in 1:j
            ii, jj, kk = sort([i, j, k]) 
            c = cubes[icube(ii, jj, kk)]
            r[k, j] += c
            r[j, k] += c
        end
    end
    return r
end

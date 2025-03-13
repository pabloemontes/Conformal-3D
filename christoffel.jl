#Christoffel symbols done by Joaquin

function christoffel!(Γ::AbstractArray,
    position::AbstractVector,
    spacetime::KerrSpacetimeKerrSchildCoordinates,
    cache::KerrChristoffelCache)
    x = position[2]
    y = position[3]
    z = position[4]    
    
    M = spacetime.M
    a = spacetime.a

    a2 = a^2
    rho2 = x^2 .+ y^2 .+ z^2
    r2 = 0.5 * (rho2 - a2) .+ sqrt(0.25 * (rho2 - a2)^2 + a2 * z^2)
    r = sqrt(r2)
    r3 = r2 * r
    r4 = r2 * r2

    #Derivatives of r(x,y,z)
    dr_dx = x * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
    dr_dy = y * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
    dr_dz = z * r * (r2 + a2)^2 / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)

    #The scalar function and the null vector of the metric

    H = M * r3 / (r4 + a2 * z^2)

    l = cache.l
    l[1] = 1.0
    l[2] = (r * x + a * y) / (r2 + a2)
    l[3] = (r * y - a * x) / (r2 + a2)
    l[4] = z / r

    #Derivatives of the null vectors ( dl[a,b]=d_a l[b]). Indexes are down.

    dl = cache.dl
    fill!(dl, 0)

    dl[2, 2] = (dr_dx * (x - 2 * r * l[2]) + r) / (r2 + a2)
    dl[2, 3] = (dr_dx * (y - 2 * r * l[3]) - a) / (r2 + a2)
    dl[2, 4] = -z / r2 * dr_dx

    dl[3, 2] = (dr_dy * (x - 2 * r * l[2]) + a) / (r2 + a2)
    dl[3, 3] = (dr_dy * (y - 2 * r * l[3]) + r) / (r2 + a2)
    dl[3, 4] = -z / r2 * dr_dy

    dl[4, 2] = dr_dz * (x - 2 * r * l[2]) / (r2 + a2)
    dl[4, 3] = dr_dz * (y - 2 * r * l[3]) / (r2 + a2)
    dl[4, 4] = 1.0 / r - z / r2 * dr_dz

    #Derivatives of the scalar function H (dH[a]=d_a H). Index is down.

    dH = cache.dH
    fill!(dH, 0)

    dH[2] = -M * r2 * (r4 - 3 * a2 * z^2) / (r4 + a2 * z^2)^2 * dr_dx
    dH[3] = -M * r2 * (r4 - 3 * a2 * z^2) / (r4 + a2 * z^2)^2 * dr_dy
    dH[4] = -M * r2 * (2 * a2 * r * z + (r4 - 3 * a2 * z^2) * dr_dz) / (r4 + a2 * z^2)^2

    # Directional derivative of H in the direction of the null vector l  (l^a d_a H)
    l_dH = -M * r2 * (r4 - a2 * z^2) / (r4 + a2 * z^2)^2

    # Tensor product of the null vector derivatives with the null vector. dlablc[a,b,c]= dl[a,b]*l[c] 
    # Derivatives of the products H*la*lb:  D[a,b,c]= d_a (H*lb*lc) (The order of fors gives the order of indexes)
    # This computation is equivalent to D[a,b,c]=dH[a]*l[b]*l[c]+H*dl[a,b]*l[c]+H*dl[a,c]*l[b]

    D = cache.D
    for i in 1:4
        for j in 1:4
            for k in 1:4
                D[i, j, k] = dH[i] * l[j] * l[k] + H * dl[i, j] * l[k] + H * dl[i, k] * l[j]
            end
        end
    end

    #Christoffel symbols

    for i in 1:4
        sign = i == 1 ? -1 : 1
        for j in 1:4
            for k in 1:4
                Γ[i, j, k] = sign * (D[j, k, i] + D[k, j, i] - D[i, j, k] +
                              2 * H * l_dH * l[i] * l[j] * l[k])
            end
        end
    end

    return nothing
end

function christoffel!(Γ::AbstractArray,
    position::AbstractVector,
    spacetime::KerrSpacetimeKerrSchildCoordinates,
    cache::KerrChristoffelCache)
    x = position[2]
    y = position[3]
    z = position[4]    
    
    M = spacetime.M
    a = spacetime.a

    a2 = a^2
    rho2 = x^2 .+ y^2 .+ z^2
    r2 = 0.5 * (rho2 - a2) .+ sqrt(0.25 * (rho2 - a2)^2 + a2 * z^2)
    r = sqrt(r2)
    r3 = r2 * r
    r4 = r2 * r2

    #Derivatives of r(x,y,z)
    dr_dx = x * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
    dr_dy = y * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
    dr_dz = z * r * (r2 + a2)^2 / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)

    #The scalar function and the null vector of the metric

    H = M * r3 / (r4 + a2 * z^2)

    l = cache.l
    l[1] = 1.0
    l[2] = (r * x + a * y) / (r2 + a2)
    l[3] = (r * y - a * x) / (r2 + a2)
    l[4] = z / r

    #Derivatives of the null vectors ( dl[a,b]=d_a l[b]). Indexes are down.

    dl = cache.dl
    fill!(dl, 0)

    dl[2, 2] = (dr_dx * (x - 2 * r * l[2]) + r) / (r2 + a2)
    dl[2, 3] = (dr_dx * (y - 2 * r * l[3]) - a) / (r2 + a2)
    dl[2, 4] = -z / r2 * dr_dx

    dl[3, 2] = (dr_dy * (x - 2 * r * l[2]) + a) / (r2 + a2)
    dl[3, 3] = (dr_dy * (y - 2 * r * l[3]) + r) / (r2 + a2)
    dl[3, 4] = -z / r2 * dr_dy

    dl[4, 2] = dr_dz * (x - 2 * r * l[2]) / (r2 + a2)
    dl[4, 3] = dr_dz * (y - 2 * r * l[3]) / (r2 + a2)
    dl[4, 4] = 1.0 / r - z / r2 * dr_dz

    #Derivatives of the scalar function H (dH[a]=d_a H). Index is down.

    dH = cache.dH
    fill!(dH, 0)

    dH[2] = -M * r2 * (r4 - 3 * a2 * z^2) / (r4 + a2 * z^2)^2 * dr_dx
    dH[3] = -M * r2 * (r4 - 3 * a2 * z^2) / (r4 + a2 * z^2)^2 * dr_dy
    dH[4] = -M * r2 * (2 * a2 * r * z + (r4 - 3 * a2 * z^2) * dr_dz) / (r4 + a2 * z^2)^2

    # Directional derivative of H in the direction of the null vector l  (l^a d_a H)
    l_dH = -M * r2 * (r4 - a2 * z^2) / (r4 + a2 * z^2)^2

    # Tensor product of the null vector derivatives with the null vector. dlablc[a,b,c]= dl[a,b]*l[c] 
    # Derivatives of the products H*la*lb:  D[a,b,c]= d_a (H*lb*lc) (The order of fors gives the order of indexes)
    # This computation is equivalent to D[a,b,c]=dH[a]*l[b]*l[c]+H*dl[a,b]*l[c]+H*dl[a,c]*l[b]

    D = cache.D
    for i in 1:4
        for j in 1:4
            for k in 1:4
                D[i, j, k] = dH[i] * l[j] * l[k] + H * dl[i, j] * l[k] + H * dl[i, k] * l[j]
            end
        end
    end

    #Christoffel symbols

    for i in 1:4
        sign = i == 1 ? -1 : 1
        for j in 1:4
            for k in 1:4
                Γ[i, j, k] = sign * (D[j, k, i] + D[k, j, i] - D[i, j, k] +
                              2 * H * l_dH * l[i] * l[j] * l[k])
            end
        end
    end

    return nothing
end

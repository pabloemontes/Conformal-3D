using LinearAlgebra

function get_pos(idx_x, idx_y, idx_z, lengthpars)
    Lx, Ly, Lz, N, M, O = lengthpars
    idx_x = idx_x-3
    idx_y = idx_y-3
    idx_z = idx_z-3
    x = -Lx/2 + (idx_x-1)*Lx/(N-1)
    y = -Ly/2 + (idx_y-1)*Ly/(M-1)
    z = -Lz/2 + (idx_z-1)*Lz/(O-1)
    return x, y, z
end



function metric_creator(lengthpars, Mass, a)
    Lx, Ly, Lz, N, M, O = lengthpars

    g_ab_table = zeros(N+6,M+6, O+6, 4,4) #empty metric
    gᵃᵇ_table = zeros(N+6,M+6, O+6,4,4) #empty metric
    
    for i in 1:N+6
        for j in 1:M+6
            for k in 1:O+6
                g_ab = @view g_ab_table[i,j,k,:,:]
                gᵃᵇ = @view gᵃᵇ_table[i,j,k,:,:]
                x,y,z = get_pos(i, j, k, lengthpars)
                r2 = x^2 + y^2 + z^2
                r = sqrt(r2)
                a2 = a*a
                
                rho2 = x^2 .+ y^2 .+ z^2
                r2 = 0.5 * (rho2 - a2) + sqrt(0.25 * (rho2 - a2)^2 + a2 * z^2)
                r = sqrt(r2)

                
                r3 = r2 * r
                r4 = r2 * r2
            
                #Derivatives of r(x,y,z)
                dr_dx = x * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
                dr_dy = y * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
                dr_dz = z * r * (r2 + a2)^2 / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
            
                #The scalar function and the null vector of the metric
            
                H2 = 2 * M * r3 / (r4 + a2 * z^2)
            
                l1 = 1.0
                l2 = (r * x + a * y) / (r2 + a2)
                l3 = (r * y - a * x) / (r2 + a2)
                l4 = z / r
            
                g_ab[1, 1] = -1.0 + H2 * l1 * l1
                g_ab[1, 2] = 0.0 + H2 * l1 * l2
                g_ab[1, 3] = 0.0 + H2 * l1 * l3
                g_ab[1, 4] = 0.0 + H2 * l1 * l4
                g_ab[2, 1] = g_ab[1, 2]
                g_ab[2, 2] = 1.0 + H2 * l2 * l2
                g_ab[2, 3] = 0.0 + H2 * l2 * l3
                g_ab[2, 4] = 0.0 + H2 * l2 * l4
                g_ab[3, 1] = g_ab[1, 3]
                g_ab[3, 2] = g_ab[2, 3]
                g_ab[3, 3] = 1.0 + H2 * l3 * l3
                g_ab[3, 4] = 0.0 + H2 * l3 * l4
                g_ab[4, 1] = g_ab[1, 4]
                g_ab[4, 2] = g_ab[2, 4]
                g_ab[4, 3] = g_ab[3, 4]
                g_ab[4, 4] = 1.0 + H2 * l4 * l4
                
                gᵃᵇ[1, 1] = -1.0 - H2 * l1 * l1
                gᵃᵇ[1, 2] = 0.0 - H2 * l1 * l2
                gᵃᵇ[1, 3] = 0.0 - H2 * l1 * l3
                gᵃᵇ[1, 4] = 0.0 - H2 * l1 * l4
                gᵃᵇ[2, 1] = gᵃᵇ[1, 2]
                gᵃᵇ[2, 2] = 1.0 - H2 * l2 * l2
                gᵃᵇ[2, 3] = 0.0 - H2 * l2 * l3
                gᵃᵇ[2, 4] = 0.0 - H2 * l2 * l4
                gᵃᵇ[3, 1] = gᵃᵇ[1, 3]
                gᵃᵇ[3, 2] = gᵃᵇ[2, 3]
                gᵃᵇ[3, 3] = 1.0 - H2 * l3 * l3
                gᵃᵇ[3, 4] = 0.0 - H2 * l3 * l4
                gᵃᵇ[4, 1] = gᵃᵇ[1, 4]
                gᵃᵇ[4, 2] = gᵃᵇ[2, 4]
                gᵃᵇ[4, 3] = gᵃᵇ[3, 4]
                gᵃᵇ[4, 4] = 1.0 - H2 * l4 * l4
                g_ab .= Diagonal([-1.0,1.0,1.0,1.0])
                gᵃᵇ .= Diagonal([-1.0,1.0,1.0,1.0])
            end
        end
    end
    return g_ab_table, gᵃᵇ_table
end


function christoffel_creator(lengthpars, Mass, a, auxtensors)
    Lx, Ly, Lz, N, M, O = lengthpars;
    l, dH, dl, D = auxtensors;
    Γ_table = zeros(N+6,M+6, O+6,4,4,4); #empty  christoffels
    for i in 1:N+6
        for j in 1:M+6
            for k in 1:O+6
                Γ = @view Γ_table[i,j,k,:,:,:]
                x,y,z = get_pos(i, j, k, lengthpars)
                r2 = x^2 + y^2 + z^2
                r = sqrt(r2)
                a2 = a*a
                #=
                if r < 2*Mass
                    x = x/(r)*2*Mass
                    z = z/(r)*2*Mass
                    r = 2*Mass
                    r2 = 4*Mass*Mass
                end
                =#
                rho2 = x^2 .+ y^2 .+ z^2
                r2 = 0.5 * (rho2 - a2) + sqrt(0.25 * (rho2 - a2)^2 + a2 * z^2)
                r = sqrt(r2)
                
                
                r3 = r2 * r
                r4 = r2 * r2
            
                #Derivatives of r(x,y,z)
                dr_dx = x * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
                dr_dy = y * r3 * (r2 + a2) / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
                dr_dz = z * r * (r2 + a2)^2 / (a2 * z^2 * (2 * r2 + a2) + r4 * rho2)
            
                #The scalar function and the null vector of the metric
            
                H = M * r3 / (r4 + a2 * z^2)
            
                l1 = 1.0
                l2 = (r * x + a * y) / (r2 + a2)
                l3 = (r * y - a * x) / (r2 + a2)
                l4 = z / r
                
                fill!(dl, 0)

                dl[2, 2] = (dr_dx * (x - 2 * r * l[2]) + r) / (r2)
                dl[2, 3] = (dr_dx * (y - 2 * r * l[3])) / (r2)
                dl[2, 4] = -z / r2 * dr_dx

                dl[3, 2] = (dr_dy * (x - 2 * r * l[2])) / (r2)
                dl[3, 3] = (dr_dy * (y - 2 * r * l[3]) + r) / (r2)
                dl[3, 4] = -z / r2 * dr_dy

                dl[4, 2] = dr_dz * (x - 2 * r * l[2]) / (r2)
                dl[4, 3] = dr_dz * (y - 2 * r * l[3]) / (r2)
                dl[4, 4] = 1.0 / r - z / r2 * dr_dz

                #Derivatives of the scalar function H (dH[a]=d_a H). Index is down.

                fill!(dH, 0)

                dH[2] = -Mass / r2 * dr_dx
                dH[3] = -Mass / r2 * dr_dy
                dH[4] = -Mass / r2 * dr_dz

                # Directional derivative of H in the direction of the null vector l  (l^a d_a H)
                l_dH = -Mass / r2

                # Tensor product of the null vector derivatives with the null vector. dlablc[a,b,c]= dl[a,b]*l[c] 
                # Derivatives of the products H*la*lb:  D[a,b,c]= d_a (H*lb*lc) (The order of fors gives the order of indexes)
                # This computation is equivalent to D[a,b,c]=dH[a]*l[b]*l[c]+H*dl[a,b]*l[c]+H*dl[a,c]*l[b]

                for i in 1:4
                    for j in 1:4
                        for k in 1:4
                            D[i, j, k] = dH[i] * l[j] * l[k] + H * dl[i, j] * l[k] + H * dl[i, k] * l[j]
                        end
                    end
                end

                #Christoffel symbols

                for i in 1:4
                    sign = (i == 1 ? -1 : 1)
                    for j in 1:4
                        for k in 1:4
                            Γ[i, j, k] =  sign * (D[j, k, i] + D[k, j, i] - D[i, j, k] +
                                        2 * H * l_dH * l[i] * l[j] * l[k])
                        end
                    end
                end
            end
        end
    end
    #Γ_table .= 0.0
    return Γ_table
end

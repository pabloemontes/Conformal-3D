
"""
    Fx!(flux_and_source, convars, absvars, χ, g_ab, gᵃᵇ, auxvectors, auxmatrices)

Compute the fluxes in the x direction and store them in flux_and_source

`flux` is an 28 element vector, the first 14 correspond to the fluxes of the conservative variables, and the last 14 are spaceholders
because we are considering the abstract variables a part of the conservative variables to have an easy place to store them.

`con_abs` is a 28 element vector which contains the conservative variables in it's first 14 indices, and the abstract variables
in it's last 14 indices.

`p` are extra parameters, and should be a list with the following elements:

`χ` is simply χ_2

`g_ab` and `gᵃᵇ` are the metric with low and upper indices

`auxvectors` and `auxmatrices` are auxiliary vectors used to avoid allocating memory inside the function
`auxvectors` contains four 4-element vectors
`auxmatrices` cointains three 4x4-element matrices.
# Examples
```julia-repl
julia> Fx!(zeros(28), zeros(28), [chi, g_ab, gᵃᵇ, auxvectors, auxmatrices])

```
"""
function Fx!(flux, con_abs, idx_x, idx_y, idx_z, par_flux)

    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = par_flux
    cflux = @view flux[1:14] #actual fluxes of the conservative variables.
    aflux = @view flux[15:end] #dummy fluxes for the abstract variables, should always be zero.
    lᵃ, l_a, ηᵃ, ξᵃ = auxvectors
    ξ_ab, ξᵃᵇ, Mᵃᵇ = auxmatrices
    convars = @view con_abs[1:14]    #conservative variables
    absvars = @view con_abs[15:end]  #abstract variables
    g_ab = @view g_ab_table[idx_x, idx_y, idx_z, :, :]
    gᵃᵇ = @view gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]
    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]

    #We fill ξ_ab
    ξ_ab_fun!(ξ_ab_vec, ξ_ab)
    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    lower_index!(l_a, g_ab, lᵃ)
    mul!(ηᵃ, ξᵃᵇ, l_a)
    #Tensors
    Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)

    #Scalars
    μ = ξᵃ'ξ_a
    ν = ξᵃ'l_a
    ψ_1 = ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_2 = lᵃ'l_a
    
    #Indices of the conservative/flux variables
    #Energy-momentum
    idT00 = 1
    idT01 = 2
    idT02 = 3
    idT03 = 4
    #A
    idA000 = 5
    idA001 = 6
    idA002 = 7
    idA003 = 8
    idA011 = 9
    idA012 = 10
    idA013 = 11
    idA022 = 12
    idA023 = 13
    idA033 = 14

    flux .= 0.0;

    cflux[idT00]  = convars[idT01] #convars[idT01] #∂t T00 = -∂x T01 + ...
    cflux[idT01]  = Tab_functions22(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T01 = -∂x T11 + ...
    cflux[idT02]  = Tab_functions23(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T02 = -∂x T12 + ...
    cflux[idT03]  = Tab_functions24(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T03 = -∂x T13 + ...
    cflux[idA000] = convars[idA001] #∂t A000 = -∂x A001 + ...
    cflux[idA001] = convars[idA011] #∂t A001 = -∂x A011 + ...
    cflux[idA002] = convars[idA012] #∂t A002 = -∂x A012 + ...
    cflux[idA003] = convars[idA013] #∂t A003 = -∂x A013 + ...
    cflux[idA011] = Aabc_functions222(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A011 = -∂x A111 + ...
    cflux[idA012] = Aabc_functions223(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A012 = -∂x A112 + ...
    cflux[idA013] = Aabc_functions224(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A013 = -∂x A113 + ...
    cflux[idA022] = Aabc_functions233(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A022 = -∂x A122 + ...
    cflux[idA023] = Aabc_functions234(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A023 = -∂x A123 + ...
    cflux[idA033] = Aabc_functions244(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A033 = -∂x A133 + ...
    
    #We set the dummy fluxes of the abstract variables to zero.
    aflux .= 0.0;
    
    return flux
end


#= Source =#
function Is!(source,absvars, T, A, idx_x, idx_y, idx_z, t,par)
    #println("Inside the source!")
    C_0, C_1, C_2, gᵃᵇ_table, auxvectors, auxmatrices, Γ_table = par
    #println("first load correct")
    csource = @view source[1:14] #actual source of the conservative variables.
    asource = @view source[15:end] #dummy sources for the abstract variables, should always be zero.

    #println("everything loaded")
    lᵃ, ξᵃ = auxvectors
    ξ_ab, ξᵃᵇ, τᵃᵇ= auxmatrices
    gᵃᵇ = gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]
    

    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]

    ξ_ab_fun!(ξ_ab_vec, ξ_ab)

    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    #Scalars
    μ = ξᵃ'ξ_a
    ν = lᵃ'ξ_a

    #Tensors
    τᵃᵇ_fun!(τᵃᵇ, μ, ν, ξᵃ, lᵃ, ξᵃᵇ, gᵃᵇ)
    
    #Indices of the conservative/flux variables
    #Energy-momentum
    idT00 = 1
    idT01 = 2
    idT02 = 3
    idT03 = 4
    #A
    idA000 = 5
    idA001 = 6
    idA002 = 7
    idA003 = 8
    idA011 = 9
    idA012 = 10
    idA013 = 11
    idA022 = 12
    idA023 = 13
    idA033 = 14
    #println("calculating source...")
    csource[idT00] = 0.0;
    csource[idT01] = 0.0;
    csource[idT02] = 0.0;
    csource[idT03] = 0.0;
    csource[idA000] = (-C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[1] - μ*gᵃᵇ[1,1]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[1]+ξᵃ[1]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[1]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,1] - (ξᵃ[1]*lᵃ[1] + ξᵃ[1]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[1] + ν*gᵃᵇ[1,1]/(3*μ)))
    csource[idA001] = (-C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[2] - μ*gᵃᵇ[1,2]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[2]+ξᵃ[2]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[2]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,2] - (ξᵃ[1]*lᵃ[2] + ξᵃ[2]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[2] + ν*gᵃᵇ[1,2]/(3*μ)))
    csource[idA002] = (-C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[3] - μ*gᵃᵇ[1,3]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[3]+ξᵃ[3]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[3]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,3] - (ξᵃ[1]*lᵃ[3] + ξᵃ[3]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[3] + ν*gᵃᵇ[1,3]/(3*μ)))
    csource[idA003] = (-C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[4] - μ*gᵃᵇ[1,4]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[4]+ξᵃ[4]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,4] - (ξᵃ[1]*lᵃ[4] + ξᵃ[4]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[4] + ν*gᵃᵇ[1,4]/(3*μ)))
    csource[idA011] = (-C_0*ν/μ^6*(ξᵃ[2]*ξᵃ[2] - μ*gᵃᵇ[2,2]/4) + 1/2*C_1/μ^5*(ξᵃ[2]*lᵃ[2]+ξᵃ[2]*lᵃ[2]-2*ξᵃ[2]*ξᵃ[2]*ν/μ) - C_2/μ^4*(ξᵃᵇ[2,2] - (ξᵃ[2]*lᵃ[2] + ξᵃ[2]*lᵃ[2])/μ + 2/3*ν/μ^2*ξᵃ[2]*ξᵃ[2] + ν*gᵃᵇ[2,2]/(3*μ)))
    csource[idA012] = (-C_0*ν/μ^6*(ξᵃ[2]*ξᵃ[3] - μ*gᵃᵇ[2,3]/4) + 1/2*C_1/μ^5*(ξᵃ[2]*lᵃ[3]+ξᵃ[3]*lᵃ[2]-2*ξᵃ[2]*ξᵃ[3]*ν/μ) - C_2/μ^4*(ξᵃᵇ[2,3] - (ξᵃ[2]*lᵃ[3] + ξᵃ[3]*lᵃ[2])/μ + 2/3*ν/μ^2*ξᵃ[2]*ξᵃ[3] + ν*gᵃᵇ[2,3]/(3*μ)))
    csource[idA013] = (-C_0*ν/μ^6*(ξᵃ[2]*ξᵃ[4] - μ*gᵃᵇ[2,4]/4) + 1/2*C_1/μ^5*(ξᵃ[2]*lᵃ[4]+ξᵃ[4]*lᵃ[2]-2*ξᵃ[2]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[2,4] - (ξᵃ[2]*lᵃ[4] + ξᵃ[4]*lᵃ[2])/μ + 2/3*ν/μ^2*ξᵃ[2]*ξᵃ[4] + ν*gᵃᵇ[2,4]/(3*μ)))
    csource[idA022] = (-C_0*ν/μ^6*(ξᵃ[3]*ξᵃ[3] - μ*gᵃᵇ[3,3]/4) + 1/2*C_1/μ^5*(ξᵃ[3]*lᵃ[3]+ξᵃ[3]*lᵃ[3]-2*ξᵃ[3]*ξᵃ[3]*ν/μ) - C_2/μ^4*(ξᵃᵇ[3,3] - (ξᵃ[3]*lᵃ[3] + ξᵃ[3]*lᵃ[3])/μ + 2/3*ν/μ^2*ξᵃ[3]*ξᵃ[3] + ν*gᵃᵇ[3,3]/(3*μ)))
    csource[idA023] = (-C_0*ν/μ^6*(ξᵃ[3]*ξᵃ[4] - μ*gᵃᵇ[3,4]/4) + 1/2*C_1/μ^5*(ξᵃ[3]*lᵃ[4]+ξᵃ[4]*lᵃ[3]-2*ξᵃ[3]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[3,4] - (ξᵃ[3]*lᵃ[4] + ξᵃ[4]*lᵃ[3])/μ + 2/3*ν/μ^2*ξᵃ[3]*ξᵃ[4] + ν*gᵃᵇ[3,4]/(3*μ)))
    csource[idA033] = (-C_0*ν/μ^6*(ξᵃ[4]*ξᵃ[4] - μ*gᵃᵇ[4,4]/4) + 1/2*C_1/μ^5*(ξᵃ[4]*lᵃ[4]+ξᵃ[4]*lᵃ[4]-2*ξᵃ[4]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[4,4] - (ξᵃ[4]*lᵃ[4] + ξᵃ[4]*lᵃ[4])/μ + 2/3*ν/μ^2*ξᵃ[4]*ξᵃ[4] + ν*gᵃᵇ[4,4]/(3*μ)))
    
    #println("About to do the christoffel symbols.")
    Γ = @view Γ_table[idx_x, idx_y, idx_z, :, :, :]
    for i in 1:4 #a
        for j in 1:4 #d
            #We have i representing the spatial coordinates (it is always the second index)
            #We have j representing the summation index (it should always be the THIRD index)
            csource[idT00]  -= Γ[1,i,j]*T[i,j] + Γ[i,i,j]*T[1,j]
            csource[idT01]  -= Γ[2,i,j]*T[i,j] + Γ[i,i,j]*T[2,j]
            csource[idT02]  -= Γ[3,i,j]*T[i,j] + Γ[i,i,j]*T[3,j]
            csource[idT03]  -= Γ[4,i,j]*T[i,j] + Γ[i,i,j]*T[4,j]
            csource[idA000] -= Γ[i,i,j]*A[j,1,1] + Γ[1,i,j]*A[1,i,j] + Γ[1,i,j]*A[1,i,j]
            csource[idA001] -= Γ[i,i,j]*A[j,1,2] + Γ[1,i,j]*A[2,i,j] + Γ[2,i,j]*A[1,i,j]
            csource[idA002] -= Γ[i,i,j]*A[j,1,3] + Γ[1,i,j]*A[3,i,j] + Γ[3,i,j]*A[1,i,j]
            csource[idA003] -= Γ[i,i,j]*A[j,1,4] + Γ[1,i,j]*A[4,i,j] + Γ[4,i,j]*A[1,i,j]
            csource[idA011] -= Γ[i,i,j]*A[j,2,2] + Γ[2,i,j]*A[2,i,j] + Γ[2,i,j]*A[2,i,j]
            csource[idA012] -= Γ[i,i,j]*A[j,2,3] + Γ[2,i,j]*A[3,i,j] + Γ[3,i,j]*A[2,i,j]
            csource[idA013] -= Γ[i,i,j]*A[j,2,4] + Γ[2,i,j]*A[4,i,j] + Γ[4,i,j]*A[2,i,j]
            csource[idA022] -= Γ[i,i,j]*A[j,3,3] + Γ[3,i,j]*A[3,i,j] + Γ[3,i,j]*A[3,i,j]
            csource[idA023] -= Γ[i,i,j]*A[j,3,4] + Γ[3,i,j]*A[4,i,j] + Γ[4,i,j]*A[3,i,j]
            csource[idA033] -= Γ[i,i,j]*A[j,4,4] + Γ[4,i,j]*A[4,i,j] + Γ[4,i,j]*A[4,i,j]
        end
    end
    


    #We set the dummy fluxes of the abstract variables to zero.
    asource .= 0.0;

    return source
end


#= Source =#
function Is_full!(source,absvars, T, A, idx_x, idx_y, t,par)
    #println("Entering source")
    #println(size(par))
    χ_1, χ_2, C_0, C_1, C_2, g_ab_table, gᵃᵇ_table, auxvectors, auxmatrices, Γ_table = par
    #println("fields loaded!")
    csource = @view source[1:14] #actual source of the conservative variables.
    asource = @view source[15:end] #dummy sources for the abstract variables, should always be zero.

    lᵃ, l_a, ηᵃ, ξᵃ = auxvectors
    τᵃᵇ, ξ_ab, ξᵃᵇ, Mᵃᵇ = auxmatrices
    
    gᵃᵇ = gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]
    g_ab = g_ab_table[idx_x, idx_y, idx_z, :, :]
    
    

    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]

    ξ_ab_fun!(ξ_ab_vec, ξ_ab)

    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    #Scalars
    μ = ξᵃ'ξ_a
    ν = lᵃ'ξ_a

    #Tensors
    τᵃᵇ_fun!(τᵃᵇ, μ, ν, ξᵃ, lᵃ, ξᵃᵇ, gᵃᵇ)
    
    lower_index!(l_a, g_ab, lᵃ)
    mul!(ηᵃ, ξᵃᵇ, l_a)
    #Tensors
    Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)

    #Scalars
    ψ_1 = ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_2 = lᵃ'l_a
    
    
    Tab_fun!(T, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T01 = -∂x T11 + ...
    Aabc_fun!(A, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T01 = -∂x T11 + ...
    
    #Indices of the conservative/flux variables
    #Energy-momentum
    idT00 = 1
    idT01 = 2
    idT02 = 3
    idT03 = 4
    #A
    idA000 = 5
    idA001 = 6
    idA002 = 7
    idA003 = 8
    idA011 = 9
    idA012 = 10
    idA013 = 11
    idA022 = 12
    idA023 = 13
    idA033 = 14

    csource[idT00] = 0.0;
    csource[idT01] = 0.0;
    csource[idT02] = 0.0;
    csource[idT03] = 0.0;
    csource[idA000] = -C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[1] - μ*gᵃᵇ[1,1]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[1]+ξᵃ[1]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[1]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,1] - (ξᵃ[1]*lᵃ[1] + ξᵃ[1]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[1] + ν*gᵃᵇ[1,1]/(3*μ))
    csource[idA001] = -C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[2] - μ*gᵃᵇ[1,2]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[2]+ξᵃ[2]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[2]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,2] - (ξᵃ[1]*lᵃ[2] + ξᵃ[2]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[2] + ν*gᵃᵇ[1,2]/(3*μ))
    csource[idA002] = -C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[3] - μ*gᵃᵇ[1,3]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[3]+ξᵃ[3]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[3]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,3] - (ξᵃ[1]*lᵃ[3] + ξᵃ[3]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[3] + ν*gᵃᵇ[1,3]/(3*μ))
    csource[idA003] = -C_0*ν/μ^6*(ξᵃ[1]*ξᵃ[4] - μ*gᵃᵇ[1,4]/4) + 1/2*C_1/μ^5*(ξᵃ[1]*lᵃ[4]+ξᵃ[4]*lᵃ[1]-2*ξᵃ[1]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[1,4] - (ξᵃ[1]*lᵃ[4] + ξᵃ[4]*lᵃ[1])/μ + 2/3*ν/μ^2*ξᵃ[1]*ξᵃ[4] + ν*gᵃᵇ[1,4]/(3*μ))
    csource[idA011] = -C_0*ν/μ^6*(ξᵃ[2]*ξᵃ[2] - μ*gᵃᵇ[2,2]/4) + 1/2*C_1/μ^5*(ξᵃ[2]*lᵃ[2]+ξᵃ[2]*lᵃ[2]-2*ξᵃ[2]*ξᵃ[2]*ν/μ) - C_2/μ^4*(ξᵃᵇ[2,2] - (ξᵃ[2]*lᵃ[2] + ξᵃ[2]*lᵃ[2])/μ + 2/3*ν/μ^2*ξᵃ[2]*ξᵃ[2] + ν*gᵃᵇ[2,2]/(3*μ))
    csource[idA012] = -C_0*ν/μ^6*(ξᵃ[2]*ξᵃ[3] - μ*gᵃᵇ[2,3]/4) + 1/2*C_1/μ^5*(ξᵃ[2]*lᵃ[3]+ξᵃ[3]*lᵃ[2]-2*ξᵃ[2]*ξᵃ[3]*ν/μ) - C_2/μ^4*(ξᵃᵇ[2,3] - (ξᵃ[2]*lᵃ[3] + ξᵃ[3]*lᵃ[2])/μ + 2/3*ν/μ^2*ξᵃ[2]*ξᵃ[3] + ν*gᵃᵇ[2,3]/(3*μ))
    csource[idA013] = -C_0*ν/μ^6*(ξᵃ[2]*ξᵃ[4] - μ*gᵃᵇ[2,4]/4) + 1/2*C_1/μ^5*(ξᵃ[2]*lᵃ[4]+ξᵃ[4]*lᵃ[2]-2*ξᵃ[2]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[2,4] - (ξᵃ[2]*lᵃ[4] + ξᵃ[4]*lᵃ[2])/μ + 2/3*ν/μ^2*ξᵃ[2]*ξᵃ[4] + ν*gᵃᵇ[2,4]/(3*μ))
    csource[idA022] = -C_0*ν/μ^6*(ξᵃ[3]*ξᵃ[3] - μ*gᵃᵇ[3,3]/4) + 1/2*C_1/μ^5*(ξᵃ[3]*lᵃ[3]+ξᵃ[3]*lᵃ[3]-2*ξᵃ[3]*ξᵃ[3]*ν/μ) - C_2/μ^4*(ξᵃᵇ[3,3] - (ξᵃ[3]*lᵃ[3] + ξᵃ[3]*lᵃ[3])/μ + 2/3*ν/μ^2*ξᵃ[3]*ξᵃ[3] + ν*gᵃᵇ[3,3]/(3*μ))
    csource[idA023] = -C_0*ν/μ^6*(ξᵃ[3]*ξᵃ[4] - μ*gᵃᵇ[3,4]/4) + 1/2*C_1/μ^5*(ξᵃ[3]*lᵃ[4]+ξᵃ[4]*lᵃ[3]-2*ξᵃ[3]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[3,4] - (ξᵃ[3]*lᵃ[4] + ξᵃ[4]*lᵃ[3])/μ + 2/3*ν/μ^2*ξᵃ[3]*ξᵃ[4] + ν*gᵃᵇ[3,4]/(3*μ))
    csource[idA033] = -C_0*ν/μ^6*(ξᵃ[4]*ξᵃ[4] - μ*gᵃᵇ[4,4]/4) + 1/2*C_1/μ^5*(ξᵃ[4]*lᵃ[4]+ξᵃ[4]*lᵃ[4]-2*ξᵃ[4]*ξᵃ[4]*ν/μ) - C_2/μ^4*(ξᵃᵇ[4,4] - (ξᵃ[4]*lᵃ[4] + ξᵃ[4]*lᵃ[4])/μ + 2/3*ν/μ^2*ξᵃ[4]*ξᵃ[4] + ν*gᵃᵇ[4,4]/(3*μ))
    
    
    Γ = @view Γ_table[idx_x, idx_y, idx_z, :, :, :]
    
    for i in 1:4 #a
        for j in 1:4 #d
            csource[idT00]  -= Γ[1,i,j]*T[i,j] + Γ[i,i,j]*T[1,j]
            csource[idT01]  -= Γ[2,i,j]*T[i,j] + Γ[i,i,j]*T[2,j]
            csource[idT02]  -= Γ[3,i,j]*T[i,j] + Γ[i,i,j]*T[3,j]
            csource[idT03]  -= Γ[4,i,j]*T[i,j] + Γ[i,i,j]*T[4,j]
            csource[idA000] -= Γ[i,i,j]*A[j,1,1] + Γ[1,i,j]*A[1,i,j] + Γ[1,i,j]*A[1,i,j]
            csource[idA001] -= Γ[i,i,j]*A[j,1,2] + Γ[1,i,j]*A[2,i,j] + Γ[2,i,j]*A[1,i,j]
            csource[idA002] -= Γ[i,i,j]*A[j,1,3] + Γ[1,i,j]*A[3,i,j] + Γ[3,i,j]*A[1,i,j]
            csource[idA003] -= Γ[i,i,j]*A[j,1,4] + Γ[1,i,j]*A[4,i,j] + Γ[4,i,j]*A[1,i,j]
            csource[idA011] -= Γ[i,i,j]*A[j,2,2] + Γ[2,i,j]*A[2,i,j] + Γ[2,i,j]*A[2,i,j]
            csource[idA012] -= Γ[i,i,j]*A[j,2,3] + Γ[2,i,j]*A[3,i,j] + Γ[3,i,j]*A[2,i,j]
            csource[idA013] -= Γ[i,i,j]*A[j,2,4] + Γ[2,i,j]*A[4,i,j] + Γ[4,i,j]*A[2,i,j]
            csource[idA022] -= Γ[i,i,j]*A[j,3,3] + Γ[3,i,j]*A[3,i,j] + Γ[3,i,j]*A[3,i,j]
            csource[idA023] -= Γ[i,i,j]*A[j,3,4] + Γ[3,i,j]*A[4,i,j] + Γ[4,i,j]*A[3,i,j]
            csource[idA033] -= Γ[i,i,j]*A[j,4,4] + Γ[4,i,j]*A[4,i,j] + Γ[4,i,j]*A[4,i,j]
        end
    end
    


    #We set the dummy fluxes of the abstract variables to zero.
    asource .= 0.0;

    return source
end


"""
    Fx!(flux_and_source, convars, absvars, χ, g_ab, gᵃᵇ, auxvectors, auxmatrices)

Compute the fluxes in the x direction and store them in flux_and_source

`flux` is an 28 element vector, the first 14 correspond to the fluxes of the conservative variables, and the last 14 are spaceholders
because we are considering the abstract variables a part of the conservative variables to have an easy place to store them.

`con_abs` is a 28 element vector which contains the conservative variables in it's first 14 indices, and the abstract variables
in it's last 14 indices.

`p` are extra parameters, and should be a list with the following elements:

`χ` is simply χ_2

`g_ab` and `gᵃᵇ` are the metric with low and upper indices

`auxvectors` and `auxmatrices` are auxiliary vectors used to avoid allocating memory inside the function
`auxvectors` contains four 4-element vectors
`auxmatrices` cointains three 4x4-element matrices.
# Examples
```julia-repl
julia> Fx!(zeros(28), zeros(28), [chi, g_ab, gᵃᵇ, auxvectors, auxmatrices])

```
"""
function Fy!(flux, con_abs, idx_x, idx_y, idx_z, par_flux)

    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = par_flux
    cflux = @view flux[1:14] #actual fluxes of the conservative variables.
    aflux = @view flux[15:end] #dummy fluxes for the abstract variables, should always be zero.
    lᵃ, l_a, ηᵃ, ξᵃ = auxvectors
    ξ_ab, ξᵃᵇ, Mᵃᵇ = auxmatrices
    convars = @view con_abs[1:14]    #conservative variables
    absvars = @view con_abs[15:end]  #abstract variables
    g_ab = @view g_ab_table[idx_x, idx_y, idx_z, :, :]
    gᵃᵇ = @view gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]
    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]

    #We fill ξ_ab
    ξ_ab_fun!(ξ_ab_vec, ξ_ab)
    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    lower_index!(l_a, g_ab, lᵃ)
    mul!(ηᵃ, ξᵃᵇ, l_a)
    #Tensors
    Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)

    #Scalars
    μ = ξᵃ'ξ_a
    ν = ξᵃ'l_a
    ψ_1 = ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_2 = lᵃ'l_a
    
    #Indices of the conservative/flux variables
    #Energy-momentum
    idT00 = 1
    idT01 = 2
    idT02 = 3
    idT03 = 4
    #A
    idA000 = 5
    idA001 = 6
    idA002 = 7
    idA003 = 8
    idA011 = 9
    idA012 = 10
    idA013 = 11
    idA022 = 12
    idA023 = 13
    idA033 = 14

    flux .= 0.0;

    cflux[idT00]  = convars[idT02] #convars[idT01] #∂t T00 = -∂x T01 + ...
    cflux[idT01]  = Tab_functions23(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T01 = -∂x T11 + ...
    cflux[idT02]  = Tab_functions33(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T02 = -∂x T12 + ...
    cflux[idT03]  = Tab_functions34(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T03 = -∂x T13 + ...
    cflux[idA000] = convars[idA002] #∂t A000 = -∂x A001 + ...
    cflux[idA001] = convars[idA012] #∂t A001 = -∂x A011 + ...
    cflux[idA002] = convars[idA022] #∂t A002 = -∂x A012 + ...
    cflux[idA003] = convars[idA023] #∂t A003 = -∂x A013 + ...
    cflux[idA011] = Aabc_functions223(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A011 = -∂x A111 + ...
    cflux[idA012] = Aabc_functions233(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A012 = -∂x A112 + ...
    cflux[idA013] = Aabc_functions234(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A013 = -∂x A113 + ...
    cflux[idA022] = Aabc_functions333(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A022 = -∂x A122 + ...
    cflux[idA023] = Aabc_functions334(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A023 = -∂x A123 + ...
    cflux[idA033] = Aabc_functions344(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A033 = -∂x A133 + ...
    
    #We set the dummy fluxes of the abstract variables to zero.
    aflux .= 0.0;
    
    return flux
end



#This is a placeholder for the z flux. It should be fairly similar to the Fx
#Flux, so we'll leave it at that for a while
function Fz!(flux, con_abs, idx_x, idx_y, idx_z, par_flux)
    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = par_flux
    cflux = @view flux[1:14] #actual fluxes of the conservative variables.
    aflux = @view flux[15:end] #dummy fluxes for the abstract variables, should always be zero.
    lᵃ, l_a, ηᵃ, ξᵃ = auxvectors
    ξ_ab, ξᵃᵇ, Mᵃᵇ = auxmatrices

    convars = @view con_abs[1:14]    #conservative variables
    absvars = @view con_abs[15:end]  #abstract variables
    g_ab = @view g_ab_table[idx_x, idx_y, idx_z, :, :]
    gᵃᵇ = @view gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]
    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]

    #We fill ξ_ab
    ξ_ab_fun!(ξ_ab_vec, ξ_ab)
    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    lower_index!(l_a, g_ab, lᵃ)
    mul!(ηᵃ, ξᵃᵇ, l_a)
    #Tensors
    Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)

    #Scalars
    μ = ξᵃ'ξ_a
    ν = ξᵃ'l_a
    ψ_1 = ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_2 = lᵃ'l_a
    
    #Indices of the conservative/flux variables
    #Energy-momentum
    idT00 = 1
    idT01 = 2
    idT02 = 3
    idT03 = 4
    #A
    idA000 = 5
    idA001 = 6
    idA002 = 7
    idA003 = 8
    idA011 = 9
    idA012 = 10
    idA013 = 11
    idA022 = 12
    idA023 = 13
    idA033 = 14

    flux .= 0.0;

    cflux[idT00]  = convars[idT03] #convars[idT01] #∂t T00 = -∂z T03 + ...
    cflux[idT01]  = Tab_functions42(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T01 = -∂z T13 + ...
    cflux[idT02]  = Tab_functions34(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T02 = -∂z T23 + ...
    cflux[idT03]  = Tab_functions44(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T03 = -∂z T33 + ...
    cflux[idA000] = convars[idA003] #∂t A000 = -∂x A003 + ...
    cflux[idA001] = convars[idA013] #∂t A001 = -∂x A013 + ...
    cflux[idA002] = convars[idA023] #∂t A002 = -∂x A023 + ...
    cflux[idA003] = convars[idA033] #∂t A003 = -∂x A033 + ...
    cflux[idA011] = Aabc_functions224(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A011 = -∂z A113 + ...
    cflux[idA012] = Aabc_functions234(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A012 = -∂z A123 + ...
    cflux[idA013] = Aabc_functions244(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A013 = -∂z A133 + ...
    cflux[idA022] = Aabc_functions334(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A022 = -∂z A223 + ...
    cflux[idA023] = Aabc_functions344(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A023 = -∂z A233 + ...
    cflux[idA033] = Aabc_functions444(ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t A033 = -∂z A333 + ...
    
    #We set the dummy fluxes of the abstract variables to zero.
    aflux .= 0.0;
    
    return flux
end

function Speed_max(u, par_flux)
    #  Here we compute the maximal propagation speed of the equation, for the cases of real eigenvalues is the spectral radious of the 
    #  Jacobian (when the roots have imaginary values I guess it is the maximal real part of the eigenvalues).
    return 1. #para revisar.... debiera ser una función de varias cosas... 
end


function TandA!(T, A, con_abs, idx_x, idx_y, idx_z, par_flux)
    
    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = par_flux
    lᵃ, l_a, ηᵃ, ξᵃ = auxvectors
    ξ_ab, ξᵃᵇ, Mᵃᵇ = auxmatrices

    convars = @view con_abs[1:14]    #conservative variables
    absvars = @view con_abs[15:end]  #abstract variables
    g_ab = @view g_ab_table[idx_x, idx_y, idx_z, :, :]
    gᵃᵇ = @view gᵃᵇ_table[idx_x, idx_y, idx_z, :, :]
    ξ_a = absvars[1:4]
    ξ_ab_vec = absvars[5:end]
    

    #We fill ξ_ab
    ξ_ab_fun!(ξ_ab_vec, ξ_ab)
    #Vectors
    rise_index!(ξᵃ, gᵃᵇ, ξ_a)
    ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    mul!(lᵃ, ξᵃᵇ, ξ_a)
    lower_index!(l_a, g_ab, lᵃ)
    mul!(ηᵃ, ξᵃᵇ, l_a)
    #Tensors
    Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)

    #Scalars
    μ = ξᵃ'ξ_a
    ν = ξᵃ'l_a
    ψ_1 = ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_2 = lᵃ'l_a
    

    Tab_fun!(T, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T01 = -∂x T11 + ...
    Aabc_fun!(A, ξᵃ, ξᵃᵇ, χ_1, χ_2, μ, ν, ψ_1, ψ_2, lᵃ, ηᵃ,Mᵃᵇ, gᵃᵇ, g_ab) #∂t T01 = -∂x T11 + ...
    
    return T, A;
end


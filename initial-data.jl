function create_initial_data(initial_data, sol, pars)
    lengthpars, par_flux = pars
    
    Lx, Ly, Lz, N, M, O, N_Fields = lengthpars;
    dx = Lx/N
    dy = Ly/M
    dz = Lz/O
    println("Lx = $Lx, Ly = $Ly, Lz = $Lz")

    A = 0.4
    ω = 5
    δ = 0.1
    println("Initializing data...")
    #This is loaded just in case it's needed for some initial data.
    χ_1, χ_2, auxvectors, auxmatrices, g_ab_table, gᵃᵇ_table = par_flux
    if initial_data == :chichon_3D
        for i in 1:N+6
            for j in 1:M+6
                for k in 1:O+6
                    convars = @view sol[i,j,k,1:N_Fields]
                    absvars = @view sol[i,j,k,N_Fields+1:end]
                    x,y,z = get_pos(i, j, k, lengthpars)
                    r2 = x^2+y^2+z^2
                    μ = -sqrt(6)/sqrt(δ)
                    if r2 < ω^2
                        r = sqrt(r2)
                        μ = -sqrt(6)/sqrt(δ + A * (ω-r)^4*(ω+r)^4/ω^8)    
                    end
                    #μ = -sqrt(6)/sqrt(δ + A*exp(-(r2)/ω^2))
                    absvars .= 0.0
                    convars .= 0.0
                    absvars[1] = -sqrt(-(μ))
                    abs2c!(convars, absvars, i, j, k, par_flux)
                end
            end
        end
    else
        println("Error, non-existing data type.")
    end

    
    println("closing initialization function...")
    
end

        

function lower_index!(v_a, g_ab, vᵃ)
    mul!(v_a, g_ab, vᵃ)
end

function rise_index!(vᵃ, gᵃᵇ, v_a)
    mul!(vᵃ, gᵃᵇ, v_a)
end

function ξᵃᵇ_fun!(ξᵃᵇ, gᵃᵇ, ξ_ab)
    for i in 1:D
        for j in 1:D
            ξᵃᵇ[i,j] = 0.0
            for k in 1:D
                for l in 1:D
                    ξᵃᵇ[i,j] = ξᵃᵇ[i,j] + ξ_ab[k,l]*gᵃᵇ[i,k]*gᵃᵇ[j,l]
                end
            end
        end
    end
end

#Mᵃᵇ = ξᵃᶜξᵇᵈg_cd = ξ_ecξ_fd gᶜᵈgᵃᵉgᵇᶠ
function Mᵃᵇ_fun!(Mᵃᵇ, g_ab, ξᵃᵇ)
    for i in 1:D
        for j in 1:D
            Mᵃᵇ[i,j] = 0.0
            for k in 1:D
                for l in 1:D
                    Mᵃᵇ[i,j] = Mᵃᵇ[i,j] + g_ab[k,l]*ξᵃᵇ[i,k]*ξᵃᵇ[j,l]
                end
            end
        end
    end
end

function ξ_ab_fun!(ξ_ab_vec, ξ_ab)
    k = 0
    for i in 1:4
        for j in i:4
            k = k+1
            ξ_ab[i,j] = ξ_ab_vec[k]
            ξ_ab[j,i] = ξ_ab_vec[k]
        end
    end
end


function ψ_1_function(ξᵃᵇ, ξ_ab)
    ψ_1 = 0.0
    @inline for i in 1:4
        for j in 1:4
            ψ_1 += ξᵃᵇ[i,j]*ξ_ab[i,j]
        end
    end
    return ψ_1
end

function τᵃᵇ_fun!(τᵃᵇ, μ, ν, ξᵃ, lᵃ, ξᵃᵇ, gᵃᵇ)
    
    @inline for i in 1:4
        for j in i:4
            τᵃᵇ[i,j] = ξᵃᵇ[i,j] - (ξᵃ[i]*lᵃ[j]+ξᵃ[j]*lᵃ[i])/μ + (2/3)*ξᵃ[i]*ξᵃ[j]*ν/(μ^2) + ν/(3*μ)*gᵃᵇ[i,j] 
        end
    end
end

    

#Mᵃᵇ = ξᵃᶜξᵇᵈg_cd = ξ_ecξ_fd gᶜᵈgᵃᵉgᵇᶠ
function Mᵃᵇ_fun_2!(Mᵃᵇ, gᵃᵇ, ξ_ab)
    for i in 1:4
        for j in 1:4
            Mᵃᵇ[i,j] = 0.0
            for k in 1:4
                for l in 1:4
                    for n in 1:4
                        for m in 1:4
                            Mᵃᵇ[i,j] = Mᵃᵇ[i,j] + gᵃᵇ[k,l]*gᵃᵇ[i,m]*gᵃᵇ[j,n]*ξ_ab[m,k]*ξ_ab[n,l]
                        end
                    end
                end
            end
        end
    end
end
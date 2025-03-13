function remove_trace(ξ_ab_vec, g_ab, gᵃᵇ)
    return;
    trace = 0.0
    k = 0
    for i in 1:4
        for j in i:4
            k = k+1
            trace += (ξ_ab_vec[k] * gᵃᵇ[i,j])*(i == j ? 1 : 2)
        end
    end
    k = 0
    for i in 1:4
        for j in i:4
            k = k+1
            ξ_ab_vec[k] -= 0.25*g_ab[i,j] * trace
        end
    end
end
function mink_metric!(i, j, lengthpars, g_ab, gᵃᵇ)
    g_ab .= 0.0
    g_ab[1,1] = -1.0
    g_ab[2,2] = 1.0
    g_ab[3,3] = 1.0
    g_ab[4,4] = 1.0
    gᵃᵇ .= inv(g_ab)
end
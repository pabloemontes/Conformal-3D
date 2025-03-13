@inline function mysign_zero(a)
    return (1.0.*(a .> 0.0) + (-1.0).* (a .< 0.0))
end

@inline function minmod(a, b)
    return 0.5*(mysign_zero(a)+mysign_zero(b))*min(abs(a), abs(b))
end

@inline function minmod(a, b, c)
    sgnbc = (mysign_zero(b)+mysign_zero(c)) #this is 2 if both are positive, -2 if both are negative, 0 otherwise
    sgnac = (mysign_zero(a)+mysign_zero(c)) #this is 2 if both are positive, -2 if both are negative, 0 otherwise
    
    return 0.25*sgnbc*abs(sgnac)*min(abs(a), abs(b),abs(c))
end

@inline function minmod(a,b,c,d)
    return 0.125*(mysign_zero(a)+mysign_zero(b))*(abs((mysign_zero(a)+mysign_zero(c))*(mysign_zero(a)+mysign_zero(d))))*min(abs(a),abs(b),abs(c), abs(d))
end

function MM3(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat, weight::AbstractFloat) #(2*D0,Dp,Dm)
    
    weight = weight * 2.
  
    if (abs(a) <= (weight*abs(b))) 
        return (abs(a) <= (weight*abs(c))) ? abs(a)*.5 : abs(c) 
    else 
        return (abs(b) <= abs(c)) ? abs(b) : abs(c)
    end 
end
#MM3(4.,-1.,4.,2.5)
function MM3N(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat) #(2*D0,Dp,Dm)
    if (abs(a) <= abs(b))
        return (abs(a) <= abs(c)) ? abs(a) : abs(c) 
    else 
        return (abs(b) <= abs(c)) ? abs(b) : abs(c)
    end 
end
#MM3(4.,-1.,4.,2.5)


function DMM(a::AbstractFloat, b::AbstractFloat)
    return 0.5 * (mysign_zero(a) + mysign_zero(b)) * minimum([abs(a),abs(b)])
end
#DMM(2.,2.)

function DM4(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat, d::AbstractFloat)
    return 0.125 * (mysign_zero(a) + mysign_zero(b)) * abs((mysign_zero(a) + mysign_zero(c)) * (mysign_zero(a) + mysign_zero(d))) * minimum([abs(a),abs(b),abs(c),abs(d)])
end
#DM4(1.,2.,3.,0.)

function MP5(F0::AbstractFloat, F1::AbstractFloat, F2::AbstractFloat, F3::AbstractFloat, F4::AbstractFloat)::AbstractFloat
    b1 = 0.0166666666667
    b2 = 1.3333333333333
    α = 4.
    ϵ = 1.e-10
    #remember that we are off 3 places, i-2 = 0 
    vor = b1 * (2. * F0 - 13. * F1 + 47. * F2 + 27. * F3 - 3. * F4)
                vmp = F2 + DMM(F3 - F2, α * (F2 - F1));
    if ((vor - F2) * (vor - vmp) <= ϵ) 
        return vor;
    else 
        djm1 = F0 - 2. * F1 + F2;
        dj = F1 - 2. * F2 + F3;
        djp1 = F2 - 2. * F3 + F4;
        dm4jph = DM4(4. * dj - djp1, 4 * djp1 - dj, dj, djp1);
        dm4jmh = DM4(4. * dj - djm1, 4 * djm1 - dj, dj, djm1);
        vul = F2 + α * (F2 - F1);
        vav = 0.5 * (F2 + F3);
        vmd = vav - 0.5 * dm4jph;
        vlc = F2 + 0.5 * (F2 - F1) + b2 * dm4jmh;
        vmin = maximum([minimum([F2,F3,vmd]), minimum([F2, vul, vlc])]);
        vmax = minimum([maximum([F2,F3,vmd]), maximum([F2, vul, vlc])]);
        return vor + DMM(vmin - vor, vmax - vor);
    end
end

function mp5(du::Array{Float64,1}, u::Array{Float64,1}, par, j, t) # j is the grid position
   
    h_1, U, N, par_flux, par_source, Fx, Speed_max, Source = par
    
    v = reshape(u,(N,U))
    dv = reshape(du,(N,U))
    
    F_RM = zeros(U)
    F_RP = zeros(U)
    F_LM = zeros(U)
    F_LP = zeros(U)
    H_p = zeros(U)
    H_m = zeros(U)
    
    u_l = v[j,:]
    u_l_p = v[mod1((j+1), N),:]
    u_l_m = v[mod1((j + N -1), N),:]
    u_l_p2 = v[mod1((j+2), N),:]
    u_l_m2 = v[mod1((j + N -2), N),:]
    u_l_p3 = v[mod1((j+3), N),:]
    u_l_m3 = v[mod1((j + N -3), N),:]
    
    S_MAX = maximum([Speed_max(u_l_p3, par_flux), Speed_max(u_l_m3, par_flux), 
            Speed_max(u_l_p2, par_flux), Speed_max(u_l_m2, par_flux), Speed_max(u_l_p, par_flux), 
            Speed_max(u_l_m, par_flux), Speed_max(u_l, par_flux)])
    
    F_Pp3 = 0.5 * (Fx(u_l_p3, par_flux) + S_MAX * u_l_p3)
    F_Mp3 = 0.5 * (Fx(u_l_p3, par_flux) - S_MAX * u_l_p3)
    F_Pp2 = 0.5 * (Fx(u_l_p2, par_flux) + S_MAX * u_l_p2)
    F_Mp2 = 0.5 * (Fx(u_l_p2, par_flux) - S_MAX * u_l_p2)
    F_Pp  = 0.5 * (Fx(u_l_p,  par_flux) + S_MAX * u_l_p)
    F_Mp  = 0.5 * (Fx(u_l_p,  par_flux) - S_MAX * u_l_p)
    F_P   = 0.5 * (Fx(u_l,    par_flux) + S_MAX * u_l)
    F_M   = 0.5 * (Fx(u_l,    par_flux) - S_MAX * u_l)
    F_Pm  = 0.5 * (Fx(u_l_m,  par_flux) + S_MAX * u_l_m)
    F_Mm  = 0.5 * (Fx(u_l_m,  par_flux) - S_MAX * u_l_m)
    F_Pm2 = 0.5 * (Fx(u_l_m2, par_flux) + S_MAX * u_l_m2)
    F_Mm2 = 0.5 * (Fx(u_l_m2, par_flux) - S_MAX * u_l_m2)
    F_Pm3 = 0.5 * (Fx(u_l_m3, par_flux) + S_MAX * u_l_m3)
    F_Mm3 = 0.5 * (Fx(u_l_m3, par_flux) - S_MAX * u_l_m3)
    
    for i in 1:U
        F_RM[i] = MP5(F_Mp2[i], F_Mp[i],  F_M[i],  F_Mm[i], F_Mm2[i])
        F_LM[i] = MP5(F_Pm3[i], F_Pm2[i], F_Pm[i], F_P[i],  F_Pp[i])
        F_LP[i] = MP5(F_Pm2[i], F_Pm[i],  F_P[i],  F_Pp[i], F_Pp2[i])
        F_RP[i] = MP5(F_Mp3[i], F_Mp2[i], F_Mp[i], F_M[i],  F_Mm[i])
        
        H_p[i] = F_LP[i] + F_RP[i]
        H_m[i] = F_LM[i] + F_RM[i]

        dv[j,i] = - h_1 * (H_p[i]-H_m[i]) + Source(u,t,par_source)[i]
    end
    return du[j,:]
end
#mp5(u,du,par,2)


function kt(du::Array{Float64,1}, u::Array{Float64,1}, par, j) # j is the grid position
    par_flux, h_1, U, N, Fx, Speed_max = par 
    theta = 1.5
    v = reshape(u,(N,U))
    dv = reshape(du,(N,U))
    
    v_x = zeros(U)
    v_xp = zeros(U)
    v_xm = zeros(U)
    
    u_l = v[j,:]
    u_l_p = v[mod1((j+1), N),:]
    u_l_m = v[mod1((j + N -1), N),:]
    u_l_p2 = v[mod1((j+2), N),:]
    u_l_m2 = v[mod1((j + N -2), N),:]
    
    Dp = u_l_p - u_l
    Dpp = u_l_p2 - u_l_p
    Dm = u_l - u_l_m
    Dmm = u_l_m - u_l_m2
    
    for i in 1:U
        v_x[i]  = 0.5*h_1*(mysign_zero(Dp[i])+mysign_zero(Dm[i]))*MM3N(0.5(Dp[i]+Dm[i]),Dp[i]*theta,Dm[i]*theta);
        v_xp[i] = 0.5*h_1*(mysign_zero(Dpp[i])+mysign_zero(Dp[i]))*MM3N(0.5(Dpp[i]+Dp[i]),Dpp[i]*theta,Dp[i]*theta);
        v_xm[i] = 0.5*h_1*(mysign_zero(Dm[i])+mysign_zero(Dmm[i]))*MM3N(0.5(Dm[i]+Dmm[i]),Dm[i]*theta,Dmm[i]*theta);
    end
    
    u_pp = u_l_p - 0.5 * v_xp / h_1
    u_pm = u_l   + 0.5 * v_x / h_1
    u_mp = u_l   - 0.5 * v_x / h_1
    u_mm = u_l_m + 0.5 * v_xm / h_1
    
    a_p = maximum([Speed_max(u_pp, par_flux), Speed_max(u_pm, par_flux)])
    a_m = maximum([Speed_max(u_mm, par_flux), Speed_max(u_mp, par_flux)])
        
    H_p = 0.5 * (Fx(u_pp, par_flux) + Fx(u_pm, par_flux)) - 0.5 * a_p * (u_pp - u_pm)
    H_m = 0.5 * (Fx(u_mm, par_flux) + Fx(u_mp, par_flux)) - 0.5 * a_m * (u_mp - u_mm)

    return du[j,:]  = - h_1 * (H_p[:]-H_m[:]) + Source(u,t,par_source)[i]
end
#mp5(u,du,par,2)


#MP5 Reconstruction
function MP5reconstruction(Vjmm, Vjm, Vj, Vjp, Vjpp)
    B1 = 0.0166666666666666667  #1/60
    B2 = 1.3333333333333333333  #4/3
    eps = 1e-10
    ALPHA = 4.0
    #=Vjmm = V[1]
    Vjm = V[2]
    Vj = V[3]
    Vjp = V[4]
    Vjpp = V[5]=#
    Vor = B1*(2.0*Vjmm - 13.0*Vjm + 47.0*Vj + 27*Vjp - 3.0*Vjpp) #=This is the original interpolation.
                                                                   All that follows is the application of 
                                                                   limiters to treat shocks=#
    Vmp = Vj + minmod(Vjp-Vj, ALPHA*(Vj-Vjm))  #mp = monotonicity preserving. It's the median between v_j, v_(j+1)
                                              #and an upper limit v^UL = v_j+ALPHA(v_j-v_(j-1))
    if ((Vor-Vj)*(Vor-Vmp)) < eps             #this condition is equivalent to asking vl in [vj, v^{MP}]
        Vl = Vor #vl = v^{L}_{j+1/2}
    else
        djm1 = Vjmm - 2.0*Vjm + Vj
        dj = Vjm - 2*Vj + Vjp
        djp1 = Vj - 2.0*Vjp + Vjpp
        dm4jph = minmod(4*dj - djp1, 4*djp1-dj, dj, djp1)  #ph = plus half (+1/2)
        dm4jmh = minmod(4*dj - djm1, 4*djm1-dj, dj, djm1)  #mh = minus half (-1/2)
        #d^{M4}_{j+1/2} = \minmod(4d_{j}-d_{j+1},4d_{j+1}-d_{j}, d_{j}, d_{j+1})
        Vul = Vj + ALPHA*(Vj - Vjm)   #upper limit
        Vav = 0.5*(Vj + Vjp)          #average
        Vmd = Vav - 0.5*dm4jph        #Vmedian
        Vlc = Vj + 0.5*(Vj-Vjm) + B2*dm4jmh
        Vmin = max(min(Vj, Vjp, Vmd), min(Vj, Vul, Vlc));
        Vmax = min(max(Vj, Vjp, Vmd), max(Vj, Vul, Vlc));
        Vl = Vor + minmod(Vmin-Vor, Vmax-Vor) #this places Vor between Vmin and Vmax
    end
    return Vl
end 


function MP5reconstruction!(Vl, Vjmm, Vjm, Vj, Vjp, Vjpp, N_Fields)
    B1 = 0.0166666666666666667  #1/60
    B2 = 1.3333333333333333333  #4/3
    eps = 1e-10
    ALPHA = 4.0
    #=Vjmm = V[1]
    Vjm = V[2]
    Vj = V[3]
    Vjp = V[4]
    Vjpp = V[5]=#
    for i in 1:N_Fields
        Vor = B1*(2.0*Vjmm[i] - 13.0*Vjm[i] + 47.0*Vj[i] + 27*Vjp[i] - 3.0*Vjpp[i]) #=This is the original interpolation.
                                                                       All that follows is the application of 
                                                                       limiters to treat shocks=#
        Vmp = Vj[i] + minmod(Vjp[i]-Vj[i], ALPHA*(Vj[i]-Vjm[i]))  #mp = monotonicity preserving. It's the median between v_j, v_(j+1)
                                                  #and an upper limit v^UL = v_j+ALPHA(v_j-v_(j-1))
        if ((Vor-Vj[i])*(Vor-Vmp)) < eps             #this condition is equivalent to asking vl in [vj, v^{MP}]
            Vl[i] = Vor #vl = v^{L}_{j+1/2}
        else
            djm1 = Vjmm[i] - 2.0*Vjm[i] + Vj[i]
            dj = Vjm[i] - 2*Vj[i] + Vjp[i]
            djp1 = Vj[i] - 2.0*Vjp[i] + Vjpp[i]
            dm4jph = minmod(4*dj - djp1, 4*djp1-dj, dj, djp1)  #ph = plus half (+1/2)
            dm4jmh = minmod(4*dj - djm1, 4*djm1-dj, dj, djm1)  #mh = minus half (-1/2)
            #d^{M4}_{j+1/2} = \minmod(4d_{j}-d_{j+1},4d_{j+1}-d_{j}, d_{j}, d_{j+1})
            Vul = Vj[i] + ALPHA*(Vj[i] - Vjm[i])   #upper limit
            Vav = 0.5*(Vj[i] + Vjp[i])          #average
            Vmd = Vav - 0.5*dm4jph        #Vmedian
            Vlc = Vj[i] + 0.5*(Vj[i]-Vjm[i]) + B2*dm4jmh
            Vmin = max(min(Vj[i], Vjp[i], Vmd), min(Vj[i], Vul, Vlc));
            Vmax = min(max(Vj[i], Vjp[i], Vmd), max(Vj[i], Vul, Vlc));
            Vl[i] = Vor + minmod(Vmin-Vor, Vmax-Vor) #this places Vor between Vmin and Vmax
        end
    end
end 



function mp5!(dv, v, par, t) # j is the grid position
    #asumimos u unidimensional por ahora
    h, N_Fields, N, par_flux, par_source, Fx!, Speed_max, Source!, par_mem = par
    #h_1, U, N, χ, par_flux, Flux, Speed_max, Source, par_mem = par_mp5 
    F_Mm3, F_Mm2, F_Mm1, F_M, F_Mp1, F_Mp2, F_Mp3, F_Pm3, F_Pm2, F_Pm1, F_P, F_Pp1, F_Pp2, F_Pp3, F_LP,   F_LM, F_RP, F_RM, H_m, H_p = par_mem
    
    fields = reshape(v,(N,N_Fields))
    dfields = reshape(dv,(N,N_Fields))
    #nota: f minuscula o u se usa para hablar de campos, F mayúscula para hablar de Flujos.
    
    for idx in 1:N
        #first we defined shifted indices
        idxm3 = mod(((idx-3) - 1),N) + 1
        idxm2 = mod(((idx-2) - 1),N) + 1
        idxm1 = mod(((idx-1) - 1),N) + 1
        idxp1 = mod(((idx+1) - 1),N) + 1
        idxp2 = mod(((idx+2) - 1),N) + 1
        idxp3 = mod(((idx+3) - 1),N) + 1
        
    
        um3 = @view fields[idxm3,:]
        um2 = @view fields[idxm2,:]
        um1 = @view fields[idxm1,:]
        u   = @view fields[idx,:]
        up1 = @view fields[idxp1,:]
        up2 = @view fields[idxp2,:]
        up3 = @view fields[idxp3,:]
        
        S_MAX = max(Speed_max(up3, par_flux), Speed_max(um3, par_flux), 
            Speed_max(up2, par_flux), Speed_max(um2, par_flux), Speed_max(up1, par_flux), 
            Speed_max(um1, par_flux), Speed_max(u, par_flux)) #maximum speed
        
        Fx!(F_Pm3, um3, par_flux)
        Fx!(F_Pm2, um2, par_flux)
        Fx!(F_Pm1, um1, par_flux)
        Fx!(F_P, u, par_flux)
        Fx!(F_Pp1, up1, par_flux)
        Fx!(F_Pp2, up2, par_flux)
        Fx!(F_Pp3, up3, par_flux)
        
        
        @. F_Mm3 = 0.5 * (F_Pm3 - S_MAX * um3)
        @. F_Mm2 = 0.5 * (F_Pm2 - S_MAX * um2)
        @. F_Mm1 = 0.5 * (F_Pm1 - S_MAX * um1)
        @. F_M   = 0.5 * (F_P   - S_MAX * u)
        @. F_Mp1 = 0.5 * (F_Pp1 - S_MAX * up1)
        @. F_Mp2 = 0.5 * (F_Pp2 - S_MAX * up2)
        @. F_Mp3 = 0.5 * (F_Pp3 - S_MAX * up3)
        @. F_Pm3 = 0.5 * (F_Pm3 + S_MAX * um3)
        @. F_Pm2 = 0.5 * (F_Pm2 + S_MAX * um2)
        @. F_Pm1 = 0.5 * (F_Pm1 + S_MAX * um1)
        @. F_P   = 0.5 * (F_P   + S_MAX * u)
        @. F_Pp1 = 0.5 * (F_Pp1 + S_MAX * up1)
        @. F_Pp2 = 0.5 * (F_Pp2 + S_MAX * up2)
        @. F_Pp3 = 0.5 * (F_Pp3 + S_MAX * up3)
    
        MP5reconstruction!(F_RM, F_Mp2, F_Mp1,  F_M,  F_Mm1, F_Mm2, N_Fields)
        MP5reconstruction!(F_LM, F_Pm3, F_Pm2, F_Pm1, F_P,  F_Pp1, N_Fields)
        MP5reconstruction!(F_LP, F_Pm2, F_Pm1,  F_P,  F_Pp1, F_Pp2, N_Fields)
        MP5reconstruction!(F_RP, F_Mp3, F_Mp2, F_Mp1, F_M,  F_Mm1, N_Fields)
        
        @. H_p = F_LP + F_RP
        @. H_m = F_LM + F_RM
        
        Source!(sourcevec, u, t, par_source)
        @. dfields[idx, :] = -h*(H_p - H_m) + sourcevec #+ Source(u,t,par_source)
        
    end
    
end


#====================WENOZ====================#

function createWENOZvectors(N_FIELDS)
    F_Mm3 = Array{Float64}(undef, N_FIELDS)
    F_Mm2 = copy(F_Mm3)
    F_Mm1 = copy(F_Mm3)
    F_M   = copy(F_Mm3)
    F_Mp1 = copy(F_Mm3)
    F_Mp2 = copy(F_Mm3)
    F_Mp3 = copy(F_Mm3)
    F_Pm3 = copy(F_Mm3)
    F_Pm2 = copy(F_Mm3)
    F_Pm1 = copy(F_Mm3)
    F_P   = copy(F_Mm3)
    F_Pp1 = copy(F_Mm3)
    F_Pp2 = copy(F_Mm3)
    F_Pp3 = copy(F_Mm3)
    F_LP  = copy(F_Mm3)
    F_LM  = copy(F_Mm3)
    F_RP  = copy(F_Mm3)
    F_RM  = copy(F_Mm3)
    H_m   = copy(F_Mm3)
    H_p   = copy(F_Mm3)
    sourcevec = copy(F_Mm3)
    
    return (F_Mm3, F_Mm2, F_Mm1, F_M, F_Mp1, F_Mp2, F_Mp3, F_Pm3, F_Pm2, F_Pm1, F_P, F_Pp1, F_Pp2, F_Pp3, F_LP, F_LM, F_RP, F_RM, H_m, H_p, sourcevec)
end

function WENOZreconstruction!(Vl, Vjmm, Vjm, Vj, Vjp, Vjpp, N_Fields)
    B1 = 1.0833333333333333333  #13/12
    B2 = 0.1666666666666666666  #1/6
    
    eps = 1e-40

    for i in 1:N_Fields
        Q0 = 2.0*Vj[i] +5.0*Vjp[i] - 1.0*Vjpp[i]
        Q1 = -Vjm[i] + 5.0*Vj[i] + 2.0*Vjp[i]
        Q2 = 2.0*Vjmm[i] - 7.0*Vjm[i] + 11* Vj[i]


        β0 =  B1*(Vj[i] - 2*Vjp[i] + Vjpp[i])^2 + 0.25*(3*Vj[i] - 4*Vjp[i]+ Vjpp[i])^2
        β1 =  B1*(Vjm[i] - 2*Vj[i] + Vjp[i])^2 + 0.25*(Vjm[i] - Vjp[i])^2
        β2 =  B1*(Vjmm[i] - 2*Vjm[i] + Vj[i])^2 + 0.25*(Vjmm[i] - 4*Vjm[i]+ 3*Vj[i])^2


        τ5 = abs(β2 - β0)

        α0 = 0.3*(1.0 + (τ5/(β0 + eps))^2)
        α1 = 0.6*(1.0 + (τ5/(β1 + eps))^2)
        α2 = 0.1*(1.0 + (τ5/(β2 + eps))^2)

        alphasum = (α0 + α1 + α2)

        Vl[i] = (α0*Q0 + α1*Q1 + α2*Q2)*B2/alphasum
    end
end

function wenoz!(dfields, fields, par, t) # j is the grid position
    #asumimos u unidimensional por ahora
    lengthpars, par_eq, par_source, Fx!, Fy!, Fz!, source!, Speed_max, auxvecs, T, A = par
    Lx, Ly, Lz, N, M, O, N_Fields = lengthpars;
    hx = N/Lx;
    hy = M/Ly;
    hz = O/Lz;
    F_Mm3, F_Mm2, F_Mm1, F_M, F_Mp1, F_Mp2, F_Mp3, F_Pm3, F_Pm2, F_Pm1, F_P, F_Pp1, F_Pp2, F_Pp3, F_LP, F_LM, F_RP, F_RM, H_m, H_p, sourcevec = auxvecs;
    #nota: f minuscula o u se usa para hablar de campos, F mayúscula para hablar de Flujos.
    for idx_x in 4:N+3
        for idx_y in 4:N+3
            for idx_z in 4:O+3
                #first we define shifted indices
                idx_xm3 = mod(((idx_x-3) - 1),N) + 1
                idx_xm2 = mod(((idx_x-2) - 1),N) + 1
                idx_xm1 = mod(((idx_x-1) - 1),N) + 1
                idx_xp1 = mod(((idx_x+1) - 1),N) + 1
                idx_xp2 = mod(((idx_x+2) - 1),N) + 1
                idx_xp3 = mod(((idx_x+3) - 1),N) + 1
                
                idx_ym3 = mod(((idx_y-3) - 1),M) + 1
                idx_ym2 = mod(((idx_y-2) - 1),M) + 1
                idx_ym1 = mod(((idx_y-1) - 1),M) + 1
                idx_yp1 = mod(((idx_y+1) - 1),M) + 1
                idx_yp2 = mod(((idx_y+2) - 1),M) + 1
                idx_yp3 = mod(((idx_y+3) - 1),M) + 1
                
                idx_zm3 = mod(((idx_z-3) - 1),O) + 1
                idx_zm2 = mod(((idx_z-2) - 1),O) + 1
                idx_zm1 = mod(((idx_z-1) - 1),O) + 1
                idx_zp1 = mod(((idx_z+1) - 1),O) + 1
                idx_zp2 = mod(((idx_z+2) - 1),O) + 1
                idx_zp3 = mod(((idx_z+3) - 1),O) + 1
            
                u_xm3 = @view fields[idx_xm3, idx_y, idx_z,:]
                u_xm2 = @view fields[idx_xm2, idx_y, idx_z,:]
                u_xm1 = @view fields[idx_xm1, idx_y, idx_z,:]
                u_x   = @view fields[idx_x  , idx_y, idx_z,:]
                u_xp1 = @view fields[idx_xp1, idx_y, idx_z,:]
                u_xp2 = @view fields[idx_xp2, idx_y, idx_z,:]
                u_xp3 = @view fields[idx_xp3, idx_y, idx_z,:]
                
                u_ym3 = @view fields[idx_x, idx_ym3, idx_z,:]
                u_ym2 = @view fields[idx_x, idx_ym2, idx_z,:]
                u_ym1 = @view fields[idx_x, idx_ym1, idx_z,:]
                u_y   = @view fields[idx_x, idx_y  , idx_z,:]
                u_yp1 = @view fields[idx_x, idx_yp1, idx_z,:]
                u_yp2 = @view fields[idx_x, idx_yp2, idx_z,:]
                u_yp3 = @view fields[idx_x, idx_yp3, idx_z,:]

                u_zm3 = @view fields[idx_x, idx_y, idx_zm3, :]
                u_zm2 = @view fields[idx_x, idx_y, idx_zm2, :]
                u_zm1 = @view fields[idx_x, idx_y, idx_zm1, :]
                u_z   = @view fields[idx_x, idx_y, idx_z  , :]
                u_zp1 = @view fields[idx_x, idx_y, idx_zp1, :]
                u_zp2 = @view fields[idx_x, idx_y, idx_zp2, :]
                u_zp3 = @view fields[idx_x, idx_y, idx_zp3, :]
                

                S_MAX_X = max(Speed_max(u_xp3, par_eq), Speed_max(u_xm3, par_eq), 
                    Speed_max(u_xp2, par_eq), Speed_max(u_xm2, par_eq), Speed_max(u_xp1, par_eq), 
                    Speed_max(u_xm1, par_eq), Speed_max(u_x, par_eq)) #maximum speed
                
                S_MAX_Y = max(Speed_max(u_yp3, par_eq), Speed_max(u_ym3, par_eq), 
                    Speed_max(u_yp2, par_eq), Speed_max(u_ym2, par_eq), Speed_max(u_yp1, par_eq), 
                    Speed_max(u_ym1, par_eq), Speed_max(u_y, par_eq)) #maximum speed
                
                S_MAX_Z = max(Speed_max(u_zp3, par_eq), Speed_max(u_zm3, par_eq), 
                    Speed_max(u_zp2, par_eq), Speed_max(u_zm2, par_eq), Speed_max(u_zp1, par_eq), 
                    Speed_max(u_zm1, par_eq), Speed_max(u_z, par_eq)) #maximum speed
                
                
                Fx!(F_Pm3, u_xm3, idx_xm3, idx_y, idx_z, par_eq)
                Fx!(F_Pm2, u_xm2, idx_xm2, idx_y, idx_z, par_eq)
                Fx!(F_Pm1, u_xm1, idx_xm1, idx_y, idx_z, par_eq)
                Fx!(F_P  , u_x  , idx_x  , idx_y, idx_z, par_eq)
                Fx!(F_Pp1, u_xp1, idx_xp1, idx_y, idx_z, par_eq)
                Fx!(F_Pp2, u_xp2, idx_xp2, idx_y, idx_z, par_eq)
                Fx!(F_Pp3, u_xp3, idx_xp3, idx_y, idx_z, par_eq)
    
    
                @. F_Mm3 = 0.5 * (F_Pm3 - S_MAX_X * u_xm3)
                @. F_Mm2 = 0.5 * (F_Pm2 - S_MAX_X * u_xm2)
                @. F_Mm1 = 0.5 * (F_Pm1 - S_MAX_X * u_xm1)
                @. F_M   = 0.5 * (F_P   - S_MAX_X * u_x)
                @. F_Mp1 = 0.5 * (F_Pp1 - S_MAX_X * u_xp1)
                @. F_Mp2 = 0.5 * (F_Pp2 - S_MAX_X * u_xp2)
                @. F_Mp3 = 0.5 * (F_Pp3 - S_MAX_X * u_xp3)
                @. F_Pm3 = 0.5 * (F_Pm3 + S_MAX_X * u_xm3)
                @. F_Pm2 = 0.5 * (F_Pm2 + S_MAX_X * u_xm2)
                @. F_Pm1 = 0.5 * (F_Pm1 + S_MAX_X * u_xm1)
                @. F_P   = 0.5 * (F_P   + S_MAX_X * u_x)
                @. F_Pp1 = 0.5 * (F_Pp1 + S_MAX_X * u_xp1)
                @. F_Pp2 = 0.5 * (F_Pp2 + S_MAX_X * u_xp2)
                @. F_Pp3 = 0.5 * (F_Pp3 + S_MAX_X * u_xp3)
                    
                WENOZreconstruction!(F_RM, F_Mp2, F_Mp1, F_M  , F_Mm1, F_Mm2, N_Fields)
                WENOZreconstruction!(F_LM, F_Pm3, F_Pm2, F_Pm1, F_P  , F_Pp1, N_Fields)
                WENOZreconstruction!(F_LP, F_Pm2, F_Pm1, F_P  , F_Pp1, F_Pp2, N_Fields)
                WENOZreconstruction!(F_RP, F_Mp3, F_Mp2, F_Mp1, F_M  , F_Mm1, N_Fields)
                
                @. H_p = F_LP + F_RP
                @. H_m = F_LM + F_RM

                  
                @. dfields[idx_x, idx_y, idx_z, :] = (-hx)*(H_p - H_m)
                
                    
                Fy!(F_Pm3, u_ym3, idx_x, idx_ym3, idx_z, par_eq)
                Fy!(F_Pm2, u_ym2, idx_x, idx_ym2, idx_z, par_eq)
                Fy!(F_Pm1, u_ym1, idx_x, idx_ym1, idx_z, par_eq)
                Fy!(F_P  , u_y  , idx_x, idx_y  , idx_z, par_eq)
                Fy!(F_Pp1, u_yp1, idx_x, idx_yp1, idx_z, par_eq)
                Fy!(F_Pp2, u_yp2, idx_x, idx_yp2, idx_z, par_eq)
                Fy!(F_Pp3, u_yp3, idx_x, idx_yp3, idx_z, par_eq)


                @. F_Mm3 = 0.5 * (F_Pm3 - S_MAX_Y * u_ym3)
                @. F_Mm2 = 0.5 * (F_Pm2 - S_MAX_Y * u_ym2)
                @. F_Mm1 = 0.5 * (F_Pm1 - S_MAX_Y * u_ym1)
                @. F_M   = 0.5 * (F_P   - S_MAX_Y * u_y)
                @. F_Mp1 = 0.5 * (F_Pp1 - S_MAX_Y * u_yp1)
                @. F_Mp2 = 0.5 * (F_Pp2 - S_MAX_Y * u_yp2)
                @. F_Mp3 = 0.5 * (F_Pp3 - S_MAX_Y * u_yp3)
                @. F_Pm3 = 0.5 * (F_Pm3 + S_MAX_Y * u_ym3)
                @. F_Pm2 = 0.5 * (F_Pm2 + S_MAX_Y * u_ym2)
                @. F_Pm1 = 0.5 * (F_Pm1 + S_MAX_Y * u_ym1)
                @. F_P   = 0.5 * (F_P   + S_MAX_Y * u_y)
                @. F_Pp1 = 0.5 * (F_Pp1 + S_MAX_Y * u_yp1)
                @. F_Pp2 = 0.5 * (F_Pp2 + S_MAX_Y * u_yp2)
                @. F_Pp3 = 0.5 * (F_Pp3 + S_MAX_Y * u_yp3)
                
                WENOZreconstruction!(F_RM, F_Mp2, F_Mp1, F_M  , F_Mm1, F_Mm2, N_Fields)
                WENOZreconstruction!(F_LM, F_Pm3, F_Pm2, F_Pm1, F_P  , F_Pp1, N_Fields)
                WENOZreconstruction!(F_LP, F_Pm2, F_Pm1, F_P  , F_Pp1, F_Pp2, N_Fields)
                WENOZreconstruction!(F_RP, F_Mp3, F_Mp2, F_Mp1, F_M  , F_Mm1, N_Fields)
                
                @. H_p = F_LP + F_RP
                @. H_m = F_LM + F_RM

                  
                @. dfields[idx_x, idx_y, idx_z, :] += (-hy)*(H_p - H_m)
                
                Fz!(F_Pm3, u_zm3, idx_x, idx_y, idx_zm3, par_eq)
                Fz!(F_Pm2, u_zm2, idx_x, idx_y, idx_zm2, par_eq)
                Fz!(F_Pm1, u_zm1, idx_x, idx_y, idx_zm1, par_eq)
                Fz!(F_P  , u_z  , idx_x, idx_y, idx_z  , par_eq)
                Fz!(F_Pp1, u_zp1, idx_x, idx_y, idx_zp1, par_eq)
                Fz!(F_Pp2, u_zp2, idx_x, idx_y, idx_zp2, par_eq)
                Fz!(F_Pp3, u_zp3, idx_x, idx_y, idx_zp3, par_eq)


                @. F_Mm3 = 0.5 * (F_Pm3 - S_MAX_Z * u_zm3)
                @. F_Mm2 = 0.5 * (F_Pm2 - S_MAX_Z * u_zm2)
                @. F_Mm1 = 0.5 * (F_Pm1 - S_MAX_Z * u_zm1)
                @. F_M   = 0.5 * (F_P   - S_MAX_Z * u_z)
                @. F_Mp1 = 0.5 * (F_Pp1 - S_MAX_Z * u_zp1)
                @. F_Mp2 = 0.5 * (F_Pp2 - S_MAX_Z * u_zp2)
                @. F_Mp3 = 0.5 * (F_Pp3 - S_MAX_Z * u_zp3)
                @. F_Pm3 = 0.5 * (F_Pm3 + S_MAX_Z * u_zm3)
                @. F_Pm2 = 0.5 * (F_Pm2 + S_MAX_Z * u_zm2)
                @. F_Pm1 = 0.5 * (F_Pm1 + S_MAX_Z * u_zm1)
                @. F_P   = 0.5 * (F_P   + S_MAX_Z * u_z)
                @. F_Pp1 = 0.5 * (F_Pp1 + S_MAX_Z * u_zp1)
                @. F_Pp2 = 0.5 * (F_Pp2 + S_MAX_Z * u_zp2)
                @. F_Pp3 = 0.5 * (F_Pp3 + S_MAX_Z * u_zp3)
                
                WENOZreconstruction!(F_RM, F_Mp2, F_Mp1, F_M  , F_Mm1, F_Mm2, N_Fields)
                WENOZreconstruction!(F_LM, F_Pm3, F_Pm2, F_Pm1, F_P  , F_Pp1, N_Fields)
                WENOZreconstruction!(F_LP, F_Pm2, F_Pm1, F_P  , F_Pp1, F_Pp2, N_Fields)
                WENOZreconstruction!(F_RP, F_Mp3, F_Mp2, F_Mp1, F_M  , F_Mm1, N_Fields)
                
                @. H_p = F_LP + F_RP
                @. H_m = F_LM + F_RM

                  
                @. dfields[idx_x, idx_y, idx_z, :] += (-hz)*(H_p - H_m)
                
                
                source!(sourcevec, u_x[15:end], T, A, idx_x, idx_y, idx_z, t, par_source) #u_x = u_y = fields[idx_x, idx_y, :]
                
                @. dfields[idx_x, idx_y, idx_z, :] += sourcevec;
            end
        end
    end
    #boundary conditions!
    for i in 1:3
        dfields[i, :, :, :].=dfields[4, :, :, :]
        dfields[N+3+i, :, :, :].=dfields[N+3, :, :, :]
        dfields[:, i, :, :].=dfields[:, :, i, 4]
        dfields[:, N+3+i, :, :].=dfields[:, N+3, :, :]
        dfields[:, :, i, :].=dfields[:, :, 4, :]
        dfields[:, :, N+3+i, :].=dfields[:, :, N+3, :]
    end
end


#=


function wenoz_axisymmetric!(dfields, fields, par, t) # j is the grid position
    
    #println("hello there")
    lengthpars, par_eq, par_source, Fx!, Fy!, source!, Speed_max, auxvecs, T, A = par
    Lx, Ly, N, M = lengthpars;
    hx = N/Lx;
    hy = M/Ly;
    Fx_Mm3, Fx_Mm2, Fx_Mm1, Fx_M, Fx_Mp1, Fx_Mp2, Fx_Mp3, Fx_Pm3, Fx_Pm2, Fx_Pm1, Fx_P, Fx_Pp1, Fx_Pp2, Fx_Pp3, Fx_LP, Fx_LM, Fx_RP, Fx_RM, Hx_m, Hx_p, Fy_Mm3, Fy_Mm2, Fy_Mm1, Fy_M, Fy_Mp1, Fy_Mp2, Fy_Mp3, Fy_Pm3, Fy_Pm2, Fy_Pm1, Fy_P, Fy_Pp1, Fy_Pp2, Fy_Pp3, Fy_LP, Fy_LM, Fy_RP, Fy_RM, Hy_m, Hy_p, sourcevec = auxvecs
    N_Fields = 14
    #nota: f minuscula o u se usa para hablar de campos, F mayúscula para hablar de Flujos.
    #println("fields loaded")
    #The first three and last three points are ghost points, and are treated differently.
    
    @inbounds for idx_x in 4:N+3
        for idx_y in 4:M+3

            #first we define shifted indices
            idx_xm3 = idx_x-3
            idx_xm2 = idx_x-2
            idx_xm1 = idx_x-1
            idx_xp1 = idx_x+1
            idx_xp2 = idx_x+2
            idx_xp3 = idx_x+3
            
            idx_ym3 = idx_y-3
            idx_ym2 = idx_y-2
            idx_ym1 = idx_y-1
            idx_yp1 = idx_y+1
            idx_yp2 = idx_y+2
            idx_yp3 = idx_y+3
            
            u_xm3 = @view fields[idx_xm3, idx_y,:]
            u_xm2 = @view fields[idx_xm2, idx_y,:]
            u_xm1 = @view fields[idx_xm1, idx_y,:]
            u_x   = @view fields[idx_x  , idx_y,:]
            u_xp1 = @view fields[idx_xp1, idx_y,:]
            u_xp2 = @view fields[idx_xp2, idx_y,:]
            u_xp3 = @view fields[idx_xp3, idx_y,:]
            

            u_ym3 = @view fields[idx_x, idx_ym3,:]
            u_ym2 = @view fields[idx_x, idx_ym2,:]
            u_ym1 = @view fields[idx_x, idx_ym1,:]
            u_y   = @view fields[idx_x, idx_y  ,:]
            u_yp1 = @view fields[idx_x, idx_yp1,:]
            u_yp2 = @view fields[idx_x, idx_yp2,:]
            u_yp3 = @view fields[idx_x, idx_yp3,:]
            
            S_MAX_X = max(Speed_max(u_xp3, par_eq), Speed_max(u_xm3, par_eq), 
                Speed_max(u_xp2, par_eq), Speed_max(u_xm2, par_eq), Speed_max(u_xp1, par_eq), 
                Speed_max(u_xm1, par_eq), Speed_max(u_x, par_eq)) #maximum speed
            
            S_MAX_Y = max(Speed_max(u_yp3, par_eq), Speed_max(u_ym3, par_eq), 
                Speed_max(u_yp2, par_eq), Speed_max(u_ym2, par_eq), Speed_max(u_yp1, par_eq), 
                Speed_max(u_ym1, par_eq), Speed_max(u_y, par_eq)) #maximum speed
            
            
            Fx!(Fx_Pm3, u_xm3, idx_xm3, idx_y, par_eq)
            Fx!(Fx_Pm2, u_xm2, idx_xm2, idx_y, par_eq)
            Fx!(Fx_Pm1, u_xm1, idx_xm1, idx_y, par_eq)
            Fx!(Fx_P  , u_x  , idx_x  , idx_y, par_eq)
            Fx!(Fx_Pp1, u_xp1, idx_xp1, idx_y, par_eq)
            Fx!(Fx_Pp2, u_xp2, idx_xp2, idx_y, par_eq)
            Fx!(Fx_Pp3, u_xp3, idx_xp3, idx_y, par_eq)
            
            Fy!(Fy_Pm3, u_ym3, idx_x, idx_ym3, par_eq)
            Fy!(Fy_Pm2, u_ym2, idx_x, idx_ym2, par_eq)
            Fy!(Fy_Pm1, u_ym1, idx_x, idx_ym1, par_eq)
            Fy!(Fy_P  , u_y  , idx_x, idx_y  , par_eq)
            Fy!(Fy_Pp1, u_yp1, idx_x, idx_yp1, par_eq)
            Fy!(Fy_Pp2, u_yp2, idx_x, idx_yp2, par_eq)
            Fy!(Fy_Pp3, u_yp3, idx_x, idx_yp3, par_eq)
            
            
            @. Fx_Mm3 = 0.5 * (Fx_Pm3 - S_MAX_X * u_xm3)
            @. Fx_Mm2 = 0.5 * (Fx_Pm2 - S_MAX_X * u_xm2)
            @. Fx_Mm1 = 0.5 * (Fx_Pm1 - S_MAX_X * u_xm1)
            @. Fx_M   = 0.5 * (Fx_P   - S_MAX_X * u_x)
            @. Fx_Mp1 = 0.5 * (Fx_Pp1 - S_MAX_X * u_xp1)
            @. Fx_Mp2 = 0.5 * (Fx_Pp2 - S_MAX_X * u_xp2)
            @. Fx_Mp3 = 0.5 * (Fx_Pp3 - S_MAX_X * u_xp3)
            @. Fx_Pm3 = 0.5 * (Fx_Pm3 + S_MAX_X * u_xm3)
            @. Fx_Pm2 = 0.5 * (Fx_Pm2 + S_MAX_X * u_xm2)
            @. Fx_Pm1 = 0.5 * (Fx_Pm1 + S_MAX_X * u_xm1)
            @. Fx_P   = 0.5 * (Fx_P   + S_MAX_X * u_x)
            @. Fx_Pp1 = 0.5 * (Fx_Pp1 + S_MAX_X * u_xp1)
            @. Fx_Pp2 = 0.5 * (Fx_Pp2 + S_MAX_X * u_xp2)
            @. Fx_Pp3 = 0.5 * (Fx_Pp3 + S_MAX_X * u_xp3)
            
            @. Fy_Mm3 = 0.5 * (Fy_Pm3 - S_MAX_Y * u_ym3)
            @. Fy_Mm2 = 0.5 * (Fy_Pm2 - S_MAX_Y * u_ym2)
            @. Fy_Mm1 = 0.5 * (Fy_Pm1 - S_MAX_Y * u_ym1)
            @. Fy_M   = 0.5 * (Fy_P   - S_MAX_Y * u_y)
            @. Fy_Mp1 = 0.5 * (Fy_Pp1 - S_MAX_Y * u_yp1)
            @. Fy_Mp2 = 0.5 * (Fy_Pp2 - S_MAX_Y * u_yp2)
            @. Fy_Mp3 = 0.5 * (Fy_Pp3 - S_MAX_Y * u_yp3)
            @. Fy_Pm3 = 0.5 * (Fy_Pm3 + S_MAX_Y * u_ym3)
            @. Fy_Pm2 = 0.5 * (Fy_Pm2 + S_MAX_Y * u_ym2)
            @. Fy_Pm1 = 0.5 * (Fy_Pm1 + S_MAX_Y * u_ym1)
            @. Fy_P   = 0.5 * (Fy_P   + S_MAX_Y * u_y)
            @. Fy_Pp1 = 0.5 * (Fy_Pp1 + S_MAX_Y * u_yp1)
            @. Fy_Pp2 = 0.5 * (Fy_Pp2 + S_MAX_Y * u_yp2)
            @. Fy_Pp3 = 0.5 * (Fy_Pp3 + S_MAX_Y * u_yp3)
            #println("doing reconstruction...")

            
            
            WENOZreconstruction!(Fx_RM, Fx_Mp2, Fx_Mp1, Fx_M  , Fx_Mm1, Fx_Mm2, N_Fields)
            WENOZreconstruction!(Fx_LM, Fx_Pm3, Fx_Pm2, Fx_Pm1, Fx_P  , Fx_Pp1, N_Fields)
            WENOZreconstruction!(Fx_LP, Fx_Pm2, Fx_Pm1, Fx_P  , Fx_Pp1, Fx_Pp2, N_Fields)
            WENOZreconstruction!(Fx_RP, Fx_Mp3, Fx_Mp2, Fx_Mp1, Fx_M  , Fx_Mm1, N_Fields)
            
            WENOZreconstruction!(Fy_RM, Fy_Mp2, Fy_Mp1, Fy_M  , Fy_Mm1, Fy_Mm2, N_Fields)
            WENOZreconstruction!(Fy_LM, Fy_Pm3, Fy_Pm2, Fy_Pm1, Fy_P  , Fy_Pp1, N_Fields)
            WENOZreconstruction!(Fy_LP, Fy_Pm2, Fy_Pm1, Fy_P  , Fy_Pp1, Fy_Pp2, N_Fields)
            WENOZreconstruction!(Fy_RP, Fy_Mp3, Fy_Mp2, Fy_Mp1, Fy_M  , Fy_Mm1, N_Fields)
            
            #println("done!")
            @. Hx_p = Fx_LP + Fx_RP
            @. Hx_m = Fx_LM + Fx_RM
            #=if idx_x == 4
                if idx_y == 10
                    
                    println("hello there")
                    println(Fx_M[12])
                    println(Hx_p[12])
                    println(Fx_LP[12])
                    println(Fx_RP[12])
                    println(Hx_m,[12])
                    println(Fx_LM[12])
                    println(Fx_RM[12])
                end
            end
            =#

            @. Hy_p = Fy_LP + Fy_RP
            @. Hy_m = Fy_LM + Fy_RM
            #println("size of T = ", size(T))
            #println("size of A = ", size(A))
            #println("calculating TandA...")
            TandA!(T, A, u_x, idx_x, idx_y, par_eq)
            #println("TandA calculated")
            #println("about to enter the source...")
            source!(sourcevec, u_x[15:end], T, A, idx_x, idx_y, t, par_source) #u_x = u_y = fields[idx_x, idx_y, :]
            #println("Source ended")
            #=
            if ((idx_x == 4) && (idx_y == 25))
                println(A[1,3,3])
                println(A[2,3,3])
                println(Fy_RM[12])
                println(Fy_LM[12])
                println(Fy_RP[12])
                println(Fy_LP[12])
            end
            =#
            @. dfields[idx_x, idx_y, :] = (-hx)*(Hx_p - Hx_m) #+ (-hy)*(Hy_p - Hy_m) + sourcevec
        end
    end
    #Wall at x = 0
    for idx_x in 1:3
        #dfields[4-idx_x,:,:] .= 0.0#
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

        
        dfields[4-idx_x,:, 1]    .=  dfields[3+idx_x,:, 1]  #Energy has the same sign in both sides
        dfields[4-idx_x,:, 2:3]  .= 0.0#-dfields[3+idx_x,:, 2:3] #Momenta has oposite signs.
        dfields[4-idx_x,:, 4]    .= 0.0#dfields[3+idx_x,:, 4] #Momenta has oposite signs.
        
        
        dfields[4-idx_x,:, idA000]    .=  dfields[3+idx_x,:, idA000] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA001]    .= -dfields[3+idx_x,:, idA001] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA002]    .= -dfields[3+idx_x,:, idA002] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA003]    .=  dfields[3+idx_x,:, idA003] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA011]    .=  dfields[3+idx_x,:, idA011] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA012]    .=  dfields[3+idx_x,:, idA012] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA013]    .= -dfields[3+idx_x,:, idA013] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA022]    .=  dfields[3+idx_x,:, idA022] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA023]    .= -dfields[3+idx_x,:, idA023] #Momenta has oposite signs.
        dfields[4-idx_x,:, idA033]    .=  dfields[3+idx_x,:, idA033] #Momenta has oposite signs.

        dfields[4-idx_x,:,15:end] .= 0.0
    end
    #Free boundary and disk at x = N
    for idx_y in 1:M+6
        for idx_x in 1:3
            x, z = get_pos(N+3+idx_x, idx_y , lengthpars)
            dfields[N+3+idx_x,idx_y, 1:14] .= dfields[N+3,idx_y, 1:14]  #Energy has the same sign in both sides
            dfields[N+3+idx_x,idx_y, 15:end] .= 0.0  #We do not touch the abstract fields
        end
    end
    #Free boundary at other two walls
    
    for idx_y in 1:3
        dfields[:,4-idx_y, :] .= dfields[:,4, :]  #Energy has the same sign in both sides
        dfields[:,3+M+idx_y, :] .= dfields[:,M+3, :]  #We do not touch the abstract fields
    end
end
=#


function evolution!(du,u,par,t)
    #println(u[50,1,1])
    #println("entering evolution")
    par_inv, par_ev, lengthpars = par
    Lx, Ly, Lz, N, M, O = lengthpars
    #println("starting inversion")
    
    for i in 1:N+6
        for j in 1:M+6
            for k in 1:O+6
                #println("inverting index i = $i, j = $j")
                convars = @view u[i,j,k,1:14]
                absvars = @view u[i,j,k,15:end]
                
                c2abs!(convars, absvars, i, j, k, par_inv)
            end
        end
    end
    #println("inversion ended")
    wenoz!(du, u, par_ev, t)
    return du[:]
end


#=
function evolution!(du,u,par,t)
    #println(u[50,1,1])
    #println("entering evolution")
    par_inv, par_ev = par
    lengthpars = par_ev[1]
    Lx, Ly, N, M, N_Fields = lengthpars;
    for i in 1:N
        for j in 1:M
            #println("inverting index i = $i, j = $j")
            convars = @view u[i,j,1:N_Fields]
            absvars = @view u[i,j,N_Fields+1:end]
            c2abs!(convars, absvars, i, j, par_inv)
        end
    end
    #println("invertion ended!")
    wenoz!(du, u, par_ev, t)
    return du[:]
end
=#

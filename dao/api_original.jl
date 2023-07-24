#generate Recfast params
function adm_params(xi, rm, rM, ra)
    p = Recfast.Params(Yp = 0.0,  T0 = 2.726* xi, Omega_M =Ωm, Omega_B = (Ωm)/ rM, h100 = h, n_eff = 7.456/xi^4)
    Recfast.Set_VFC_params!(p, aS = ra, mS = rm, mxS = rM)
    return p 
end


##halo mass function, retrun dn/dlnM (comoving h^3/Mpc^3)
function hmf(p::Recfast.Params,  z)
    zrec = zRec(p)
    sol = Recfast.Evaluate_recombination_h2(p, logzstart = log10(zrec * 2), logzend = log10(max(z, zrec/1e3)), dt = 1e-2, dtmin = 1e-30)
    kdscale = 1/(kd(p,sol) * mtoMpc *h) #need to insert h, and compare to dao scale
    
    sm = get_linear_power(p, 0:z/10:z)
    println("kd ", 1/kdscale)
    return m -> dndmsharpk.(p, x-> sm.(z,x).* diffusion_damp.(x, kdscale), m)
end


#for cdm (no damping)
function hmfcdm(p::Recfast.Params, z)
    zrec = zRec(p)
    sol = Recfast.Evaluate_recombination_h2(p, logzstart = log10(zrec * 2), logzend = log10(max(z, zrec/1e3)), dt = 1e-2, dtmin = 1e-30)
    
    sm = get_linear_power(p, p.Omega_B, 0:z/10:z)
    return m -> dndmsharpk.(p, x-> sm.(z,x), m)
end

#specify pk, instead of calculating it
function hmf(p::Recfast.Params, pk, z)
    zrec = zRec(p)
    sol = Recfast.Evaluate_recombination_h2(p, logzstart = log10(zrec * 2), logzend = log10(max(z, zrec/1e3)), dt = 1e-2, dtmin = 1e-30)
    println(kd(p,sol))
    kdscale = 1/(h * kd(p,sol) *mtoMpc)
    
            
    println("kd ", 1/kdscale)
    return m -> dndmsharpk.(p, x-> pk.(z,x).* diffusion_damp.(x, kdscale), m)
end



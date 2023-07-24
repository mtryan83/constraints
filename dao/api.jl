import .DAOHMF_mod as dh
include("../../RecfastJulia/src/recfast.jl")
using .Recfast

#generate Recfast params
function adm_params(xi, rm, rM, ra)
    p = Recfast.Params(Yp = 0.0,  T0 = 2.726* xi, Omega_M =dh.Ωm, Omega_B = (dh.Ωm)/ rM, h100 = dh.h, n_eff = 7.456/xi^4)
    Recfast.Set_VFC_params!(p, aS = ra, mS = rm, mxS = rM)
    return p 
end

# MR: Added version to use ADM struct. Note that this version is preferred
function adm_params(adm::ADM)
    return adm_params(adm.ξ,adm.rm,adm.rM,adm.rα)
end

##halo mass function, return dn/dlnM (comoving h^3/Mpc^3)
# function hmf(p::Recfast.Params,  z)
#    zrec = zRec(p)
#    sol = Recfast.Evaluate_recombination_h2(p, logzstart = log10(zrec * 2), logzend = log10(max(z, zrec/1e3)), dt = 1e-2, dtmin = 1e-30)
#    kdscale = 1/(kd(p,sol) * mtoMpc *h) #need to insert h, and compare to dao scale
   
#    sm = get_linear_power(p, 0:z/10:z)
#    println("kd ", 1/kdscale)
#    return m -> dndmsharpk.(p, x-> sm.(z,x).* diffusion_damp.(x, kdscale), m)
# end

#specify pk, instead of calculating it
# MR: This has been changed from the email to add default for pk - means this function is identical to one above
function hmf(p::Recfast.Params, z; pk=nothing)
    zrec = dh.zRec(p)
    sol = Recfast.Evaluate_recombination_h2(p, logzstart = log10(zrec * 2), logzend = log10(max(z, zrec/1e3)), dt = 1e-2, dtmin = 1e-30)
    println("kd(p,sol): ",dh.kd(p,sol) * dh.mtoMpc)
    kdscale = 1/(dh.h * dh.kd(p,sol) *dh.mtoMpc)
    
    if isnothing(pk)
        pk = dh.get_linear_power(p, 0:z/10:z)
    end
    println("kd ", 1/kdscale)
    return m -> dh.dndmsharpk.(p, x-> pk.(z,x).* dh.diffusion_damp.(x, kdscale), m)
end

# MR: Added version to use ADM struct
function hmf(adm::ADM, z; pk=nothing)
    return hmf(adm_params(adm),z;pk)
end


#for cdm (no damping)
function hmfcdm(p::Recfast.Params, z)
    zrec = dh.zRec(p)
    sol = Recfast.Evaluate_recombination_h2(p, logzstart = log10(zrec * 2), logzend = log10(max(z, zrec/1e3)), dt = 1e-2, dtmin = 1e-30)
    
    sm = dh.get_linear_power(p, 0:z/10:z, p.Omega_B)
    return m -> dh.dndmsharpk.(p, x-> sm.(z,x), m)
end

# MR: Added version to use ADM struct.
function hmfcdm(adm::ADM, z)
    return hmfcdm(adm_params(adm),z)
end




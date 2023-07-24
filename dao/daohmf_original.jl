using PyCall
using Dierckx
camb = pyimport("camb")
using FastGaussQuadrature
using LinearAlgebra
using QuadGK
include("../RecfastJulia/src/Recfast.jl")
using .Recfast
using Roots
const h = .6774
const Ωb = 0.0223/h^2
const Ωm = .3089
const Ωr = Ωm / ( 3372)
const TCMB = 2.72
const G = 4.30092e-9 #Mpc solar mass^-1 (km/s)^2
const kb = 8.62e-5#eV/K
const kbsi   =1.3806503e-23 #joule/K
const hPlanck = 6.62606876e-34 ## Planck constant [Js]
const mpev = 1e9
const ckms = 3e5
const mtoMpc = 3.24078e-23 #meters to Mpc
const cLight  = 2.99792458E+08
const sigmaT = 6.6524616e-29
const ρc = 2.775e11 #h^2 M_\odot Mpc^-3

const cLight  = 2.99792458e+08
const sigmaT = 6.6524616e-29

#calc power spectrum for dark parameters and specified \Omega_b
function get_linear_power(p, zvals)
    pars = camb.CAMBparams()
    pars.set_cosmology(H0 = 100*p.h100, ombh2 =Ωb* p.h100^2, omch2 = (p.Omega_M - Ωb) * p.h100^2, omk=0)
    pars.NonLinear = camb.model.NonLinear_none
    pars.set_matter_power(redshifts=zvals, kmax=1000.0)
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1000, npoints = 1000)

    println(results.get_sigma8())
    println(results.get_sigmaR(R= 8, hubble_units=false))
    return Spline2D(z, kh   , pk)
end

function diffusion_damp(k, kd)
    return exp(-2*(k/kd)^2)
end

#diffusion scale, in comoving meters
function kd(p, sol)
    R(a) = 2.32e10*0.02^4 * (p.T0/2.72)^-4 *a #The inertia term. Eq 3.3 of dark recombination
    dsigmaT = p.rescaling.RF_mScale^-2 * p.rescaling.RF_aScale^2 * sigmaT
    zdec = zDec(p, sol)
    adecp = 1/(zdec + 1)
    xc(a) = sol(1/a-1)[2]
    nz(a) = Recfast.NH(p, 1/a-1)
    integrand(a) = cLight/(Recfast.H_z_loc(p, 1/a - 1)
        * a^2 * 6 * (1 + R(a)) * nz(a)*xc(a) * a *dsigmaT) * (16/15 + R(a)^2/(1 + R(a)))
    return quadgk(integrand, 0, adecp)[1]^0.5
end

#dark acoustic scale, in comoving meters
function dao(p, sol)
    R(a) = 2.32e10*0.02^4 * a * (p.T0/2.72)^-4
    adec = 1/(zDec(p, sol) + 1)
    cs(a) = cLight/(3 * (1 + R(a)))^(1/2)
    return quadgk(a -> (cs(a) /(a^2 * Recfast.H_z_loc(p, (1/a -1)))), 0, adec)[1]
end

xcSaha(p, z) = (- saha(p, z) + sqrt(saha(p, z)^2 + 4 * saha(p, z)))/2

#recombination redshift, saha
function zRec(p)
    lb = 1100 * p.rescaling.fl * (2.72/p.T0) * 1e-3
    ub = 1100 * p.rescaling.fl * (2.72/p.T0) * 1e3
    zrec = find_zero(x->xcSaha(p, x) .- .5, (lb, ub))
    return zrec
end



#the saha ionization ratio
function saha(p, z)
    mc = p.rescaling.RF_mScale * 9.10938188E-31

    rydberg = 1.096787737e+7 * cLight * hPlanck * p.rescaling.fl

    T = (1 + z) * p.T0
    sh= (Recfast.NH(p, z))^(-1) * (2 * pi * mc * kbsi * T/(hPlanck^2))^(3/2) * exp(-rydberg/(kbsi * T))
    return sh
end


#thermal decoupling redshift
function zDec(p, sol)
        rec = zRec(p)

        xc(z) = sol(z)[2]
        dsigmaT = p.rescaling.fC * sigmaT #the dark one
        nz(z) = Recfast.NH(p,z)
        g(z) = cLight * dsigmaT * nz(z) * xc(z)/((1 + z) * Recfast.H_z_loc(p, z)) * exp(-dsigmaT * cLight
            * quadgk(x -> (nz(x) * xc(x)/((1 + x) * Recfast.H_z_loc(p, x))), 0, z)[1])

        zs = 10 .^ range(log10(rec/10), stop = log10(1000 * rec), length = 20)

        for i = 1:10
            guess =  argmax(g.(zs))
            zminind = max(1, guess-1)
            zmaxind = min(20, guess+1)
            zs = 10 .^ range(log10(zs[zminind]), log10(zs[zmaxind]), length = 20)
        end

        return zs[argmax(g.(zs))]
    end

#calculate linear power spectrum. p.Omega_b is set to get the correct adm number density, so the real baryon density must be separatelhy specified 
function get_linear_power(p, zvals)
    pars = camb.CAMBparams()
    pars.set_cosmology(H0 = 100*p.h100, ombh2 =Ωb* p.h100^2, omch2 = (p.Omega_M - Ωb) * p.h100^2, omk=0)
    pars.NonLinear = camb.model.NonLinear_none
    pars.set_matter_power(redshifts=zvals, kmax=1000.0)
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1000, npoints = 1000)

    println(results.get_sigma8())
    println(results.get_sigmaR(R= 8, hubble_units=false))
    return Spline2D(z, kh   , pk)
end

#gaussian cutoff
function diffusion_damp(k, kd)
    return exp(-2*(k/kd)^2)
end


#for the window function
j1(x) = (sin(x) - x * cos(x))/x^2
w(k,r) = 3 * j1(k * r)/(k * r)
#the mass/radius assigment from (I think) Benson 2013
rsk(ρ, m) = (3/(4 * pi * ρ)* m)^(1/3)/2.7


#sharp k filter
function sharpk(k, r)
    if (k < 1. /r)
        return 1.
    end
    return 0.
end

function σ(pk, r; filter = w)
    nodes, weights = gausslegendre(2000)
    integrand(k) = 1/(2 * pi^2) *k^2 *pk(k) * filter(k, r)^2
    kmin = 0
    kmax = 10 ./ r
    transformedInt(k) = integrand(kmax/2 * (k + 1)) * kmax/2
    return dot(weights, transformedInt.(nodes) )^0.5

end

function σ(pk, r, kmax, filter = w)
    nodes, weights = gausslegendre(2000)
    integrand(k) = 1/(2 * pi^2) *k^2 *pk(k) * filter(k, r)^2
    transformedInt(x) = integrand((kmax-kmin)/2 * (x+ (2*kmin/(kmax-kmin)) + 1)) * (kmax-kmin)/2
    return dot(weights, transformedInt.(nodes) )^0.5
end




function shethtormensk(σ)
    δc = 1.69
    A = 0.3222
    a = 1.
    p = 0.3
    #f(σ) = sqrt(2/pi) * δc/σ * exp(-(δc/σ)^2/2)
    return A *sqrt(2*a/pi) * (1 + (σ^2/(a*δc^2))^p) * (δc/σ)* exp(-a * δc^2/2σ^2)
end




function dndmsharpk(p, pk, M)
    ρm = ρc * p.Omega_M  #background matter density at z=0, M_sol/Mpc^3. The value differs from Donghui's class notes


    radius = rsk(ρm, M)
    #This is d\sqrt{\sigma^2} = 1/2 * sigma^2^-1/2 * d sigma^2/d(1/r) * d(1/r) * dr/dM
    dσ = 0.5 / σ(pk, radius, filter = sharpk) * (1/(2 * pi^2) *(1/radius)^2 *pk(1/radius)) * ((1/3) * radius/M)/ ( radius^2)
    ld = M * abs(dσ)/σ(pk, radius, filter = sharpk) #dlog \sigma / dlog m
  
    return ld * shethtormensk(σ(pk, radius, filter = sharpk)) *ρm/M #dn/d ln M

end
using LazyGrids
using PyPlot
using Statistics
using Interpolations

# Someday it would be nice to use the following for the Properties
#using Measurements

import FromFile: @from

@from "adm.jl" using ADM_mod
@from "Get_Lambda.jl" using Get_Lambda
#include("Get_Lambda.jl")
#using .Get_Lambda

# Defining constants
    
# Defining global constants
# Move constants from above if they end up being used in subfunctions
# Hopefully, these should just be fundamental physical constants
const G_cgs = 6.67259e-8 # cgs 
const kb_keV = 8.617e-8 # keV/K
const kb_geV = 8.617e-14 # GeV/K
const kb_erg = 1.38065e-16 # erg/K
const planck = 4.135667696e-18 # keV*s
const age_of_universe = 13.797 # Gyr
const c_cm = 2.998e10 # cm/s
const c_km = 2.998e5 #  km/s
const kev_to_erg = 1.602e-9

# Using Planck 2018 data for cosmological constants
h = 0.67019
H0_base = 100*h # in (km/s)/Mpc
H0 = H0_base*3.24078e-20 # 1/s
OmegaM = 0.12029/h^2
OmegaB = 0.0222/h^2
OmegaR = 2.488e-5/h^2 
OmegaDM = OmegaM-OmegaB
OmegaADM(epsilon) = epsilon * OmegaDM
OmegaL = 0.6889
Omega = OmegaM


delta = 18*pi^2


const s_to_gy = 3.17098e-17 
const gy_to_s = 3.154e16
const m_sol_to_g = 1.9891e33 
const g_to_gev = 5.6096e23
const gev_to_g = 1/g_to_gev
    
# Proton and electron masses
#mP = 0.938 # GeV
#mE = 511   # keV
#alpha = 0.00729
    
    
H(a) = H0*sqrt(OmegaL + OmegaM./a.^3 + OmegaR./a.^4 + (1-OmegaL-OmegaM-OmegaR)./a.^2)
        
rho_cr(z) = 3*H(1.0./(1+z)).^2/(8*pi*G_cgs)
    
delta_c = 18*pi^2
rho_0 = 3*H0^2/(8*pi*G_cgs)
rho_v(z) = OmegaDM*rho_cr(z)


# Mass-loss type constraint 
# Can provide either E_tot and x_DM or E_dens and n_DM
# returns Lambda
function m_loss_constraint(f_lost,v_coll,Sigma_s,m_H;E_tot=nothing,E_dens=nothing,n_DM=nothing,x_DM=nothing)
    if isnothing(E_dens) && isnothing(E_tot)
        throw(ArgumentError("must provide either E_dens or E_tot"))
    elseif isnothing(E_dens)
        Lambda = ((f_lost .* E_tot .* v_coll) ./ Sigma_s) .* (m_H ./ x_DM)
    else
        Lambda = ((f_lost .* E_dens .* v_coll) ./ Sigma_s) .* (m_H ./ n_DM)
    end
    return Lambda
end

# Structure formation type constraint
# returns Lambda
function struc_form_constraint(E_dens,tff,nDM)
    Lambda = ((E_dens) ./ tff) .* (1 ./ nDM) .^ 2
    return Lambda
end

# Calculate freefall time
# returns tff
function tff(rho)
    t_ff = sqrt(1 ./(G_cgs .* rho))
    return t_ff
end

# Calculate total energy from temperature
# return E_tot in ergs
function E_tot(T;gamma=nothing)
    if isnothing(gamma)
        gamma = 5/3
    end
    return kb_erg * T / (gamma - 1)
end

# Calculate total energy density from temperature and number density
# return E_dens in ergs/cm^3
function E_dens(T,n;gamma=nothing)
    if isnothing(gamma)
        gamma = 5/3
    end
    return n .* kb_erg .* T / (gamma - 1)
end

# Calculate virial temp from halo mass
# Ideally, this would be replaced by a library call
# Note that rho and mh need to be in CGS, Mv in Msol
# returns T_virial
function Tv(rho,Mv,mh;dens_prof=3/5)
    Mv = Mv * m_sol_to_g # convert to g
    Tv = (4 / 3 * pi * rho)^(1 / 3) * 1 / 2 * dens_prof * mh * G_cgs * Mv^(2 / 3) / kb_erg
end

# Calculate number density as a function of z assuming 3 sigma overdensity
# Comes from the paragraph above eq 8.92 in MBW
# inputs: mh in GeV, redshift
# returns n in cm^-3
function nofz(mh,z;f=OmegaADM(1),mu=1,delta=delta)
    # f = M_halo_DM/M_halo * epsilon
    return 1.05367e-5 .* f ./ (mu .* mh) .* (1 .+ delta).*(Omega .* h.^2) .* (1 .+ z).^3
end
    
@enum cons_dir above below # direction of constraint - above means constrains from above=upper limit on Lambda
EmptyParam=Union{Float64,Array{Float64},Some{Array{Float64}},Nothing}
EmptyFun=Union{Function,Some{Function},Nothing}
@kwdef struct Properties
    Mhalo = nothing
    z = nothing
    Sigma_T = nothing
    v_c=nothing
    f_lost=nothing
    direction::cons_dir=above
    x_D::EmptyFun=nothing
    allow_any::Bool=false
end

#Properties()


# Get LSS properties from Hennig et al 2017
# returns M_scale, z, Sigma_T, v_c, f_lost, direction, nothing
function getLSSProperties()
    M=[6.3000,5.2000,5.9000,6.1000,17.5000,10.2000,25.7000,5.3000,11.7000,
            6.1000,13.2000,5.4000,6.5000,5.2000,4.5000,8.8000,5.5000,5.2000,
            6.3000,7.6000,4.3000,5.8000,5.4000,18.7000,9.3000,6.5000,8.2000,
            8.8000,4.6000,5.7000,5.2000,7.2000,6.9000,9.1000,5.6000,5.2000,6.3000,
            5.1000,5.5000,4.8000,6.0000,7.2000,5.2000,5.4000,6.6000,8.5000,5.8000,
            9.1000,11.4000,12.3000,5.8000,5.9000,7.9000,6.3000,6.2000,5.2000,5.1000,7.1000,
            9.5000,5.2000,21.1000,6.2000,13.2000,7.0000,28.0000,28.9000,4.7000,5.8000,6.5000,
            5.0000,9.3000,7.7000,6.8000,6.1000]*1e14
    #M = [M,7e13]
    a0 = delta_c*rho_0
    a1 = 5.77
    function truncatedSigma(a0,a1,x,c)
        if x<1
            F = (sqrt(c^2-x^2)/(1+c) - acosh((c+x^2)/(x+x*c))/sqrt(1-x^2))/(x^2-1)
        else
            F = (sqrt(c^2-x^2)/(1+c) - acos((c+x^2)/(x+x*c))/sqrt(x^2-1))/(x^2-1)
        end
        f = 2*a0*a1*F
        return f
    end
    #println("truncated sigma:$(truncatedSigma(a0,a1,1e-30,5.77))")
    z=[0.89,0.39,0.4,0.87,0.39,0.36,0.88,0.6,0.45,0.28,0.6,0.67,0.6,
            0.63,1.05,0.66,0.75,0.53,0.23,0.07,1.12,0.7,0.28,0.42,0.39,
            0.43,0.4,0.8,0.95,0.58,0.99,1.11,0.42,0.8,0.41,0.44,0.88,0.65,
            0.43,0.2,0.24,0.15,0.83,0.8,0.22,0.29,0.43,0.46,0.4,0.29,0.91,
            0.18,0.81,0.33,0.75,0.58,0.85,0.75,0.48,0.81,0.42,0.74,0.4,
            0.31,0.33,0.37,0.75,0.45,0.17,0.76,0.42,0.98,0.43,0.55]
    #z = [z,0.296]
    # Total observed surface mass density of the Bullet Cluster(?)
    Sigma_T = 0.25 # cm^-2 g, 
    v_c = 1000 * 1000 * 100 # 1000 km/s
    
    f_lost = [0.01,0.1,0.3]
    
    direction = above
    
    # yes there's a better way to do this and I can't figure it out at the moment
    M, = ndgrid(M,f_lost)
    z,f_lost = ndgrid(z,f_lost)
    
    props = Properties(;Mhalo=M,z,Sigma_T,v_c,f_lost,direction)
    
    return props
end
    
# Get DwGal properties from McConnachie 2012 and Herrmann 2016
# returns M_scale, z, Sigma_T, v_c, f_lost, direction, nothing
function getDwGalProperties()
    M=sort([1.9e+08, 3.9e+06, 3.3e+06, 230000, 270000, 940000, 810000,
        1.1e+07, 9.5e+06, 1.4e+07, 2.5e+07, 1.1e+07, 6.3e+06, 2.6e+06,
        5.6e+07, 1.3e+06, 910000, 1.1e+06, 1.9e+07, 4.6e+06, 1.2e+07,
        5.4e+08, 6.5e+06, 4.2e+08, 4.4e+07, 6.1e+06, 2.3e+06, 1.2e+06,
        9.3e+07, 6.1e+06, 1.6e+07, 1.1e+07, 3.6e+07, 1.5e+08, 4.2e+07,
        1.7e+07, 3.9e+06, 3.8e+08, 2.5e+07, 4.1e+07,])
    z=zeros(size(M))
    Sigma_T = 1e-3  # Approximation from Herrmann 2016
    v_c = 100 * 1000 * 100# km/s
    
    f_lost = [0.5]
    
    direction = below
    
    M, = ndgrid(M,f_lost)
    z,f_lost = ndgrid(z,f_lost)
    
    props = Properties(;Mhalo=M,z,Sigma_T,v_c,f_lost,direction)
    return props
end
    
# Get DBC properties from Lee 2021
# returns M_scale, z, Sigma_T, v_c, f_lost, direction, x_D(rM)
function getDBCProperties()
    Md = [3, 3.88, 5]*1e9
    z = [0.02]
    Sigma_T = [1,10] # Msol/pc^2
    Sigma_T = Sigma_T*2.088e-4 # g/cm^2
    vs = [100:100:900;] * 1000 * 100 # Note the ; converts this into an Array!
    Mgas = [0.6, 1.68, 2.2]*1e9
    f_lost = [0.01,0.1,0.3]
    
    # Still looking for a better way of doing this. Maybe just list comprehensions?
    Md,zSvMgf = ndgrid(Md,vec(collect(Iterators.product(z,Sigma_T,vs,Mgas,f_lost))))
    
    M = Md
    z=getindex.(zSvMgf,1)
    Sigma_T=getindex.(zSvMgf,2)
    vs=getindex.(zSvMgf,3)
    Mgas = getindex.(zSvMgf,4)
    f_lost = getindex.(zSvMgf,5)
    xDf(rX) = Md ./ (Md .+ rX .* Mgas)
    
    direction = above
    
    props = Properties(M,z,Sigma_T,vs,f_lost,direction,xDf,false)
    return props
end

# Check that Milky Way doesn't collapse
function getMWProperties()
    # Only need Mhalo, z, and direction for struc formation constraint
    # from 1804.11348
    Mhalo = 1.54e12
    # interested in structure formation era, so z=5,10,15?
    z = [5 10 15]
    direction = above
    return Properties(;Mhalo,z,direction)
end

# Based off of 2009.05209
function getDCOProperties()
    # 10^6 - 10^9 halos collapse
    Mhalo = 10 .^ range(6,9,length=10)
    # structure formation era
    z = [5 10]
    direction = below
    # Constraint passed if any halo passed
    allow_any=true
    return Properties(;Mhalo,z,direction,allow_any)
end

# Calculate the mass-loss constraint specified by properties
# returns Lambda_cons, direction::cons_dir, n, T
function calcMassLossConstraint(properties::Properties,adm::ADM)
    # properties() will return Mhalo, z, Sigma_T, v_c, f_lost, direction, and optionally x_D
    (;Mhalo, z, Sigma_T, v_c, f_lost, direction, x_D) = properties
    #println("Properties: Mhalo:$Mhalo z:$z Sigma_T:$Sigma_T v_c:$v_c f_lost:$f_lost direction:$direction x_D:$x_D")
    n = nofz(adm.M,z)
    T = Tv.(adm.M_mass .* n, Mhalo, adm.M_mass)
    # Check if we have x_D - if so, use E_tot+x_d, otherwise E_dens+n_dm
    #f_lost, v_coll, Sigma_s, m_H; E_tot+x_DM,E_dens+n_DM
    has_x = !isnothing(x_D)
    Lambda_cons = m_loss_constraint(f_lost,v_c,Sigma_T,adm.M_mass;
        E_tot=has_x ? E_tot(T) : nothing,
        E_dens=!has_x ? E_dens(n,T) : nothing,
        x_DM=has_x ? x_D(adm.rM) : nothing,
        n_DM=!has_x ? n : nothing         # Revisit this. It's supposed to be n_DM
    )
    
    # This currently can produce upwards of 100 lines for certain properties 
    # (Looking at you getDBCProperties). Which takes forever to run 
    # get_Lambda on. So we'll reduce to 3 lines - max at every M, min at every
    # M, and median
    if size(Lambda_cons,2)>3
        max_inds = argmax(Lambda_cons,dims=2)
        med = median(Lambda_cons,dims=2)
        med_inds = argmin(abs.(Lambda_cons.-med),dims=2)
        min_inds = argmin(Lambda_cons,dims=2)
        inds = [min_inds med_inds max_inds]
        return @views Lambda_cons[inds], direction, n[inds], T[inds]
    end
    
    return Lambda_cons,direction,n,T
end

function calcMassLossConstraint(properties::Function,adm::ADM)
    return calcMassLossConstraint(properties(),adm)
end

# Check if provided cooling violates the mass-loss constraint specified
# by properties
function checkMassLossConstraint(Lambda,properties,adm::ADM;n_test=nothing)
    Lambda_cons, direction, n, T = calcMassLossConstraint(properties,adm)
    if isnothing(n_test)
        Lambda_test = Lambda.(n,T;adm)
    else
        Lambda_test = Lambda.(n_test,T;adm)
    end
    # Could write this with short-circuiting booleans, but this is cleaner
    check = direction==above ? Lambda_test.<=Lambda_cons : Lambda_test.>=Lambda_cons
    if ndims(check)==2
        acheck = []
        for i = 1:size(check)[2]
            push!(acheck,all(check[:,i]))
        end
        return acheck
    elseif ndims(check)>2
        throw(ErrorException("I don't know how to handle $(ndims(check)) dimensions"))
    else
        return all(check)
    end
end

# Calculate the structure-formation constraint specified by properties
# returns Lambda_cons, direction::cons_dir, n, T
function calcStrucFormConstraint(properties::Properties,adm::ADM)
    (;Mhalo, z, direction) = properties
    n = nofz(adm.M,z)
    T = Tv.(adm.M_mass .* n,Mhalo, adm.M_mass;dens_prof=2) # dens_prof=2 is equivalent to Buckley version
    rho = adm.M_mass * n  # note this should include a factor of mu, but we'll assume it's 1
    
    Lambda_cons = struc_form_constraint.(E_dens.(n,T),tff.(rho),n)
    if size(Lambda_cons,2)>3
        max_inds = argmax(Lambda_cons,dims=2)
        med = median(Lambda_cons,dims=2)
        med_inds = argmin(abs.(Lambda_cons.-med),dims=2)
        min_inds = argmin(Lambda_cons,dims=2)
        inds = [min_inds med_inds max_inds]
        return @views Lambda_cons[inds], direction, n[inds], T[inds]
    end
    
    return Lambda_cons,direction,n,T
end

function calcStrucFormConstraint(properties::Function,adm::ADM)
    return calcStrucFormConstraint(properties(),adm)
end

# We'll duplicate for now, but we can probably combine with checkMassLossConstraint
function checkStrucFormConstraint(Lambda,properties::Properties,adm::ADM;n_test=nothing)
    Lambda_cons, direction, n, T = calcStrucFormConstraint(properties,adm)
    if isnothing(n_test)
        Lambda_test = Lambda.(n,T,adm=adm)
    else
        Lambda_test = Lambda.(n_test,T,adm=adm)
    end
    check = direction==above ? Lambda_test.<=Lambda_cons : Lambda_test.>=Lambda_cons
    if properties.allow_any
        any_all = any
    else
        any_all = all
    end
    if ndims(check)==2
        acheck = []
        for i = 1:size(check)[2]
            push!(acheck,any_all(check[:,i]))
        end
        return acheck
    elseif ndims(check)>2
        throw(ErrorException("I don't know how to handle $(ndims(check)) dimensions"))
    else
        return any_all(check)
    end
end

function checkStrucFormConstraint(Lambda,properties::Function,adm::ADM;n_test=nothing)
    return checkStrucFormConstraint(Lambda,properties(),adm;n_test)
end

# Common function for constraint processing
# Input: Lambda(n,T,dark_params), constraint_type (0=mass-loss,1=structure-formation)
#        properties() - function that returns the necessary properties for the constraint
#            i.e., for a mass-loss, returns Mhalo, z (for n), Sigma_T, v_c, x_D (optional)
#        dark is a dark_params struct
@enum constraint_type massloss strucform
function checkConstraint(Lambda,ct::constraint_type,properties,adm::ADM;n_test=nothing)
    # Julia doesn't have a built-in switch statement and adding in the Match.jl package is
    # way overkill. Short-circuiting is ok, but switch is way better
    validate_ADM(adm)
    ct==massloss && return checkMassLossConstraint(Lambda,properties,adm,;n_test)
    ct==strucform && return checkStrucFormConstraint(Lambda,properties,adm;n_test)
    throw(ArgumentError("$ct is not a valid constraint-type"))
end

function checkAllConstraints(Lambda,adm::ADM;n_test=nothing)
    mlc = [getLSSProperties,getDwGalProperties,getDBCProperties,getMWProperties,getDCOProperties]
    names = ["LSS","DwG","DBC","MWG","DCO"]
    types = [massloss,massloss,massloss,strucform,strucform]
    checks = Dict()
    for i in eachindex(mlc)
        checks[names[i]] = checkConstraint(Lambda,types[i],mlc[i],adm;n_test)
    end
    return checks
end


function plotConstraint(properties,adm::ADM;fig=nothing,name=nothing,ctype::constraint_type=massloss)
    calcConstraint = ctype==massloss ? calcMassLossConstraint : calcStrucFormConstraint
    Lambda_cons, direction, n, T = calcConstraint(properties,adm)

    if isnothing(name)
        name="Constraint"
    end
    if size(Lambda_cons,1)==1
        marker="*"
        linestyle="none"
    else
        marker="none"
        linestyle=nothing
    end
    h = loglog(T, Lambda_cons;marker,linestyle,label=name)
    midT = mean(log10.([minimum(T),maximum(T)]))
    if size(T)[2]>1
        x = log10.(vec(T[:,1]))
    else
        x = log10.(vec(T))
    end
    y = vec(minimum(Lambda_cons,dims=2))
    if length(x)==1
        # we have a 1xn vector, just use everything
        x = log10.(vec(T))
        y = vec(Lambda_cons)
        #println("We have a 1x$(size(T)[2]) vector. Switching to using everything:\nx=$x\ny=$y")
    end
    if length(x)>1
        inds = sortperm(x)
        x = x[inds]
        y = y[inds]
        Interpolations.deduplicate_knots!(x,move_knots=true)
        Lambda_fit = linear_interpolation(x,y)
        y1 = Lambda_fit(midT)
    else
        y1 = y
    end
    if direction==above # allowed arrow points down
        # Needs to be done in logspace - makes this trickier
        y2 = y1/10
    else
        y2 = y1*10
    end
    col = h[1].get_color()
    PyPlot.annotate("",xytext=(10^midT,y1),xy=(10^midT,y2),arrowprops=Dict("arrowstyle"=>"->","facecolor"=>col,"edgecolor"=>col))
    #println("x=$midT\ny=$(Lambda_fit[midT])\nx+dx=$midT\ny+dy=$(y)")
    if size(h)[1]>1
        for h1 in h
            h1.update_from(h[1])
        end
    end
    if isnothing(fig)
        ylim(1e-35, 1e-18)
        xlim(1e0, 1e10)
        grid()
        return gcf,h
    end
    return fig,h
end


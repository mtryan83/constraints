module ADM_mod

export ADM, set_m!, set_M!, set_α!, validate_ADM

using Exceptions

mutable struct ADM
    # rm::Float64 = 1.0
    # rM::Float64 = 1.0
    # rα::Float64 = 1.0
    # p_mass::Float64 = 1.673e-24 # mass in grams
    # e_mass::Float64 = 9.109e-28
    # m_p::Float64 = 0.938        # mass in GeV
    # m_e::Float64 = 511          # mass in keV
    # α::Float64 = 7.2973525693e-3
    # ϵ::Float64 = 1.0
    # ξ::Float64 = 1.0
    # z::Float64 = 0.0
    # M_mass::Float64 = rM * p_mass # mass in grams
    # m_mass::Float64 = rm * e_mass
    # M::Float64 = rM * m_p         # mass in GeV
    # m::Float64 = rm * m_e         # mass in keV
    # dα::Float64 = rα * α
    p_mass::Float64 # mass in grams
    e_mass::Float64
    m_p::Float64    # mass in GeV
    m_e::Float64    # mass in keV
    α::Float64
    ϵ::Float64
    ξ::Float64
    z::Float64
    rm::Float64
    rM::Float64
    rα::Float64
    M_mass::Float64 # mass in grams
    m_mass::Float64
    M::Float64      # mass in GeV
    m::Float64      # mass in keV
    dα::Float64
    function ADM(;rm = nothing,rM = nothing,rα = nothing,p_mass = 1.673e-24,e_mass = 9.109e-28,m_p = 0.938,m_e = 511,α = 7.2973525693e-3,
            ϵ = 1.0,ξ = 1.0,z = 0.0,M_mass = nothing,m_mass = nothing,M = nothing,m = nothing,dα = nothing) 
        #adm = new(rm,rM,rα,p_mass,e_mass,m_p,m_e,α,ϵ,ξ,z,M_mass,m_mass,M,m,dα)
        adm = new(p_mass,e_mass,m_p,m_e,α,ϵ,ξ,z)
        set_m!(adm,rm=rm,m_mass=m_mass,m=m)
        set_M!(adm,rM=rM,M_mass=M_mass,M=M)
        set_α!(adm,rα=rα,dα=dα)
        validate_ADM(adm)
    end
end

EmptyParam=Union{Number,Some{Number},Nothing}
function set_m!(adm::ADM;rm::EmptyParam=nothing,m_mass::EmptyParam=nothing,m::EmptyParam=nothing)
    if isnothing(rm) && isnothing(m_mass) && isnothing(m) 
        adm.m = adm.m_e
        adm.m_mass = adm.e_mass
        adm.rm = 1
    elseif (isnothing(rm) ⊻ isnothing(m_mass) ⊻ isnothing(m)) ||
        (!isnothing(rm) && !isnothing(m_mass) && !isnothing(m))
        throw(ArgumentError("Must provide one and only one of rm($rm), m_mass($m_mass), or m($m)"))
    elseif !isnothing(rm)
        adm.rm = rm
        adm.m_mass = rm * adm.e_mass
        adm.m = rm * adm.m_e
    elseif !isnothing(m_mass)
        adm.m_mass = m_mass
        adm.rm = m_mass / adm.e_mass
        adm.m = adm.rm * adm.m_e
    else
        adm.m = m
        adm.rm = m / adm.m_e
        adm.m_mass = adm.rm * adm.e_mass
    end
    return adm
end

function set_M!(adm::ADM;rM=nothing,M_mass=nothing,M=nothing)
    if isnothing(rM) && isnothing(M_mass) && isnothing(M) 
        adm.M = adm.m_p
        adm.M_mass = adm.p_mass
        adm.rM = 1
    elseif (isnothing(rM) ⊻ isnothing(M_mass) ⊻ isnothing(M)) ||
        (!isnothing(rM) && !isnothing(M_mass) && !isnothing(M))
        throw(ArgumentError("Must provide one and only one of rM, M_mass, or M"))
    elseif !isnothing(rM)
        adm.rM = rM
        adm.M_mass = rM * adm.p_mass
        adm.M = rM * adm.m_p
    elseif !isnothing(M_mass)
        adm.M_mass = M_mass
        adm.rM = M_mass / adm.p_mass
        adm.M = adm.rM * adm.m_p
    else
        adm.M = M
        adm.rM = M / adm.m_p
        adm.M_mass = adm.rM * adm.p_mass
    end
    return adm
end

function set_α!(adm::ADM;rα=nothing,dα=nothing)
    if isnothing(rα) && isnothing(dα)
        adm.dα = adm.α
        adm.rα = 1
    elseif !(isnothing(rα) ⊻ isnothing(dα))
        throw(ArgumentError("Must provide one and only one of rα and dα"))
    elseif !isnothing(dα)
        adm.dα = dα
        adm.rα = dα/adm.α
    else
        adm.rα = rα
        adm.dα = rα * adm.α
    end
    return adm
end

# Check that all ADM parameters are consistent and pass basic
# constraints
@exception msgerror msg::String
@msgerror ADMException "ADM struct error" e.msg
function validate_ADM(adm::ADM)
    if !(adm.rm ≈ adm.m_mass/adm.e_mass && adm.rm ≈ adm.m/adm.m_e)
        throw(ADMException("e_d mismatch: rm=$(adm.rm) m_mass/e_mass=$(adm.m_mass/adm.e_mass) m/m_e=$(adm.m/adm.m_e)"))
    end
    if !(adm.rM ≈ adm.M_mass/adm.p_mass && adm.rM ≈ adm.M/adm.m_p)
        throw(ADMException("p_d mismatch: rM=$(adm.rM) M_mass/p_mass=$(adm.M_mass/adm.p_mass) M/m_p=$(adm.M/adm.m_p)"))
    end
    if !(0 <= adm.ϵ <= 1)
        throw(ADMException("ϵ outside of [0,1]!"))
    end
    if adm.z < 0
        throw(ADMException("z cannot be negative"))
    end
    if !(0 <= adm.ξ <= 1)
        throw(ADMException("ξ outside of [0,1]!"))
    end
    return adm
end

end
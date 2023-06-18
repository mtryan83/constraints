module DK_modded_files

export generate_react
export generate_test

using DelimitedFiles
using DataFrames
using PyPlot
using QuadGK
using LaTeXStrings
using Statistics
using Printf
using PyCall 

function generate_react(r_m,r_M,r_α, ξ)
    m_SM = 9.109e-28        # g
    M_SM = 1.673e-24        # g
    α_SM = 7.2973525693e-3  #
    M = r_M * M_SM
    m = r_m * m_SM
    α = r_α * α_SM
    xi=replace(@sprintf("%e", ξ ), "e"=>"d")
    qe_mass=replace(@sprintf("%e", m ), "e"=>"d")
    qp_mass=replace(@sprintf("%e", M ), "e"=>"d")
    Dalpha=replace(@sprintf("%e", α ), "e"=>"d")

    testfile= """
    @dark: qe_mass = $qe_mass !9.10938188d-28 !g
    @dark: qp_mass = $qp_mass !1.67262158d-24 !g
    @dark: Dalpha = $Dalpha !7.2973525664d-3 !dark fine structure constant
    @dark: xi = $xi
    # Cen92 uses 1.5
    @darkGauntFF: 1.5d0 ! Cen92 Gaunt free-free factor
    #@darkGauntFF: 1.1d0 + 0.34d0 * exp(-1.d0/3.d0*(5.5d0-log10(temp))**2.d0) ! Dark gaunt factor


    # Reaction network used in the DARK test
    # This network has the same chemical rates as react_dark_minimal, but also includes
    # ordinary matter rates that are turned off using the variable below. Toggling this
    # flag and adjusting the abundances in the dark/test.f90 file allows comparisons 
    # between ordinary matter and dark matter runs without significant changes to the
    # other files
    @var: rateOff = 0 ! Flag to turn unused rates on or off

    # Example custom cooling function
    #cooling block START
    #@cooling_start
    #
    # #additional variables in F90 format
    # # note that the @var tokens inside a @cooling block do not interfere with the other ones
    # @var:a=1.35d0-0.3d0*log10(Tgas*1d-3)
    # @var:lambda3=1.32d-23*(Tgas*1d-3)**a
    #
    # @var: x = log10(Tgas)
    # @var: logCool = 0.d0
    # if (x<7.5) then
    # 	logCool = -22.866077d0 + 0.625805d0*cos(x*1.339969) + 1.977851d0*sin(x*1.339969) + &
    #               1.404537d0*cos(2*x*1.339969d0) + 0.234920d0*sin(2*x*1.339969d0) + 0.091744d0*cos(3*x*1.339969d0) + -0.937203d0*sin(3*x*1.339969d0) + &
    #               -0.736857d0*cos(4*x*1.339969d0) + -0.044308d0*sin(4*x*1.339969d0) + -0.165526d0*cos(5*x*1.339969d0) + 0.387575d0*sin(5*x*1.339969d0) + &
    #               0.208907d0*cos(6*x*1.339969d0) + 0.079093d0*sin(6*x*1.339969d0) + 0.075171d0*cos(7*x*1.339969d0) + -0.009931d0*sin(7*x*1.339969d0) + &
    #               0.029750d0*cos(8*x*1.339969d0) + -0.002096d0*sin(8*x*1.339969d0)
    # else
    #	logCool = 0.361329d0*x + -25.469652d0
    # end if
    # @var: cooling = 10**(logCool)*n(idx_QH2)
    #
    #@cooling_stop
    #cooling block END
    # 

    #################################################################
    # Ordinary Sector
    #
    # Duplicate of react_primordial3
    # These rates should generally be turned off using the rateOff flag
    #
    #################################################################
    @var:T=Tgas

    #Janev 1987
    1,H,E,,H+,E,E,,NONE,NONE,rateOff*exp(-32.71396786d0+13.5365560d0*lnTe-5.73932875d0*(lnTe**2)+1.56315498d0*(lnTe**3)-0.28770560d0*(lnTe**4)+3.48255977d-2*(lnTe**5)-2.63197617d-3*(lnTe**6)+1.11954395d-4*(lnTe**7)-2.03914985d-6*(lnTe**8))

    #Abel et al. 1997, fit by data from Ferland et al. 1992
    2,H+,E,,H,g,,,NONE,.LE.5.5e3,rateOff*3.92d-13*invTe**0.6353d0
    3,H+,E,,H,g,,,>5.5e3,NONE,rateOff*exp(-28.61303380689232d0-0.7241125657826851d0*lnTe-0.02026044731984691d0*lnTe**2-0.002380861877349834d0*lnTe**3-0.0003212605213188796d0*lnTe**4-0.00001421502914054107d0*lnTe**5+4.989108920299513d-6*lnTe**6+5.755614137575758d-7*lnTe**7-1.856767039775261d-8*lnTe**8-3.071135243196595d-9*lnTe**9)

    #Janev 1987
    4,HE,E,,HE+,E,E,,NONE,NONE,rateOff*exp(-44.09864886d0+23.91596563d0*lnTe-10.7532302d0*(lnTe**2)+3.05803875d0*(lnTe**3)-0.56851189d0*(lnTe**4)+6.79539123d-2*(lnTe**5)-5.00905610d-3*(lnTe**6)+2.06723616d-4*(lnTe**7)-3.64916141d-6*(lnTe**8))

    #Radiative+dielectronic from Cen 1992, Aldrovandi & Pequignot 1973
    5,HE+,E,,HE,g,,,NONE,.LE.9.28e3,rateOff*3.92d-13*invTe**0.6353d0
    6,HE+,E,,HE,g,,,>9.28e3,NONE,rateOff*1.54d-9*(1.d0+0.3d0/exp(8.099328789667d0*invTe))/(exp(40.49664394833662d0*invTe)*Te**1.5d0)+3.92d-13/Te**0.6353d0

    #Aladdin database 1989 (see Abel et al. 1997)
    7,HE+,E,,HE++,E,E,,NONE,NONE,rateOff*exp(-68.71040990212001d0+43.93347632635d0*lnTe-18.48066993568d0*lnTe**2+4.701626486759002d0*lnTe**3-0.7692466334492d0*lnTe**4+0.08113042097303d0*lnTe**5-0.005324020628287001d0*lnTe**6+0.0001975705312221d0*lnTe**7-3.165581065665d-6*lnTe**8)

    # Verner & Ferland 1996 !NEW!
    8,HE++,E,,HE+,,,,NONE,NONE,rateOff*1.891d-10/(sqrt(Tgas/9.37)*(1.+sqrt(Tgas/9.37))**0.2476*(1.+sqrt(Tgas/2.774d6))**1.7524)

    # De Jong (1972) !NEW!
    9,H,E,,H-,,,,NONE,NONE,rateOff*1.4d-18*Tgas**0.928*exp(-Tgas/16200.)

    # fit by Kreckel et al. 2010 !NEW!
    @var:a1=1.3500e-09
    @var:a2=9.8493e-02
    @var:a3=3.2852e-01
    @var:a4=5.5610e-01
    @var:a5=2.7710e-07
    @var:a6=2.1826e+00
    @var:a7=6.1910e-03
    @var:a8=1.0461e+00
    @var:a9=8.9712e-11
    @var:a10=3.0424e+00
    @var:a11=3.2576e-14
    @var:a12=3.7741e+00
    10,H-,H,,H2,E,,,NONE,NONE,rateOff*a1*(Tgas**a2+a3*Tgas**a4+a5*Tgas**a6)/(1.+a7*Tgas**a8+a9*Tgas**a10+a11*Tgas**a12)

    # fit to Ramaker & Peek 1976, corrected by Coppola !NEW!
    11,H,H+,,H2+,,,,NONE,.LT.30.0d0,rateOff*2.10e-20*(Tgas/30.)**(-0.15)
    12,H,H+,,H2+,,,,.GE.30.0d0,NONE,rateOff*10**(-18.20-3.194*log10(Tgas)+1.786*(log10(Tgas))**2-0.2072*(log10(Tgas))**3)

    #Karpas 1979
    13,H2+,H,,H2,H+,,,NONE,NONE,rateOff*6.0d-10

    # fit by Savin et al. 2004, see also Glover et al. 2010  !NEW!
    @var:asav = 2.1237150d4
    @var:bsav1=-3.3232183d-7
    @var:bsav2= 3.3735382d-7
    @var:bsav3=-1.4491368d-7 
    @var:bsav4= 3.4172805d-8 
    @var:bsav5=-4.7813728d-9 
    @var:bsav6= 3.9731542d-10
    @var:bsav7=-1.8171411d-11
    @var:bsav8= 3.5311932d-13
    @var:sumsav=bsav1+bsav2*log(Tgas)+bsav3*(log(Tgas))**2+bsav4*(log(Tgas))**3+bsav5*(log(Tgas))**4+bsav6*(log(Tgas))**5+bsav7*(log(Tgas))**6+bsav8*(log(Tgas))**7
    14,H2,H+,,H2+,H,,,.GE.1.d2,.LE.3.d4,rateOff*sumsav*exp(-asav*invT)

    # Capitelli et al. 2007 ! NEW REACTION! 
    15,H2,E,,H,H-,,,NONE,NONE,rateOff*3.55d1*Tgas**(-2.28)*exp(-46707./Tgas)

    # fit by Mitchell & Deveau 1983 of data by Corrigan 1965 !NEW!
    16,H2,E,,H,H,E,,NONE,NONE,rateOff*4.38d-10*exp(-102000./Tgas)*Tgas**(0.35)

    #Martin et al 1996
    @noTabNext
    17,H2,H,,H,H,H,,NONE,NONE,rateOff*dissH2_Martin96(n, Tgas)

    #Janev 1987
    18,H-,E,,H,E,E,,NONE,NONE,rateOff*exp(-18.01849334273d0+2.360852208681d0*lnTe-0.2827443061704d0*lnTe**2+0.01623316639567d0*lnTe**3-0.03365012031362999d0*lnTe**4+0.01178329782711d0*lnTe**5-0.001656194699504d0*lnTe**6+0.0001068275202678d0*lnTe**7-2.631285809207d-6*lnTe**8)

    #Abel et al. 1997, based on Janev 1987
    19,H-,H,,H,H,E,,NONE,.LE.1.16e3,rateOff*2.56d-9*Te**1.78186d0
    20,H-,H,,H,H,E,,>1.16e3,NONE,rateOff*exp(-20.37260896533324d0+1.139449335841631d0*lnTe-0.1421013521554148d0*lnTe**2+0.00846445538663d0*lnTe**3-0.0014327641212992d0*lnTe**4+0.0002012250284791d0*lnTe**5+0.0000866396324309d0*lnTe**6-0.00002585009680264d0*lnTe**7+2.4555011970392d-6*lnTe**8-8.06838246118d-8*lnTe**9)

    # Stenrup et al. 2009 !NEW!
    21,H-,H+,,H,H,,,.GE.1e1,.LE.1e5,rateOff*(2.96d-6/sqrt(Tgas)-1.73d-9+2.50d-10*sqrt(Tgas)-7.77d-13*Tgas)

    #Poulart 1978
    22,H-,H+,,H2+,E,,,NONE,NONE,rateOff*1.d-8*Tgas**(-0.4d0)

    # fit by Coppola et al. (2011) !NEW!
    23,H2+,E,,H,H,,,NONE,.LE.1e4,rateOff*1.d6*(4.2278d-14-2.3088d-17*Tgas+7.3428d-21*Tgas**2-7.5474d-25*Tgas**3+3.3468d-29*Tgas**4-5.528d-34*Tgas**5)

    #Dalgarno & Lepp 1987
    24,H2+,H-,,H,H2,,,NONE,NONE,rateOff*5.d-7*sqrt(1.d2*invT)

    #Forrey 2013 !NEW! 
    25,H,H,H,H2,H,,,NONE,NONE,rateOff*(6.d-32*Tgas**(-0.25d0)+2.d-31*Tgas**(-0.5d0))

    #Glover&Abel 2008 
    26,H2,H,H,H2,H2,,,NONE,NONE,rateOff*(6.d-32*Tgas**(-0.25d0)+2.d-31*Tgas**(-0.5d0))/8.d0

    #Omukai 2001
    @var:invT = 1d0/Tgas
    @var:Hnuclei = get_Hnuclei(n(:))
    @var:kl21 = 1.18d-10*exp(-6.95d4*invT)
    @var:kh21 = 8.125d-8*T**(-0.5d0)*exp(-5.2d4*invT)*(1.d0-exp(-6d3*invT))
    @var:ncr21 = 1d1**(4.845d0-1.3d0*log10(T*1d-4)+1.62d0*(log10(T*1d-4))**2)
    @var:a21=1.d0/(1.d0+(Hnuclei/ncr21))
    @noTabNext
    27,H2,H2,,H,H,H2,,NONE,NONE,rateOff*kh21**(1.-a21)*kl21**a21


    #################################################################
    # Dark Sector
    #
    # Rescaled version of the ordinary sector, using the rescalings
    # from 2106.13245
    #
    #################################################################

    #@var:T=Tgas #Already defined
    @var:ra = Dalpha/fine_structure_constant
    @var:rP = qp_mass/p_mass
    @var:rE = qe_mass/e_mass
    @var:ra0 = 1.d0/(ra*rE)
    @var:rh = rE * ra**2
    @var:rrot = (ra*rE)**2/rP
    @var:rvib = ra**2 * rE**1.5d0 / (rP**5.d-1)

    @var:Th = T/rh
    @var:Tr = T/rrot
    @var:Tv = T/rvib

    #CMB Photon temperature
    @var:Tp = (1.d0+phys_zredshift)*xi*phys_Tcmb ! K
    @var:Tph = Tp/rh
    @var:Tpr = Tp/rrot
    @var:Tpv = Tp/rvib


    #Rosenberg and Fan 2017
    @var:DRyd=qe_mass*g_to_KeV*iboltzmann_eV*1.d3*Dalpha**2/2
    @var:y2=DRyd/T
    @var:logy2=log(y2)
    28,qH,qE,,qH+,qE,qE,,NONE,NONE,2.2d-7 * (e_mass/qe_mass)**(3.d0/2.d0) * sqrt(1.d5/T) * (exp(-y2)/2.d0+(y2*ei(-y2))/2)

    #Rosenberg and Fan 2017
    @var:Tcutoff=DRyd*16
    29,qH+,qE,,qH,qg,,,NONE,.LE.Tcutoff, 8.4d-14*(Dalpha/1.d-2)**3 * (e_mass/qe_mass)**(3.d0/2.d0) * sqrt(1.d5/T) * (1.744d0+logy2+1.d0/(6.d0*y2))
    30,qH+,qE,,qH,qg,,,>Tcutoff,NONE, 1.3d-15*(Dalpha/1.d-2)**5 * sqrt(e_mass/qe_mass) * (1.d6/T)**(3.d0/2.d0) * (-4.66d0-1.5d1*logy2+y2*(5.42d0-1.4d1*logy2))

    # De Jong (1972) !NEW! (Actually using the fit from Galli & Palla '98)
    31,qH,qE,,qH-,qg,,,NONE,NONE, (ra**2.d0)*(rE**(-2.d0)) * 1.4d-18*Th**0.928*exp(-Th/16200.)

    # fit by Kreckel et al. 2010 !NEW!
    #@var:a1=1.3500e-09 # Already defined
    #@var:a2=9.8493e-02
    #@var:a3=3.2852e-01
    #@var:a4=5.5610e-01
    #@var:a5=2.7710e-07
    #@var:a6=2.1826e+00
    #@var:a7=6.1910e-03
    #@var:a8=1.0461e+00
    #@var:a9=8.9712e-11
    #@var:a10=3.0424e+00
    #@var:a11=3.2576e-14
    #@var:a12=3.7741e+00
    32,qH-,qH,,qH2,qE,,,NONE,NONE, (ra**(-1.d0))*(rE**(-1.5d0))*(rP**(-5d-1)) *a1*(Th**a2+a3*Th**a4+a5*Th**a6)/(1.+a7*Th**a8+a9*Th**a10+a11*Th**a12)

    # fit to Ramaker & Peek 1976, corrected by Coppola !NEW!
    @var:scaleHHp = (ra**(2.d0))*(rE**(-1.0))*(rP**(-1.d0))
    33,qH,qH+,,qH2+,qg,,,NONE,.LT.30.0d0*rh,scaleHHp * 2.10e-20*(Th/30.)**(-0.15)
    34,qH,qH+,,qH2+,qg,,,.GE.30.0d0*rh,NONE,scaleHHp * 10**(-18.20-3.194*log10(Th)+1.786*(log10(Th))**2-0.2072*(log10(Th))**3)

    #Karpas 1979
    @var:chargeExchange = (ra**(-1.d0))*(rE**(-1.5d0))*(rP**(-0.5d0))
    35,qH2+,qH,,qH2,qH+,,,NONE,NONE,chargeExchange*6.0d-10

    # fit by Savin et al. 2004, see also Glover et al. 2010  !NEW!
    #@var:asav = 2.1237150d4 # Already defined
    #@var:bsav1=-3.3232183d-7
    #@var:bsav2= 3.3735382d-7
    #@var:bsav3=-1.4491368d-7 
    #@var:bsav4= 3.4172805d-8 
    #@var:bsav5=-4.7813728d-9 
    #@var:bsav6= 3.9731542d-10
    #@var:bsav7=-1.8171411d-11
    #@var:bsav8= 3.5311932d-13
    # Need to use new variable name, since we're using Th scaling here
    @var:sumsavD=bsav1+bsav2*log(Th)+bsav3*(log(Th))**2+bsav4*(log(Th))**3+bsav5*(log(Th))**4+bsav6*(log(Th))**5+bsav7*(log(Th))**6+bsav8*(log(Th))**7
    36,qH2,qH+,,qH2+,qH,,,.GE.1.d2*rh,.LE.3.d4*rh,chargeExchange*sumsavD*exp(-asav/Th)


    #Martin et al 1996
    #@noTabNext # Already defined
    37,qH2,qH,,qH,qH,qH,,NONE,NONE,dark_dissH2_Martin96(n, T)

    # Stenrup et al. 2009 !NEW!
    38,qH-,qH+,,qH,qH,,,.GE.1e1*rh,.LE.1e5*rh,(ra**(-3.d0))*(rE**(-3.d0))*(2.96d-6/sqrt(Th)-1.73d-9+2.50d-10*sqrt(Th)-7.77d-13*Th)

    # fit by Coppola et al. (2011) !NEW!
    # This rate has been left in so that the network passes the recombination
    # test, but since it hasn't been rescaled, it's set to 0
    39,qH2+,qE,,qH,qH,,,NONE,.LE.1e4,0 * 1.d6*(4.2278d-14-2.3088d-17*Tgas+7.3428d-21*Tgas**2-7.5474d-25*Tgas**3+3.3468d-29*Tgas**4-5.528d-34*Tgas**5)


    #Forrey 2013 !NEW! 
    @var:threeBody = (ra**(-4.d0))*(rE**(-4d0))*(rP**(-1d0))
    40,qH,qH,qH,qH2,qH,,,NONE,NONE,threeBody*(6.d-32*Th**(-0.25d0)+2.d-31*Th**(-0.5d0))

    #Glover&Abel 2008 
    41,qH2,qH,qH,qH2,qH2,,,NONE,NONE,threeBody*(6.d-32*Th**(-0.25d0)+2.d-31*Th**(-0.5d0))/8.d0

    #Yoshida 2006
    42,qH,qH,qH+,qH2,qH+,,,NONE,NONE,threeBody*(6.d-32*Th**(-0.25d0)+2.d-31*Th**(-0.5d0))

    #Yoshida 2006
    43,qH,qH,qH+,qH2+,qH,,,NONE,NONE,threeBody*(6.d-32*Th**(-0.25d0)+2.d-31*Th**(-0.5d0))/8.d0

    #Galli&Palla 1998
    44,qH2+,qH2,,qH3+,qH,,,NONE,NONE,(ra**(-1.d0))*(rE**(-1.5d0))*(rP**(-0.5d0))*2.d-9

    #Galli&Palla 1998
    45,qH3+,qE,,qH,qH2,,,NONE,NONE,(ra**(-1.d0))*(rE**(-2.d0))*(rP**(1.d0))*4.6d-6*Th**(-0.65d0)

    # Currently assuming transparent clouds
    # # Photo Reactions
    # @photo_start
    # @format:idx,R,R,P,P,Tmin,Tmax,rate
    # @var:ibeV=iboltzmann_eV
    # # Need to repeat these for photoreactions
    # @var:Tcutoff=DRyd
    # @var:ra = Dalpha/fine_structure_constant
    # @var:rP = qp_mass/p_mass
    # @var:rE = qe_mass/e_mass
    # @var:ra0 = 1.d0/(ra*rE)
    # @var:rh = rE * ra**2
    # @var:rrot = (ra*rE)**2/rP
    # @var:rvib = ra**2 * rE**1.5d0 / (rP**5.d-1)
    #
    # #Reaction #2 in Galli&Palla 98, inverse of 29/30
    # # Photoionization rate is computed in krome_subs using a (not great)
    # # method based on Rosenberg & Fan 2017 (great)
    # 46,qH,qg,qH+,qE,NONE,NONE,darkPhotoionizationRate()
    # 
    # #Reaction #4 in Galli&Palla 98, inverse of 36
    # 47,qH-,qg,qH,qE,NONE,NONE,(ra**5.d0)*(rE**(1.d0))*(rP**(0.d0)) * 1.1d-1 * (ibev*energy_eV/rh)**(2.13d0) * exp(-8823.d0/(ibev*energy_eV/rh))
    # 
    # #Reaction #9 in Galli&Palla 98, inverse of 38/39
    # 48,qH2+,qg,qH,qH+,NONE,NONE,(ra**5.d0)*(rE**(0.5d0))*(rP**(0.5d0)) * 2.e1*(ibev*energy_eV/rvib)**(1.59)*exp(-8.2d4/(ibev*energy_eV/rvib))
    # 
    # #Reaction #18 in Galli&Pall 98, photoionization of H2_dot
    # 49,qH2,qg,qH2+,qE,NONE,NONE,(rvib**2.44) * (rE**-5.d0) * (ra**-4.d0) * 2.9d2*(ibev*energy_eV/rvib)**1.56d0*exp(-178500/(ibev*energy_eV/rvib))
    # @photo_end
    """
    open("react_dark", "w") do io
        write(io, testfile)
    end
end
function generate_test(n,T,ξ,ϵ,z,r_m,r_M,r_α, r_rot_E, r_bin_E)
    xi = replace( @sprintf("%e", ξ), "e"=>"d")
    epsilon = replace( @sprintf("%e", ϵ ), "e"=>"d")
    n_s = @sprintf("%.2f",log10(n))
    T_s = @sprintf("%.2f",log10(T))
    rRotE = @sprintf("%.2f",log10(r_rot_E))
    rBinE = @sprintf("%.2f",log10(r_bin_E))
    T = replace( @sprintf("%e", T), "e"=>"d")
    n = replace( @sprintf("%e", n), "e"=>"d")
    zred = replace( @sprintf("%e", z), "e"=>"d")
    r_m = @sprintf("%.2f",log10(r_m))
    r_M = @sprintf("%.2f", log10(r_M))
    r_a = @sprintf("%.2f", log10(r_α))
    
    #name = "run_n" * n_s * "_T" * T_s * "_RE" * rRotE * "_BE" * rBinE 
                
    testfile = """
        program test_equillibrium_darkkrome

            use krome_main
            use krome_user
            use krome_user_commons
            use krome_cooling
            use krome_constants
            use starting_densities
            implicit none   

            ! Interface !
            interface
                subroutine my_dump_flux(n,T,nfile_in)
                  real*8, intent(in) :: n(*)
                  real*8, intent(in) :: T
                  integer,optional, intent(in) :: nfile_in
                end subroutine my_dump_flux
                subroutine open_out_file(nfile_in,fileName)
                  integer, intent(inout) :: nfile_in
                  character*40, intent(in) :: fileName
                  integer :: ios
                end subroutine open_out_file
            end interface

            ! Type Declerations !
            integer,parameter::rstep = 10000
            integer::i,j,unit,ios,numargs,dynamicDensity,relAbunError=0
            real*8::dtH,deldd
            real*8::tff,dd,dd1
            real*8::x(krome_nmols),Tgas,dt
            real*8::ntot,t_init,tuni
            real*8::rm,rMH,ra
            real*8::rho,rhocrit,rho_nddm,rhotot,rho_b
            real*8::Ttmp,h,time,Tvir,gamma,krome_gamma,tsFlag,mvircheck
            real*8::mjeans,tcool,mvir,R,tsound,tH,tlook,heat_array(krome_nheats)
            real*8::G,mp,me,alpha,kb,kbeV,cs,Told,xe,xHe,xH2,Msun,nmax,mu
            real*8::omega_m,omega_b,omega_dm,omega_ddm,omega_mz,epsilon,d,delta_c,delta_t,Mpcincm,Myrinsec,zred
            real*8::y2
            real*8::f_approx,ei,coll_ion
            real*8::rlarge,rsmall,k_one,x_e,t_rec,rec,CI,nion,t_lim,t_step,timee,dtt,c_sp,y2_lim
            real*8::rRotE, rBinE
            real*8::log10T,log10n

            character*40 arg
            character*40 fileName,datFileName





            ! Constants !
            G=6.6743d-08  ! G_Newton in cgs
            mp=p_mass  ! Proton mass in cgs
            me=9.10938188d-28  ! Electron mass in cgs
            c_sp=2.998e10
            alpha=7.3d-03
            kb=boltzmann_erg  ! k_B in cgs
            kbeV=boltzmann_ev
            h=0.6770  ! Reduced Hubble constant
            omega_m=0.3107
            omega_b=0 !.022447/h/h
            omega_dm=omega_m !*(1-omega_b/omega_m) ! amount of dark matter
            Mpcincm = 3.085678d24  ! one Megaparsec in cm
            delta_c = 18.d0*pi**2 !high redshift, no contribution from lambda
            delta_t = 9.d0*pi**2/16.d0! density contrast at turnaround
            Myrinsec = 60*60*24*365*1d6 
            Msun = 1.989d33  ! solar mass in g
            nmax = 1d22  ! Maximum number density we evolve to
            tH = 3.09e17/h ! Inverse Hubble in seconds

            rhocrit=(3 / (8*pi*G)) * (h*1.d7/Mpcincm)**2  !rho_crit today in g/cm^3

            ! Input custom variable definitions (i.e. from Julia) !
            Tgas = $T
            ntot = $n
            zred = $zred
            rm = $r_m
            rMH = $r_M
            ra = $r_a
            xi = $xi
            epsilon = $epsilon
            ! t_lim =         ! Must be in cgs !

            ! Variable definitions from command line arguments !
            numargs = iargc()
            if (numargs.GE.1) then 
                call getarg(1,arg)
            read(arg,*)Tgas
            !print *,"Running TEST with starting temp: ",Tgas
            if (numargs.GE.2) then
                call getarg(2,arg)
                read(arg,*)ntot
                    !print *, "Running TEST at density: ", ntot
                if (numargs.GE.3) then
                call getarg(3,arg)
                        read(arg,*)zred
                !print *,"Running TEST at redshift: ",zred
                if (numargs.GE.4) then
                            call getarg(4,arg)
                            read(arg,*)epsilon
                            !print *,"Running TEST at epsilon: ",epsilon
                    if (numargs.GE.5) then
                    !print *,"Error: wrong number of arguments: ",numargs
                    stop
                    end if
                endif
                endif
            endif
            endif

            ! Initial Conditions !
            omega_ddm=epsilon*omega_dm
            call krome_set_zredshift(zred)
            call krome_set_Tcmb(2.725d0)    ! This should not be multiplied by xi!
            call krome_set_Tfloor(2.725d0*xi)

            rhotot = omega_m * rhocrit*(1+krome_redshift)**3*delta_c !This is the total matter density

            rho_b = rhotot * omega_b/omega_m
            rho_nddm = rhotot * omega_dm/omega_m * (1-epsilon) ! Non-Dissipative Dark Matter density


            !initialize KROME (mandatory) !
            call krome_init()


            ! Use James' data to initialize abundances to compute mu and gamma !
            x(:) = get_n(1.d0,qe_mass*g_to_keV,qp_mass*g_to_keV/1.e6,Dalpha,xi,epsilon)
            mu = krome_get_mu(x) ! For ordinary matter this would be 1.22
            gamma = krome_get_gamma_x(x(:),Tgas)

            !print *,"Running with mu=",mu
            !print *,"Running with gamma=",gamma

            ! Thermally Averaged Atomic Binding Energy ! Rosenberg & Fan 2017
            y2 = ((10 ** rm) * 511000) * (((10 ** ra) * alpha)**2.d0) / (2.d0 * kbeV * Tgas) 

            ! Collisional Ionization Rate ! Rosenberg & Fan 2017
            CI =  2.2d-7 * (1/(10 ** rm))**(3.d0/2.d0) * sqrt(1.d5/Tgas) * (exp(-y2)/2.d0+(y2*ei(-y2))/2)

            ! Recombination Rate ! Rosenberg & Fan 2017
            rlarge = 8.4d-14 * (((10 ** ra) * alpha)/1.d-2)**3.d0 * (1.d0/ (10 ** rm) )**(3/2) * sqrt(1.d5 / Tgas) * ( 1.744d0 + log(y2) + 1.d0 / (6.d0 * y2) )
            rsmall = 1.3d-15 * (((10 ** ra) * alpha)/1.d-2)**5.d0 * sqrt(1.d0/(10 ** rm) ) * (1.d6 / Tgas)**(1.5d0) * ( -4.66 - 15d0 * log(y2) + y2 * (5.42d0 - 1.4d1 * log(y2) ) )
            y2_lim=16*Tgas*y2
            if (y2_lim.gt.0.25d0) then
              REC = rlarge
            else
              REC = rsmall
            end if

            ! Species default, cm-3 !
            x(:) = 0
            x(:) = get_n(ntot,qe_mass*g_to_keV,qp_mass*g_to_keV/1.e6,Dalpha,xi,epsilon)

            ! Species Abundances !
            x(KROME_idx_QE) = 1.d-6*ntot ! Dark e- 
            !x(KROME_idx_QE) = (CI/rec)*ntot ! Dark e- ! CIE Fraction
            x(KROME_idx_QHj) = x(KROME_idx_QE) ! Dark p ! Charge Neutrality
            x(KROME_idx_QH2) = 1.d-6*ntot ! Dark H2
            !x(KROME_idx_QH) = 1 ! Dark H ! Estimate
            x(KROME_idx_QH)  = (ntot-x(KROME_idx_QE) - 2.0d0 * x(KROME_idx_QH2)) ! Dark H ! Tegmark et al. 1997 Expression
            x(KROME_idx_QHk) = 0.0 ! Dark H-
            x(KROME_idx_QH2j) = 0.0 ! Dark H2+
            x(KROME_idx_QH3j) = 0.0 ! Darl H3+

            ! Recombination Timescale !
            t_rec = 1.d0/(1.0d-6*ntot * rec)
            t_lim = 2.d0 * t_rec

            ! Free-Fall Timescale !
            rho = krome_get_rho(x(:))
            tff = sqrt(3*pi/32/G/(rho+rho_nddm+rho_b))
            !t_lim = tff

            ! Sound-Crossing Timescale !
                !! mvir assumes a spherical, uniform density profile and gamma = 5/3. May try to change gamma value used !!
            mvir = (4 * pi * rhotot / 3.d0)**(-0.5) * (5 * kb * Tgas / (G*mp*krome_get_mu(x)))**1.5
            cs=sqrt(krome_get_gamma_x(x(:),Tgas)*kb*Tgas / (mp * krome_get_mu(x)))
            R=((mvir*omega_ddm/omega_m)/(4*pi*rho/3))**(1./3.)  ! Radius
            tsound = R/cs
            !t_lim = tsound

            ! Cooling Timescale !
            tcool = 3./2*kb*Tgas*sum(x)/cooling(x,Tgas)
            if(tcool<0) tcool=1.d38
            !t_lim = tcool

            ! Age of the Universe at zred !
            tuni = 2.05712d19/(100*h*sqrt(omega_m)*(1+zred)**(1.5))
            !t_lim = tuni

            ! Initializing Time and Timesteps !
            t_init = 0
            time = t_init
            dt = t_lim/rstep

            ! Loop on time steps !
            do i = 1,rstep

            ! Solve the Chemistry !
                call krome(x(:),Tgas,dt)


                ! Compute New Redshift !
                    !! If new redshift differs from old redshift by more than 5%, update
                    !! We do this because some redshift dependent calculations are slow and 
                    !! the results cached. So update those as infrequently as possible.
            zred = 7.50769d12/(time**2 * (100*h)**2 * omega_m)**(1.d0/3.d0) - 1.d0
                krome_redshift = krome_get_zredshift()
                if((zred.GE.0.d0).AND.((2*abs(zred-krome_redshift)/(zred+krome_redshift)).GE.0.015)) then
                    call krome_set_zredshift(zred)
                end if

                ! Update Time !
            time = time + dt

            ! Check if t_lim has been reached !
            if (time.GE.t_lim) then
                print '(E0.4)',cooling(x(:),Tgas)
                exit
                end if

            end do

        end program test_equillibrium_darkkrome   
        """
        open("tests/dark/test.f90", "w") do io
            write(io, testfile)
        end
        open("build/test.f90", "w") do io
            write(io, testfile)
        end
    end


end



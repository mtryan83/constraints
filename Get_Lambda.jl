module Get_Lambda

export get_Lambda

using DelimitedFiles
using DataFrames
using PyPlot
using QuadGK
using LaTeXStrings
using Statistics
using Printf
using PyCall

#### Class Structure ####

#### - aDM params
#### - reaction network
#### - atomic cooling
#### - molecular cooling (If there are problems with rescaling THIS IS WHERE IT WILL BE)
#### - heating
#### - master function

py"""



import numpy as np
from scipy.integrate import quad

HBAR = 1.054571817e-27  # erg s
H = 6.62607015e-27      # erg s
C = 2.9979245800e10     # cm/s
G_N = 6.6743e-8         # cm^3/g s^2
E = 1.602176634e-12  
K_B = 1.380649e-16      # erg/K
K_B_EV = 8.617333262e-5
m_SM = 9.109e-28        # g
M_SM = 1.673e-24        # g
a_SM = 7.2973525693e-3  #
m_SM_eV = 5.11e8

class aDM:
    def __init__(self, r_m, r_M, r_a, xi, epsilon, z):
        self.r_a = r_a
        self.r_m = r_m
        self.r_M = r_M
        self.xi = xi
        self.epsilon = epsilon
        self.z = z


class react_network:
    def __init__(self, n, r_m, r_M, r_a, xi, epsilon, z):
        ############################################
        ################ Parameters
        ############################################
        ADM = aDM(r_m, r_M, r_a, xi, epsilon, z)
        
        ## Dark Parameter Rescaling ##
        self.r_m = r_m
        self.r_M = r_M
        self.r_a = r_a
        self.xi = xi
        self.epsilon = epsilon

        ## Dark Parameter Values ##
        self.m = m_SM * ADM.r_m
        self.m_eV = m_SM_eV * ADM.r_m
        self.M = M_SM * ADM.r_M
        self.a = a_SM * ADM.r_a
        self.xi = ADM.xi
        self.epsilon = ADM.epsilon
        self.z = z
        
        ## Parameters Specific to the Halo ##
        self.n = n

        ## Thermally Avg. Atomic Binding ##
    def atom_bind(self, T):
        return (((self.m * (C ** 2)) * (self.a ** 2)) / (2 * K_B * T)) ** (1 / 2)
#### Rate Coefficients ####


    ## H + e -> p + e + e ## Rosenberg & Fan 2017 ##
    def k6(self,T):
        def argument(u):
        #    return ((u * (np.exp(-(u ** 2)))) / (1 + ((2 * (self.atom_bind(T) ** 2)) / (u ** 2)))) * (1 - ((self.atom_bind(T) ** 2) / (u ** 2)) - ((1 - ((self.atom_bind(T) ** 4) / (u ** 4))) / 2) * np.log((self.atom_bind(T) ** 2) / (u ** 2)) + ((((self.atom_bind(T) ** 2) / (u ** 2)) * np.log((self.atom_bind(T) ** 2) / (u ** 2))) / (1 + ((self.atom_bind(T) ** 2) / (u ** 2)))))
            return ((1 - ((self.atom_bind(T) ** 2)/(u ** 2))) * (u * (np.exp(-(u ** 2)))))
        return ((HBAR ** 2) * (C ** 0)) * ((((2 ** 7) * np.pi) / ((self.m ** 3) * (K_B * T))) ** (1/2)) * quad(argument, self.atom_bind(T), np.Infinity)[0]



    ## p + e -> H ## Rosenberg & Fan 2017 ##
    def k1(self, T):
        def argument(u):
            lst = np.logspace(0, 100, 50, base=10)
            def sum_argument(n):
               return ((u) * np.exp(- (u ** 2))) / (((u ** 2) * (n ** 3)) + ((self.atom_bind(T) ** 2) * n))
            return sum(list(map(sum_argument, lst)))
        return ((HBAR ** 2) * (C ** 2)) * (self.a ** 5) * ((((2 ** 11) * (np.pi)) / ((3 ** 3) * (self.m) * ((K_B * T) ** 3))) ** (1 / 2)) * quad(argument, 0, np.Infinity)[0]
        #y2 = ((511000 * self.r_m)*((a_SM * self.r_a) ** 2) / (2 * K_B_EV * T))
        #y2lim = 16*T*y2
        #return y2lim
        #if y2lim > 0.25:
        #    return 8.4e-14 * (((self.r_a * a_SM)/(1e-2)) ** 3) * ((self.r_m) ** -1.5) * ((1e5/T) ** 0.5) * (1.744 + np.log10(y2) + 1) 
        #else:
        #    return 1.3e-15 * (((self.r_a * a_SM)/1.e-2)**5) * ((1/(self.r_m))**0.5 ) * (1e6 / T)**(1.5) * ( -4.66 - 15 * np.log10(y2) + y2 * (5.42 - 14 * np.log10(y2) ) )



    ## H + e -> H^- ## Galli & Palla 1998 ##
    def k2(self, T):
        return ((self.r_a ** 2) * (self.r_m ** -2)) * 1.4e-18 * ((T/(self.r_m * (self.r_a ** 2))) ** 0.928) * np.exp(-T/(16200 * self.r_m * (self.r_a ** 2)))  
        
        ## H + e -> H^- ## Glover & Abel 2008 / Wishart 1979 ##
        # I don't know if this is scaled right #
        # def k_3(self, T):
        #    if T<=6000:
        #        return ((self.r_a ** 2) * ( r_m ** -2)) * ( 10 ** ( -17.845 + 0.762 * np.log10(T) + 0.1523 * (np.log10(T) ** 2) - 0.03274 * (np.log10(T) ** 3)))
        #    else:
        #        return ((self.r_a ** 2) * ( r_m ** -2)) * ( 10 ** ( -16.4199 + 0.1998 * (np.log10(T) ** 2) - 5.447e-3 * (np.log10(T) ** 4) + 4.0415e-5 * (np.log10(T) ** 6)))
        
    ## H^- + H -> H2 + e ## Kreckel et al. 2010 ##
    def k3(self, T):
        Th = T/(self.r_m * (self.r_a ** 2))
        a1 = 1.35e-9
        a2 = 9.8493e-2
        a3 = 3.2852e-1
        a4 = 5.5610e-1
        a5 = 2.771e-7
        a6 = 2.1826
        a7 = 6.191e-3
        a8 = 1.0461
        a9 = 8.9712e-11
        a10 = 3.0424
        a11 = 3.2576e-14
        a12 = 3.7741
        return (((self.r_a) * (self.r_m ** 1.5) * (self.r_M ** 0.5)) ** -1) * a1 * ((Th ** a2) + a3 * (Th ** a4) + a5 * (Th ** a6))/(1 + a7 * (Th ** a8) + a9 * (Th ** a10) + a11 * (Th ** a12))

    ## H + p -> H2^+ ## Ramaker & Peek 1976 / Coppola et al. 2011 ##
    def k4(self, T):
        Th = T/(self.r_m * (self.r_a ** 2))
        Tcrit =  30 * (self.r_m * (self.r_a ** 2))
        if T < Tcrit:
            return ((self.r_a ** 2)/(self.r_m * self.r_M)) * 2.1e-20 * ((Th/30) ** -0.15)
        elif T >= Tcrit:
            return ((self.r_a ** 2)/(self.r_m * self.r_M)) * 10 ** (-18.20 - 3.194 * np.log10(Th) + 1.786 * (np.log10(Th) ** 2) - 0.2072 * (np.log10(Th) ** 3))

    ## H2^+ + H -> H2 + p ## Galli & Palla 1998 ##
    def k5(self, T):
        return (6.0e-10)/((self.r_a) * (self.r_a ** 1.5) * (self.r_M ** 0.5))

    ## H2 + p -> H2^+ + H ## Savin et al. 2004 / Glover et al. 2010 ##
    def k7(self, T):
        Th = T/(self.r_m * (self.r_a ** 2))
        Tcrit1 = 100 * (self.r_m * (self.r_a ** 2))
        Tcrit2 = 30000 * (self.r_m * (self.r_a ** 2))
        a0 = 2.1237150e4
        a1 = -3.3232183e-7
        a2 = 3.3735382e-7
        a3 = -1.4491368e-7
        a4 = 3.4172805e-8
        a5 = -4.7813728e-9
        a6 = 3.9731542e-10
        a7 = -1.8171411e-11
        a8 = 3.5311932e-13
        if T >= Tcrit1 and T <= Tcrit2:
            return (1/((self.r_a) * (self.r_a ** 1.5) * (self.r_M ** 0.5))) * (np.exp(-a0/Th)) * (a1 + a2 * np.log(Th) + a3 * (np.log(Th) ** 2) + a4 * (np.log(Th) ** 3) + a5 * (np.log(Th) ** 4) + a6 * (np.log(Th) ** 5) + a7 * (np.log(Th) ** 6) + a8 * (np.log(Th) ** 7))
        else:
            return 0

    ## H2 + H -> 3H ## Martin et al. 1996 ##
    def k8(self, T, n_H):
        rrot = (self.r_a ** 2) * (self.r_m ** 2) * (self.r_M ** -1)
        rvib = (self.r_a ** 2) * (self.r_m ** 1.5) * (self.r_M ** -0.5)
        Th = np.minimum(T/(self.r_m * (self.r_a ** 2)), 4.5e4)
        iTh = 1/Th
        logTh = np.log10(Th)
        ilogTh = logTh ** -1
        logTh2 = logTh ** 2
        logTh3 = logTh ** 3
        Tv = Th * (self.r_m * (self.r_a ** 2)) * (rvib ** -1)
        iTv = Tv ** -1
        logTv = np.log10(Tv)
        ilogTv = logTv ** -1
        logTv2 = logTv ** 2
        logTv3 = logTv ** 3
        coeff_CD = [-178.4239, -68.42243, 43.20243, -4.633167, 69.70086, 40870.38, -23705.70, 128.8953, -53.91334, 5.315517, -19.73427, 16780.95, -25786.11, 14.82123, -4.890915, 0.4749030, -133.8283, -1.164408, 0.8227443, 0.5864073, -2.056313]
        coeff_DT = [-142.7664, 42.70741, -2.027365, -0.2582097, 21.36094, 27535.31, -21467.79, 60.34928, -27.43096, 2.676150, -11.28215, 14254.55, -23125.20, 9.305564, -2.464009, 0.1985955, 743.0600, -1.174242, 0.7502286, 0.2358848, 2.937507]
        n_c_rescale = ((self.r_a ** 8) * (self.r_m ** 4.75) * (self.r_M ** -1.75))
        
        ## Collisional Dissociative ##
        CD_kh1 = coeff_CD[0] + coeff_CD[1] * logTh + coeff_CD[2] * logTh2 + coeff_CD[3] * logTh3 + coeff_CD[4] * np.log10(1 + coeff_CD[5] * iTh)
        CD_kh2 = coeff_CD[6] * iTh
        CD_kl1 = coeff_CD[7] + coeff_CD[8] * logTh + coeff_CD[9] * logTh2 + coeff_CD[10] * np.log10(1 + coeff_CD[11] * iTh)
        CD_kl2 = coeff_CD[12] * iTh

        CD_nc1 = coeff_CD[13] + coeff_CD[14] * logTv + coeff_CD[15] * logTv2 + coeff_CD[16] * iTv
        CD_nc2 = CD_nc1 + coeff_CD[17]
        
        p = coeff_CD[18] + coeff_CD[19] * np.exp(-Th/1.85e3) + coeff_CD[20] * np.exp(-Th/4.4e2)

        CD_nc_1 = (10 ** CD_nc1 * (n_c_rescale) )
        CD_nc_2 = (10 ** CD_nc2 * (n_c_rescale) )

        k_CD = CD_kh1 - ((CD_kh1 - CD_kl1)/(1 + ((n_H / CD_nc_1) ** p))) + CD_kh2 - ((CD_kh2 - CD_kl2)/(1 + ((n_H / CD_nc_2) ** p)))


        ## Dissociative Tunneling ##
        DT_kh1 = coeff_DT[0] + coeff_DT[1] * logTh + coeff_DT[2] * logTh2 + coeff_DT[3] * logTh3 + coeff_DT[4] * np.log10(1 + coeff_DT[5] * iTh)
        DT_kh2 = coeff_DT[6] * iTh
        DT_kl1 = coeff_DT[7] + coeff_DT[8] * logTh + coeff_DT[9] * logTh2 + coeff_DT[10] * np.log10(1 + coeff_DT[11] * iTh)
        DT_kl2 = coeff_DT[12] * iTh

        DT_nc1 = coeff_DT[13] + coeff_DT[14] * logTv + coeff_DT[15] * logTv2 + coeff_DT[16] * iTv
        DT_nc2 = DT_nc1 + coeff_DT[17]
        
        p = coeff_DT[18] + coeff_DT[19] * np.exp(-Th/1.85e3) + coeff_DT[20] * np.exp(-Th/4.4e2)

        DT_nc_1 = (10 ** DT_nc1 * (n_c_rescale) )
        DT_nc_2 = (10 ** DT_nc2 * (n_c_rescale) )

        k_DT = DT_kh1 - ((DT_kh1 - DT_kl1)/(1 + ((n_H / DT_nc_1) ** p))) + DT_kh2 - ((DT_kh2 - DT_kl2)/(1 + ((n_H / DT_nc_2) ** p)))


        k8 = (10**k_CD) + (10**k_DT)

        if np.isfinite(k8) == True and np.isnan(k8) != True:
            return ((self.r_a ** -1) * (self.r_m ** -1.5) * (self.r_M ** -0.5)) * k8
        else:
            return 0

    ## H^- + p -> 2H ## Stenrup et al. 2009 ##
    def k9(self, T):
        Th = T/((self.r_m) * (self.r_a ** 2))
        Tcrit1 = 10 * ((self.r_m) * (self.r_a ** 2))
        Tcrit2 = 1e5 * ((self.r_m) * (self.r_a ** 2))
        
        if T <= Tcrit2 and T >= Tcrit1:
            return ((self.r_a ** -3) * (self.r_m ** -3)) * (((2.96e-6)/(Th ** 0.5)) - 1.73e-9 + (2.5e-10 * (Th ** 0.5)) - 7.77e-13 * Th)
        else:
            return 0
    ## H2^+ + e -> 2H ## Coppola et al. 2011 ##
    def k10(self, T):
        #if T <= 1e4:
        #    return 1e6 * (4.2278e-14 - (2.3088e-17 * T) + (7.342e-21 * (T ** 2)) - (7.5474e-25 * (T ** 3)) + (3.3468e-29 * (T ** 4)) - (5.528e-34  * (T ** 5)))
        #else:
        #    return 0
        # Set to 0 because we don't know rescaling #
        return 0

    ## 3H -> H2 + H ## Forrey 2013 ##
    def k11(self, T):
        Th = T/((self.r_m) * (self.r_a ** 2))
        return ((self.r_a ** -4) * (self.r_m ** -4) * (self.r_M ** -1)) * ((6e-32 * (Th ** -0.25)) + (2e-31 * (Th ** -0.5)))

    ## 2H + H2 -> H2 + H2 ## Glover & Abel 2008 ##
    def k12(self, T):
        return self.k11(T)/8

    ## 2H + H^+ -> H2 + H^+ ## Yoshida 2006 ##
    def k13(self, T):
        return self.k11(T)

    ## 2H + H^+ -> H2^+ + H ## Yoshida 2006 ##
    def k14(self, T):
        return self.k11(T)/8

    ## H2^+ + H2 -> H3^+ + H ## Galli & Palla 1998 ##
    def k15(self, T):
        return ((self.r_a ** -1) * (self.r_m ** -1.5) * (self.r_M ** -0.5)) * 2e-9

    ## H3^+ + e -> H + H2 ## Galli & Palla 1998 ##
    def k16(self, T):
        Th = T/((self.r_m) * (self.r_a ** 2))
        return ((self.r_a ** -1) * (self.r_m ** -1.5) * (self.r_M ** -0.5)) * 4.6e-6 * (Th ** -0.65)


class abundances:
    def __init__(self, n, T, x_H2_pre, r_m, r_M, r_a, xi, epsilon, z):
        
        ############################################
        ################ Parameters
        ############################################
        ADM = aDM(r_m, r_M, r_a, xi, epsilon, z)
        RN = react_network(n, r_m, r_M, r_a, xi, epsilon, z)

        
        ## Dark Parameter Rescaling ##
        self.r_m = r_m
        self.r_M = r_M
        self.r_a = r_a
        self.xi = xi
        self.epsilon = epsilon

        ## Dark Parameter Values ##
        self.m = m_SM * ADM.r_m
        self.m_eV = m_SM_eV * ADM.r_m
        self.M = M_SM * ADM.r_M
        self.a = a_SM * ADM.r_a
        self.xi = ADM.xi
        self.epsilon = ADM.epsilon
        self.z = z
        
        ## Parameters Specific to the Halo ##
        self.n = n
        self.T = T

        ## Rates ##
        self.Rec = RN.k1(self.T)
        self.CI = RN.k6(self.T)
        self.k3 = RN.k3(self.T)
        self.k4 = RN.k4(self.T)
        self.k5 = RN.k5(self.T)
        self.k7 = RN.k7(self.T)
        self.k9 = RN.k9(self.T)
        self.k10 = RN.k10(self.T)
        self.k14 = RN.k14(self.T)
        self.k15 = RN.k15(self.T)
        self.k16 = RN.k16(self.T)


        self.ionfrac = self.CI / self.Rec
        
        self.x_e = self.ionfrac/(1 + self.ionfrac)
        self.x_p = self.x_e
        self.x_H = 1/(1+self.ionfrac)
        self.x_H2 = x_H2_pre * 0.5 * np.exp((- K_B_EV * self.T)/((self.r_a ** 2) * (self.r_m) * 4480))
        self.x_Hneg = (self.x_e * self.x_p * self.Rec)/(self.x_H * self.k3 + self.x_p * self.k9)
        self.x_H2pl = (self.x_H * self.x_p * self.k4 + self.x_H2 * self.x_p * self.k7 + self.x_H * self.x_H * self.x_p * self.n * self.k14)/(self.x_H * self.k5 + self.x_e * self.k10 + self.x_H2 * self.k15)
        #self.x_H3pl_temp = (self.x_H2pl * self.x_H2 * self.k15)/(self.x_e * self.k16)
        if self.x_H2pl * self.x_H2 < 1e-120:
            self.x_H3pl = 0
        else:
            self.x_H3pl = (self.x_H2pl * self.x_H2 * self.k15)/(self.x_e * self.k16)
            
        
    


class atomic:
    def __init__(self, n, T, x_H2_pre, r_m, r_M, r_a, xi, epsilon, z):



############################################
################ Parameters
############################################
        ADM = aDM(r_m, r_M, r_a, xi, epsilon, z)
        abund = abundances(n, T, x_H2_pre, r_m, r_M, r_a, xi, epsilon, z)
        
        ## Dark Parameter Values ##
        self.m = m_SM * ADM.r_m
        self.m_eV = m_SM_eV * ADM.r_m
        self.M = M_SM * ADM.r_M
        self.a = a_SM * ADM.r_a
        self.xi = ADM.xi
        self.epsilon = ADM.epsilon
        self.z = z

        ## Dark Parameter Rescaling ##
        self.r_m = r_m
        self.r_M = r_M
        self.r_a = r_a

        ## Parameters Specific to the Halo ##
        self.n = n

        ## Population Fractions ##
        self.x_e = abund.x_e
        self.x_H2 = abund.x_H2
        self.x_H = abund.x_H
        self.x_p = abund.x_p

        ## Populations ##
        self.n_H2 = self.x_H2 * self.n
        self.n_H = self.x_H * self.n
        self.n_p = self.x_p * self.n
        self.n_e = self.x_e * self.n

        ## Background Temp ##
        self.T_D_gamma = (1 + self.z) * self.xi * 2.7

        ## Thermally Avg. Atomic Binding ##
    def atom_bind(self, T):
        return (((self.m * (C ** 2)) * (self.a ** 2)) / (2 * K_B * T)) ** (1 / 2)




############################################
############################################



############################################
################ Atomic Cooling Rates
################ From Rosenberg and Fan 2017
############################################

    def compton_scat(self, T):

        return ((K_B) ** (5)) * ((HBAR) ** (-1)) * ((C) ** -6) * (4 * (T-self.T_D_gamma)/self.m) * ((8 * np.pi)/3) * (((self.a)/(self.m)) ** 2) * ((np.pi ** 2)/15) * (self.T_D_gamma ** 4)



    def atomic_collisional_ion(self, T):
        def argument(u):
            return ((u * (np.exp(-(u ** 2)))) / (1 + ((2 * (self.atom_bind(T) ** 2)) / (u ** 2)))) * (1 - ((self.atom_bind(T) ** 2) / (u ** 2)) - ((1 - ((self.atom_bind(T) ** 4) / (u ** 4))) / 2) * np.log((self.atom_bind(T) ** 2) / (u ** 2)) + ((((self.atom_bind(T) ** 2) / (u ** 2)) * np.log((self.atom_bind(T) ** 2) / (u ** 2))) / (1 + ((self.atom_bind(T) ** 2) / (u ** 2)))))
            return ((1 - ((self.atom_bind(T) ** 2)/(u ** 2))) * (u * (np.exp(-(u ** 2)))))
        return ((HBAR ** 2) * (C ** 2)) * (self.a ** 2) * (self.m) * ((((2 ** 5) * np.pi) / ((self.m ** 3) * (K_B * T))) ** (1/2)) * quad(argument, self.atom_bind(T), np.Infinity)[0]


    def recombination(self, T):
        def argument(u):
            lst = np.logspace(0, 100, 50, base=10)

            def sum_argument(n):
                return ((u ** 3) * np.exp(- (u ** 2))) / (((u ** 2) * (n ** 3)) + ((self.atom_bind(T) ** 2) * n))

            return sum(list(map(sum_argument, lst)))

        return ((HBAR ** 2) * (C ** 2)) * (self.a ** 5) * ((((2 ** 11) * (np.pi)) / ((3 ** 3) * (self.m) * (K_B * T))) ** (1 / 2)) * quad(argument, 0, np.Infinity)[0]


    def bremstrahlung(self, T):
        return (HBAR ** 2) * (16/3) * (((2 * np.pi)/(3)) ** (1/2)) * ((self.a ** 3)/(self.m ** 2)) * ((K_B * T) ** (1/2)) * (self.m ** (1/2))


    def atomic_collisional_excitation(self, T):
        lowerbound = ((3 ** (1/2))/(2)) * self.atom_bind(T)
        def argument(u):
            return (np.log((4 * u)/(self.atom_bind(T)))) * ((u * np.exp(-(u ** 2)))/(1 + ((7 * (self.atom_bind(T) ** (2))/ (4 * (u ** 2))))))

        return (HBAR ** 2) * (C ** 3) * (((2 ** 16)/(3 ** 9)) * ((2 * np.pi)/((self.m * (C ** 2)) * K_B * T)) ** (1/2)) * ((self.a) ** 2) * quad(argument, lowerbound, np.Infinity)[0]

############################################
############################################



############################################
################ Total Atomic Cooling Rate
############################################

    def total_dark_atomic_cooling(self, n, T):  ## ergs/cm*s ##
        return (self.x_H * self.x_e * self.atomic_collisional_ion(T)) + (self.x_e * self.x_p * self.recombination(T)) + (self.x_e * self.x_p * self.bremstrahlung(T)) + (self.x_H * self.x_e * self.atomic_collisional_excitation(T)) + (self.x_e * self.compton_scat(T))
############################################
############################################




class molecular:
    def __init__(self, n, T, x_H2_pre , r_m, r_M, r_a, xi, epsilon, z):




############################################
################ Parameters
############################################
        ADM = aDM(r_m, r_M, r_a, xi, epsilon, z)
        abund = abundances(n, T, x_H2_pre, r_m, r_M, r_a, xi, epsilon, z)
        RN = react_network(n, r_m, r_M, r_a, xi, epsilon, z)


        ## Dark Parameter Values ##
        self.m = m_SM * ADM.r_m
        self.m_eV = m_SM_eV * ADM.r_m
        self.M = M_SM * ADM.r_M
        self.a = a_SM * ADM.r_a
        self.xi = ADM.xi
        self.epsilon = ADM.epsilon
        self.z = z

        ## Dark Parameter Rescaling ##
        self.r_m = r_m
        self.r_M = r_M
        self.r_a = r_a

        ## Parameters Specific to the Halo ##
        self.n = n
        self.T = T

        ## Population Fractions ##
        self.x_H2 = abund.x_H2
        self.x_Hneg = abund.x_Hneg
        self.x_H2pl = abund.x_H2pl
        self.x_H3pl = abund.x_H3pl
        self.x_e = abund.x_e
        self.x_p = abund.x_p
        self.x_H = abund.x_H


        ## Populations ##
        self.n_H2 = self.x_H2 * self.n
        self.n_H2pl = self.x_H2pl * self.n
        self.n_Hneg = self.x_Hneg * self.n
        self.n_H3pl = self.x_H3pl * self.n
        self.n_H = self.x_H * self.n
        self.n_p = self.x_p * self.n
        self.n_e = self.x_e * self.n
    
        ## Reactions ##
        self.k4 = RN.k4(self.T)
        self.k8 = RN.k8(self.T, self.n_H)


        ## Coefficients ##
        ## From Table 8 in Glover and Abel 2008 and Appendix A of Glover 2015##
        self.a_i_H_10_100 = [-16.818342, 37.383713, 58.145166, 48.656103, 20.159831, 3.847961, 0,0,0]
        self.a_i_H_100_1000 = [-24.311209, 3.5692468, -11.332860, -27.850082, -21.328264, -4.2519023, 0,0,0]
        self.a_i_H_1000_6000 = [-24.311209, 4.6450521, -3.7209846, 5.9369081, -5.5108047, 1.5538288, 0,0,0]
        self.a_i_p_10_10000 = [-22.089523,  1.5714711, 0.015391166, -0.23619985, -0.51002221, 0.32168730, 0,0,0]
        self.a_i_e_10_200 = [-21.928796,  16.815730,96.743155,  343.19180, 734.71651,  983.67576, 801.81247, 364.14446,  70.609154]
        self.a_i_e_200_10000 = [-22.921189,  1.6802758, 0.93310622, 4.0406627, -4.7274036, -8.8077017, 8.9167183, 6.4380698, -6.3701156]
        self.a_i_h2_100_6000 = [-23.962112, 2.09433740, -0.77151436, 0.43693353, -0.14913216, -0.033638326, 0,0,0]
        self.a_i_HDL = [-20.584225, 5.0194035, -1.5739905, -4.7155769, 2.4714161,5.4710750, -3.9467356, -2.2148338, 1.8161874]

        


        ## Limiting Temperatures for Coefficients##
        self.T_H_lim_1 = 100
        self.T_H_lim_2 = 1000
        self.T_H_lim_3 = 6000
        self.T_p_lim = 1000000000000
        self.T_e_lim_1 = 500
        self.T_e_lim_2 = 100000000000
        self.T_H2_lim = 6000
        self.global_lower = 10

############################################
############################################





############################################
######## Low Density Molecular Cooling
############################################



############################################
######## Fitting Function
######## From Eq(40) in Glover and Abel 2008
############################################

    def fitting_function(self, coeff, T, type):


        if type == 1:

            T3 = T / 1000
            logt3 = np.log10(T3)
            logt32 = logt3 ** 2
            logt33 = logt3 ** 3
            logt34 = logt3 ** 4
            logt35 = logt3 ** 5
            logt36 = logt3 ** 6
            logt37 = logt3 ** 7
            logt38 = logt3 ** 8
            a1 = coeff[0]
            a2 = coeff[1]
            a3 = coeff[2]
            a4 = coeff[3]
            a5 = coeff[4]
            a6 = coeff[5]
            a7= coeff[6]
            a8 = coeff[7]
            a9 = coeff[8]
            return ((a1) + (a2 * logt3) + (a3 * logt32) + (a4 * logt33) + (a5 * logt34) + (a6 * logt35) + (a7 * logt36) + (a8 * logt37) + (a9 * logt38))
        else:
            T3 = T / 1000
            logt3 = np.log10(T3)
            logt32 = logt3 ** 2
            logt33 = logt3 ** 3
            logt34 = logt3 ** 4
            logt35 = logt3 ** 5
            a1 = coeff[0]
            a2 = coeff[1]
            a3 = coeff[2]
            a4 = coeff[3]
            a5 = coeff[4]
            a6 = coeff[5]

            curve = (a1) + (a2 * logt3) + (a3 * logt32) + (a4 * logt33) + (a5 * logt34) + (a6 * logt35)
            return curve
############################################
############################################


    def sigmoid(self, x, x_0, s):
        return 10/(10+np.exp(-s*(x-x_0)))

    def wCool(self, logTgas, logTmin, logTmax):
        x=(logTgas-logTmin)/(logTmax-logTmin)
        wcool = 10 ** (200 * (self.sigmoid(x, -0.2, 50)*self.sigmoid(-x, -1.2, 50) - 1))
        if wcool < 1e-199:
            return 0
        else:
            return wcool


############################################
######## SM Molecular Cooling Rates
######## Glover and Abel 2008
############################################

    def H_H2_coll(self, T, winL, winH):
        if T <= self.T_H_lim_1:
            coefficient = self.a_i_H_10_100
            return 10 ** self.fitting_function(coefficient, T, 0)
        elif T > self.T_H_lim_1 and T <= self.T_H_lim_2:
            coefficient = self.a_i_H_100_1000
            return 10 ** self.fitting_function(coefficient, T, 0)
        elif T > self.T_H_lim_2 and T <= self.T_H_lim_3:
            coefficient = self.a_i_H_1000_6000
            return 10 ** self.fitting_function(coefficient, T, 0)
        else:
            return 1.862314467912518e-22 * self.wCool(np.log10(T), 1 + winL, np.log10(6000) + winH)


    def H2_H2_coll(self, T, winL, winH):
        #if T <= self.T_H2_lim:
        coefficient = self.a_i_h2_100_6000
        H2_H2 = (10 ** self.fitting_function(coefficient, T, 0)) * self.wCool(np.log10(T), 2 + winL, 4 + winH)
        if np.isfinite(H2_H2) == True:
            return H2_H2
        else:
            return 0
        #elif T> self.T_H2_lim:
        #    return 0

    def e_H2_coll(self, T, winL, winH):
        if T <= self.T_e_lim_1:
            coefficient = self.a_i_e_10_200
            return (10 ** min((self.fitting_function(coefficient, T, 1)), 30)) * self.wCool(np.log10(T), 2 + winL, 4 + winH)
        elif T > self.T_e_lim_1: # and T <= self.T_e_lim_2:
            coefficient = self.a_i_e_200_10000
            return (10 ** self.fitting_function(coefficient, T, 1)) * self.wCool(np.log10(T), 2 + winL, 4 + winH)
        #elif T > self.T_e_lim_2:
        #    return 0

    def p_H2_coll(self, T, winL, winH):
        if T <= 1e4:
            coefficient = self.a_i_p_10_10000
            return 10 ** self.fitting_function(coefficient, T, 0)
        else:
            return (1.182509139382060e-21) * self.wCool(np.log10(T), 1 + winL, 4 + winH)

############################################
############################################

############################################
######## DM Molecular Cooling Rates
######## Ryan et al 2022, Glover and Abel 2008
############################################

    def H2_H2_DM_coll(self, T):
        TscaleR = (self.r_a ** -2) * (self.r_m ** -2) * (self.r_M ** 1)
        TscaleV = (self.r_a ** -2) * (self.r_m ** -1.5) * (self.r_M ** 0.5)
        TscaleA = (self.r_a ** -2) * (self.r_m ** -1) * (self.r_M ** 0)


        ## The temperature at which the SM rotational rate equals the SM vibrational rate ##
        T_0_SM = 5.40244e3

        ## The temperature at which the DM rotational rate equals the DM vibrational rate ##
        T_0_DM = (((self.r_a) ** 2) * ((self.r_m) ** 1.5882) * ((self.r_M) ** -0.5861)) * T_0_SM

        ## The limitations of the DM rotational and vibrational rates ##
                    ## From Appendix D in Ryan et al. 2022
        T_0_r = (((self.r_a) ** 2) * ((self.r_m) ** 2) * ((self.r_M) ** -1)) * T_0_DM
        T_0_v = (((self.r_a) ** 2) * ((self.r_m) ** 1.5) * ((self.r_M) ** -0.5)) * T_0_DM

        ## The rescaled temperature values ##
                    ## Ryan et al. 2022
        T_r = TscaleR * T
        T_v = TscaleV * T

        ## The rate at the rescaled temperature values ##
        H2_H2_DM_coll_r = self.H2_H2_coll(T_r, 0, np.log10((TscaleR/TscaleA)))
        H2_H2_DM_coll_v = self.H2_H2_coll(T_v, np.log10((TscaleV/TscaleR)), np.log10((TscaleV/TscaleA)))


        #### The following conditionals are from DarkKROME, Ryan et al. 2022 ####
                ## (The computeMultiCase function in the krome_cooling file) ##



        ### There are 5 cases ###
        ## Case 1 ## r_M <= r_m -> then the rovib rate is just the sum of the rescaled rotational and the rescaled vibrational rates
        ## Case 2 ## r_M > r_m && T<=T_0_r  -> then the rovib rate is just the rescaled rotational rate
        ## Case 3 ## r_M > r_m && T_0_r < T < T_0_DM  -> then the rovib cooling is the linearly extrapolated rescaled rotational rate from T_0_r
        ## Case 4 ## r_M > r_m && T_0_DM < T < T_0_v  -> then the rovib cooling is the linearly extrapolated rescaled vibrational rate from T_0_v
        ## Case 5 ## r_M > r_m && T >= T_0_v  -> then the rovib cooling is just the rescaled vibrational rate




        ## Case 1 ##
        if self.r_M <= self.r_m:
            if T <= T_0_r and T < T_0_v:
                return (((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H2_H2_DM_coll_r)

            elif T >= T_0_v and T > T_0_r:
                return (((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -0.75)) * H2_H2_DM_coll_v)

            elif T <= T_0_r and T >= T_0_v:
                return (((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H2_H2_DM_coll_r) + (((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -1.25)) * H2_H2_DM_coll_v)

        ## Case 2 ##
        elif self.r_M > self.r_m and T <= T_0_r:
            return ((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H2_H2_DM_coll_r

        ## Case 3 ##
        elif self.r_M > self.r_m and T > T_0_r and T <= T_0_DM:
            x_2 = T_0_r
            x_1 = 0.98 * x_2
            y_2 = self.H2_H2_coll(x_2, 0,0)
            y_1 = self.H2_H2_coll(x_1, 0,0)
            slope = (y_2 - y_1)/(x_2 - x_1)
            H2_H2_lin = slope * (T_r - T_0_r) + y_2
            return ((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H2_H2_lin

        ## Case 4 ##
        elif self.r_M > self.r_m and T > T_0_DM and T <= T_0_v:
            x_2 = T_0_v
            x_1 = 1.01 * x_2
            y_2 = self.H2_H2_coll(x_2, 0,0)
            y_1 = self.H2_H2_coll(x_1, 0,0)
            slope = (y_2 - y_1)/(x_2 - x_1)
            H2_H2_lin = slope * (T_v - T_0_v) + y_2
            return ((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -1.25)) * H2_H2_lin

        ## Case 5 ##
        elif self.r_M > self.r_m and T > T_0_v:
            return ((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -1.25)) * H2_H2_DM_coll_v





    def H_H2_DM_coll(self, T):
        TscaleR = (self.r_a ** -2) * (self.r_m ** -2) * (self.r_M ** 1)
        TscaleV = (self.r_a ** -2) * (self.r_m ** -1.5) * (self.r_M ** 0.5)
        TscaleA = (self.r_a ** -2) * (self.r_m ** -1) * (self.r_M ** 0)

        ## The temperature at which the SM rotational rate equals the SM vibrational rate ##
        T_0_SM = 855.833

        ## The temperature at which the DM rotational rate equals the DM vibrational rate ##
        T_0_DM = (((self.r_a) ** 2) * ((self.r_m) ** 1.5) * ((self.r_M) ** -0.5)) * T_0_SM

        ## The limitations of the DM rotational and vibrational rates ##
                    ## From Appendix D in Ryan et al. 2022
        T_0_r = (((self.r_a) ** 2) * ((self.r_m) ** 2) * ((self.r_M) ** -1)) * T_0_DM
        T_0_v = (((self.r_a) ** 2) * ((self.r_m) ** 1.5) * ((self.r_M) ** -0.5)) * T_0_DM

        ## The rescaled temperature values ##
                    ## Ryan et al. 2022
        T_r = TscaleR * T
        T_v = TscaleV * T

        ## The rate at the rescaled temperature values ##
        H_H2_DM_coll_r = self.H_H2_coll(T_r, 0, np.log10((TscaleR/TscaleA)))
        H_H2_DM_coll_v = self.H_H2_coll(T_v, np.log10((TscaleV/TscaleR)), np.log10((TscaleV/TscaleA)))




        #### The following conditionals are from DarkKROME, Ryan et al. 2022 ####
                    ## (The computeMultiCase function in the krome_cooling file) ##


        ### There are 5 cases ###
        ## Case 1 ## r_M <= r_m -> then the rovib rate is just the sum of the rescaled rotational and the rescaled vibrational rates
        ## Case 2 ## r_M > r_m && T<=T_0_r  -> then the rovib rate is just the rescaled rotational rate
        ## Case 3 ## r_M > r_m && T_0_r < T < T_0_DM  -> then the rovib cooling is the linearly extrapolated rescaled rotational rate from T_0_r
        ## Case 4 ## r_M > r_m && T_0_DM < T < T_0_v  -> then the rovib cooling is the linearly extrapolated rescaled vibrational rate from T_0_v
        ## Case 5 ## r_M > r_m && T >= T_0_v  -> then the rovib cooling is just the rescaled vibrational rate





        ## Case 1 ##
        if self.r_M <= self.r_m:
            if T <= T_0_r and T < T_0_v:
                return (((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H_H2_DM_coll_r)

            elif T >= T_0_v and T > T_0_r:
                return (((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -0.75)) * H_H2_DM_coll_r)

            elif T <= T_0_r and T >= T_0_v:
                return (((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H_H2_DM_coll_r) + (((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -1.25)) * H_H2_DM_coll_v)

        ## Case 2 ##
        elif self.r_M > self.r_m and T <= T_0_r:
            return ((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H_H2_DM_coll_r

        ## Case 3 ##
        elif self.r_M > self.r_m and T > T_0_r and T <= T_0_DM:
            x_2 = T_0_r
            x_1 = 0.98 * x_2
            y_2 = self.H_H2_coll(x_2, 0,0)
            y_1 = self.H_H2_coll(x_1, 0,0)
            slope = (y_2 - y_1) / (x_2 - x_1)
            H_H2_lin = slope * (T_r - T_0_r) + y_2
            return ((self.r_a) * (self.r_m) * ((self.r_M) ** -2)) * H_H2_lin

        ## Case 4 ##
        elif self.r_M > self.r_m and T > T_0_DM and T <= T_0_v:
            x_2 = T_0_v
            x_1 = 1.01 * x_2
            y_2 = self.H_H2_coll(x_2, 0,0)
            y_1 = self.H_H2_coll(x_1, 0,0)
            slope = (y_2 - y_1) / (x_2 - x_1)
            H_H2_lin = slope * (T_v - T_0_v) + y_2
            return ((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -1.25)) * H_H2_lin

        ## Case 5 ##
        elif self.r_M > self.r_m and T > T_0_v:
            return ((self.r_a) * ((self.r_m) ** 0.25) * ((self.r_M) ** -1.25)) * H_H2_DM_coll_v



    def p_H2_DM_coll(self, T):
        TscaleR = (self.r_a ** -2) * (self.r_m ** -2) * (self.r_M ** 1)
        TscaleV = (self.r_a ** -2) * (self.r_m ** -1.5) * (self.r_M ** 0.5)
        TscaleA = (self.r_a ** -2) * (self.r_m ** -1) * (self.r_M ** 0)
        T_r = (((self.r_a) ** -2) * ((self.r_m) ** -2) * ((self.r_M) ** 1)) * T
        winL = 0
        winH = np.log10((TscaleR/TscaleA))
        
        return (((self.r_a) ** 1) * ((self.r_m) ** 0.5) * ((self.r_M) ** -1.5)) * self.p_H2_coll(T_r, winL, winH)


    def e_H2_DM_coll(self, T):
        TscaleR = (self.r_a ** -2) * (self.r_m ** -2) * (self.r_M ** 1)
        TscaleV = (self.r_a ** -2) * (self.r_m ** -1.5) * (self.r_M ** 0.5)
        TscaleA = (self.r_a ** -2) * (self.r_m ** -1) * (self.r_M ** 0)
        T_r = (((self.r_a) ** -2) * ((self.r_m) ** -2) * ((self.r_M) ** 1)) * T
        winL = 0
        winH = np.log10((TscaleR/TscaleA))
        return (((self.r_a) ** 1) * ((self.r_m) ** 0) * ((self.r_M) ** -1)) * self.e_H2_coll(T_r, winL, winH)

############################################
############################################


############################################
######## Total Molecular Cooling Rates
############################################
            ## With populations #

    ## Total DM LD Molecular Cooling ##
    def total_DM_LD_cool(self, n, T):
        return ((self.x_H2 * self.H2_H2_DM_coll(T))+ (self.x_H *self.H_H2_DM_coll(T)) + (self.x_e * self.e_H2_DM_coll(T)) + (self.x_p * self.p_H2_DM_coll(T)))


############################################
############################################





############################################
############################################





############################################
######## Local Thermal Equilibrium Molecular Cooling
######## From KROME, Grassi et al. 2013
############################################



############################################
######## SM LTE Rates
############################################

    def R_SM(self, T):
        T_3 = T/1000
        return (((9.5 * 10 ** -22 ) * (T_3 ** 3.76))/(1 + 0.12 * (T_3 ** 2.1))) * (np.exp(-(0.13/T_3) ** 3)) + (3 * 10 ** -24) * np.exp(-0.51/T_3)

    def V_SM(self, T):
        T_3 = T/1000
        return (6.7 * 10 ** -19) * np.exp(-5.86/T_3) + (1.6 * 10 ** -18) * np.exp(-11.7/T_3)

############################################
############################################



############################################
######## DM LTE Rates
############################################

    def R_DM(self, T):
        return (self.r_m ** 8) * (self.r_M ** -6) * (self.r_a ** 9) * self.R_SM(((self.r_a) ** (-2)) * ((self.r_m) ** -2) * ((self.r_M) ** (1)) * T)

    def V_DM(self, T):
        return (self.r_m ** 5) * (self.r_M ** -3) * (self.r_a ** 9) * self.V_SM(((self.r_a) ** (-2)) * ((self.r_m) ** (-3/2)) * ((self.r_M) ** (1/2)) * T)

############################################
############################################

    def HDL(self, T):
        #if T < 2000:
        return self.R_DM(T) + self.V_DM(T)

############################################
############################################

    def Chem_Cooling(self, n, T):
        return (((self.r_a ** 2) * (self.r_m) * (4.48 * E)) * self.n_H2 * self.n_H * self.k8) + (((self.r_a ** 2) * (self.r_m) * (1.83 * E)) * self.n_Hneg * self.n_H * self.k4)


############################################
######## Total DM Molecular Cooling Rates
######## Ryan et al 2022
############################################

    def total_DM_Molecular_Cooling(self, n, T): ## ergs/cm*s ##
        
        R_LTE_DM = self.R_DM(T)
        V_LTE_DM = self.V_DM(T)
        tot_LTE_DM = self.HDL(T)
        tot_LD_DM = self.total_DM_LD_cool(n, T)
        if tot_LD_DM <= 1e-215:
            return 0 * self.x_H2 * tot_LTE_DM
        else:
            return ((self.x_H2 * (R_LTE_DM + V_LTE_DM))/(1 + ((R_LTE_DM + V_LTE_DM)/(tot_LD_DM)))) + self.Chem_Cooling(n, T)

############################################
############################################




class heating:
    def __init__(self, n,T, x_H2_pre , r_m, r_M, r_a, xi, epsilon, z):
        ############################################
        ################ Parameters
        ############################################
        ADM = aDM(r_m, r_M, r_a, xi, epsilon, z)
        abund = abundances(n, T, x_H2_pre, r_m, r_M, r_a, xi, epsilon, z)
        RN = react_network(n, r_m, r_M, r_a, xi, epsilon, z)

        
        ## Dark Parameter Rescaling ##
        self.r_m = r_m
        self.r_M = r_M
        self.r_a = r_a

        ## Dark Parameter Values ##
        self.m = m_SM * ADM.r_m
        self.M = M_SM * ADM.r_M
        self.a = a_SM * ADM.r_a
        self.xi = xi
        self.epsilon = epsilon


        ## Parameters Specific to the Halo ##
        self.n = n
        self.T = T

        ## Population Fractions ##
        self.x_H2 = abund.x_H2
        self.x_Hneg = abund.x_Hneg
        self.x_H2pl = abund.x_H2pl
        self.x_H3pl = abund.x_H3pl
        self.x_e = abund.x_e
        self.x_p = abund.x_p
        self.x_H = abund.x_H


        ## Populations ##
        self.n_H2 = self.x_H2 * self.n
        self.n_H2pl = self.x_H2pl * self.n
        self.n_Hneg = self.x_Hneg * self.n
        self.n_H3pl = self.x_H3pl * self.n
        self.n_H = self.x_H * self.n
        self.n_p = self.x_p * self.n
        self.n_e = self.x_e * self.n

        ## Rates ##
        self.k3 = RN.k3(self.T)
        self.k7 = RN.k7(self.T)
        self.k11 = RN.k11(self.T)
        self.k12 = RN.k12(self.T)

    def mu(self, n):
        E = self.n_e * self.m
        P = self.n_p * self.M
        H = self.n_H * (self.m + self.M)
        Hneg = self.n_Hneg * (2 * self.m + self.M)
        H2 = self.n_H2 * (2 * self.m + 2 * self.M)
        H2pl = self.n_H2pl * (self.m + 2 * self.M)

        return (E + P + H + Hneg + H2 + H2pl)/(n * self.M)

    def t_ff(self, n):
        rho = self.mu(n) * self.M * n
        return ((3 * np.pi)/(32 * G_N * rho)) ** (1/2)

    def Compressional_Heating(self, n, T):
        return ((K_B * T)/ (n * self.t_ff(n)))

    def Chemical_Heating(self, n, T):
        return (((self.r_m) * (self.r_a ** 2) * 3.53 * E) * self.x_Hneg * self.x_H * self.k3) + (((self.r_m) * (self.r_a ** 2) * 1.83 * E) * self.x_H2 * self.x_p * self.k7) + (((self.r_m) * (self.r_a ** 2) *4.48 * E) * n * (self.x_H ** 3) * self.k11) + (((self.r_m) * (self.r_a ** 2) * 4.48 * E) * n * (self.x_H ** 2) * self.x_H2 * self.k12) 

    def Heating(self, n, T):
        return self.Chemical_Heating(n, T) + self.Compressional_Heating(n, T)


############################################
############################################

def get_Lambda(n,T, x_H2_pre , r_m, r_M, r_a, xi, epsilon, z):
    # Make Lambda return heating #
    #return heating(n,T, x_H2_pre , r_m, r_M, r_a, xi, epsilon, z).Heating(n, T)
    # Make Lambda return cooling #
    return molecular(n,T, x_H2_pre , r_m, r_M, r_a, xi, epsilon, z).total_DM_Molecular_Cooling(n, T) + atomic(n, T, x_H2_pre, r_m, r_M, r_a, xi, epsilon, z).total_dark_atomic_cooling(n, T)


"""



function get_Lambda(n, T; adm=ADM())
    x_H2_pre = 1e-6
    r_m = adm.rm
    r_M = adm.rM
    r_α = adm.rα
    ϵ = adm.ϵ
    ξ = adm.ξ
    z = adm.z
    py"get_Lambda($n, $T, $x_H2_pre ,$r_m , $r_M, $r_α, $ξ, $ϵ, $z)"
end

#T = 10 .^ range(3, stop=7, length=20)
#PP = abs.(get_Lambda.(1e0, T))
#loglog(T, PP)
#ylim(1e-35, 1e-20)
#xlim(7e3, 1e7)
#grid()


end
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
    integer,parameter::rstep = 100000
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
    ! Tgas = $
    ! ntot = $
    ! zred = $
    ! rm = $
    ! rMH = $
    ! ra = $
    ! xi = $
    ! epsilon = $
    ! t_lim = $        ! Must be in cgs !

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
    dt = t_lim/1000

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
	    print '(e0.4)',cooling(x(:),Tgas)
	    exit
        end if

    end do

end program test_equillibrium_darkkrome
!############### MODULE ##############
module starting_densities

  public abundance

  type abundance
    real*8,dimension(:,:,:,:),allocatable :: data
    real*8,dimension(:),allocatable :: parE
    real*8,dimension(:),allocatable :: parP
    real*8,dimension(:),allocatable :: parA
    real*8,dimension(:),allocatable :: parX
  end type
  
  type bounds
    integer*8::x1,x2
    real*8::x
  end type

contains

  !**************************
  function get_n(ntot,mE,mP,alphaD,xi,epsilon)
    use  krome_user
    implicit none
    real*8::get_n(krome_nmols),n(krome_nmols),ntot,mE,mP,alphaD,xi,epsilon

    n(:) = 1d-40/ntot
    n(KROME_idx_QH2) = interpolate(mE,mP,alphaD,xi,getSpeciesData("h2"))
    n(KROME_idx_QH2) = 10**(n(KROME_idx_QH2))
    n(KROME_idx_QH2j) = interpolate(mE,mP,alphaD,xi,getSpeciesData("h2p"))
    n(KROME_idx_QH2j) = 10**(n(KROME_idx_QH2j))
    n(KROME_idx_QH3j) = interpolate(mE,mP,alphaD,xi,getSpeciesData("h3p"))
    n(KROME_idx_QH3j) = 10**(n(KROME_idx_QH3j))
    n(KROME_idx_QHk) = interpolate(mE,mP,alphaD,xi,getSpeciesData("hm"))
    n(KROME_idx_QHk) = 10**(n(KROME_idx_QHk))
    n(KROME_idx_QHj) = interpolate(mE,mP,alphaD,xi,getSpeciesData("p"))
    n(KROME_idx_QHj) = 10**(n(KROME_idx_QHj))/epsilon
    n(KROME_idx_QE) = interpolate(mE,mP,alphaD,xi,getSpeciesData("e"))
    n(KROME_idx_QE) = 10**(n(KROME_idx_QE))/epsilon
    n(KROME_idx_QH) = 1-(2*n(KROME_idx_QH2)+2*n(KROME_idx_QH2j) &
    		 +3*n(KROME_idx_QH3j)+1*n(KROME_idx_QHk) &
    		 +1*n(KROME_idx_QHj)+0*n(KROME_idx_QE))

    n = ntot*n
    get_n = n

  end function get_n

  !**************************
  function get_nCE(ntot,T,mE,mP,alphaD,xi,epsilon)
    use  krome_user
    implicit none
    real*8::get_nCE(krome_nmols),n(krome_nmols),T,ntot,mE,mP,alphaD,xi,epsilon
    real*8::jxe,xe,x2,xp,xH
    
    ! Get primordial abundances from James data
    n(:) = get_n(ntot,mE,mP,alphaD,xi,epsilon)
    
    ! Compute James ionization rate
    jxe = n(KROME_idx_QE)/n(KROME_idx_QH)
    
    ! Compute ionization rate as n_e/n_H0 = coll_ion/recomb
    ! Note we're ignoring H2+ and H- reactions as subdominant
    xe = ionization_rate(T,mE,alphaD)
    
    ! Compute molecularization rate as n_H2/n_H2_primordial
    ! Assume H2 will only be destroyed, not created
    x2 = molecularization_rate(T,mE,alphaD)
    
    xp = xe/(1+xe)
    xH = 1/(xe+1)
    
    if (xe.LE.jxe) then
        ! Use primordial abundances
        ! scale H2 by mole rate
        n(KROME_idx_QH2) = n(KROME_idx_QH2)*x2
        get_nCE = n
        return
    else
        ! Predominantly ionised regime
        n(KROME_idx_QH) = ntot * xH
        n(KROME_idx_QE) = ntot * xp
        n(KROME_idx_QHj) = n(KROME_idx_QE)
        n(KROME_idx_QH2) = x2 * n(KROME_idx_QH2)
    end if
    
    
    ! We'll neglect H2+ and H- for now
    n(KROME_idx_QH2j) = 1d-40/ntot
    n(KROME_idx_QHk) = 1d-40/ntot
    
    
    get_nCE = n
  end function get_nCE
  
  !**************************
  function ionization_rate(T,mE,alphaD)
    use krome_constants
    implicit none
    real*8::ionization_rate,T,mE,alphaD
    real*8::y2
    real*8::f_approx,ei,coll_ion
    real*8::rlarge,rsmall,recomb
    
    y2 = mE * (alphaD**2.d0) / (2.d3 * boltzmann_eV * T) 
    
    ! Collisional Ionization rate
    f_approx = exp(-y2)/2 + (y2*ei(-y2))/2
    coll_ion = 2.2d-7 * (5.11d2 / mE)**(1.5d0) * sqrt(1.d5 / T) * f_approx
    
    ! Recombination rate
    rlarge = 8.4d-14 * (alphaD/1.d-2)**3.d0 * (511./mE)**(3/2) *&
             sqrt(1.d5 / T) * ( 1.744d0 + log(y2) + 1.d0 / (6.d0 * y2) )
    rsmall = 1.3d-15 * (alphaD/1.d-2)**5.d0 * sqrt(511./mE) *&
             (1.d6 / T)**(1.5d0) * ( -4.66 - 15d0 * log(y2) +&
              y2 * (5.42d0 - 1.4d1 * log(y2) ) )
    if (y2 > 0.25d0) then
      recomb = rlarge
    else
      recomb = rsmall
    end if
    
    ionization_rate = coll_ion/recomb
  
  end function ionization_rate

  !**************************
  function molecularization_rate(T,mE,alphaD)
    use krome_constants
    implicit none
    real*8::molecularization_rate,T,mE,alphaD
    real*8::deltaE
    
    deltaE = 4.48 * (mE/5.11d2) * (alphaD/0.00729)**2.d0
    
    molecularization_rate = 1.d0 / 2.d0 * exp(-(boltzmann_eV*T)/deltaE)
    
  end function molecularization_rate

  !**************************
  function getSpeciesData(specName)
    use  HDF5
    implicit none
    character(*) :: specName
    character(13) :: filename = "abundances.h5"
    character(16) :: abunddataset
    character(16) :: Edataset = "pars/mvals"
    character(16) :: Pdataset = "pars/Mvals"
    character(16) :: Adataset = "pars/avals"
    character(17) :: Xdataset = "pars/xivals"
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: space_id       ! Dataspace identifier
    integer(HID_T) :: dtype_id       ! Dataspace identifier
    integer(HSIZE_T),dimension(4) :: abund_data_dims,abund_max_dims
    integer(HSIZE_T),dimension(1) :: data_dims,max_dims
    integer :: error
    real*8,dimension(:),allocatable :: parE
    real*8,dimension(:),allocatable :: parP
    real*8,dimension(:),allocatable :: parA
    real*8,dimension(:),allocatable :: parX
    real*8,dimension(:,:,:,:),allocatable :: data
    type(abundance)::getSpeciesData

    abunddataset = "abundances/x" // trim(specName)

    call h5open_f(error)

    call h5fopen_f(filename,H5F_ACC_RDONLY_F, file_id, error)

    ! read abundance data
    call h5dopen_f(file_id,abunddataset,dset_id,error)
    call h5dget_space_f(dset_id,space_id,error)
    call h5sget_simple_extent_dims_f(space_id,abund_data_dims, abund_max_dims, error)

    ALLOCATE(data(abund_data_dims(1), abund_data_dims(2), abund_data_dims(3), abund_data_dims(4)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, abund_data_dims, error)

    getSpeciesData%data = log10(data+1.d-300)

    
    ! read mE parameter data
    call h5dopen_f(file_id,Edataset,dset_id,error)
    call h5dget_space_f(dset_id,space_id,error)
    call h5sget_simple_extent_dims_f(space_id,data_dims, max_dims, error)

    ALLOCATE(parE(data_dims(1)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, parE, data_dims, error)

    getSpeciesData%parE = parE

    
    ! read mP parameter data
    call h5dopen_f(file_id,Pdataset,dset_id,error)
    call h5dget_space_f(dset_id,space_id,error)
    call h5sget_simple_extent_dims_f(space_id,data_dims, max_dims, error)

    ALLOCATE(parP(data_dims(1)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, parP, data_dims, error)

    getSpeciesData%parP = parP

    
    ! read alpha parameter data
    call h5dopen_f(file_id,Adataset,dset_id,error)
    call h5dget_space_f(dset_id,space_id,error)
    call h5sget_simple_extent_dims_f(space_id,data_dims, max_dims, error)

    ALLOCATE(parA(data_dims(1)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, parA, data_dims, error)

    getSpeciesData%parA = parA

    
    ! read xi parameter data
    call h5dopen_f(file_id,Xdataset,dset_id,error)
    call h5dget_space_f(dset_id,space_id,error)
    call h5sget_simple_extent_dims_f(space_id,data_dims, max_dims, error)

    ALLOCATE(parX(data_dims(1)))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, parX, data_dims, error)

    getSpeciesData%parX = parX


    call h5close_f(error)

    end function getSpeciesData

  function interpolate(mE,mP,alphaD,xi,abun)
    integer:: i,ND1,ND2
    real*8::interpolate,mE,mP,alphaD,rE,rP,rA,xi
    real*8::v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16
    type(abundance) :: abun
    type(bounds):: boundsX,boundsY,boundsZ,boundsW
    real*8,dimension(:,:,:,:),allocatable :: data

    data = abun%data

    rE = mE/511
    rP = mP/0.938
    rA = alphaD/0.00729

    boundsZ = getBounds(rE,abun%parE)
    boundsW = getBounds(rP,abun%parP)
    boundsY = getBounds(rA,abun%parA)
    boundsX = getBounds(xi,abun%parX)

    v1 = data(boundsX%x1,boundsY%x1,boundsZ%x1,boundsW%x1) ! 0000
    v2 = data(boundsX%x1,boundsY%x1,boundsZ%x1,boundsW%x2) ! 0001
    v3 = data(boundsX%x1,boundsY%x1,boundsZ%x2,boundsW%x1) ! 0010
    v4 = data(boundsX%x1,boundsY%x1,boundsZ%x2,boundsW%x2) ! 0011
    v5 = data(boundsX%x1,boundsY%x2,boundsZ%x1,boundsW%x1) ! 0100
    v6 = data(boundsX%x1,boundsY%x2,boundsZ%x1,boundsW%x2) ! 0101
    v7 = data(boundsX%x1,boundsY%x2,boundsZ%x2,boundsW%x1) ! 0110
    v8 = data(boundsX%x1,boundsY%x2,boundsZ%x2,boundsW%x2) ! 0111
    v9 = data(boundsX%x2,boundsY%x1,boundsZ%x1,boundsW%x1) ! 1000
    v10 = data(boundsX%x2,boundsY%x1,boundsZ%x1,boundsW%x2) ! 1001
    v11 = data(boundsX%x2,boundsY%x1,boundsZ%x2,boundsW%x1) ! 1010
    v12 = data(boundsX%x2,boundsY%x1,boundsZ%x2,boundsW%x2) ! 1011
    v13 = data(boundsX%x2,boundsY%x2,boundsZ%x1,boundsW%x1) ! 1100
    v14 = data(boundsX%x2,boundsY%x2,boundsZ%x1,boundsW%x2) ! 1101
    v15 = data(boundsX%x2,boundsY%x2,boundsZ%x2,boundsW%x1) ! 1110
    v16 = data(boundsX%x2,boundsY%x2,boundsZ%x2,boundsW%x2) ! 1111

    interpolate = interpolate4D(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,boundsW%x,boundsZ%x,boundsY%x,boundsX%x)

  end function interpolate

  function getBounds(val,par)
    integer*8::x1,x2,ND
    real*8 :: x,val
    real*8,dimension(:) :: par
    type(bounds) :: getBounds

    ND = SIZE(par,1)

    x1=1
    x2=1
    do i=2,ND
        if (val >= par(i)) then
            x1=i
        end if
        if (val > par(i-1)) then
            x2=i
        end if
    end do
    if (x2.NE.x1) then
        x = (val-par(x1))/(par(x2)-par(x1))
    else if(val<par(x1)) then
        x2 = x1+1
        x = (val-par(x1))/(par(x2)-par(x1))
    else
        x = 0
    end if

    getBounds%x1 = x1
    getBounds%x2 = x2
    getBounds%x = x
  end function getBounds

  function interpolate1D(v1,v2,x)
    real*8::interpolate1D,v1,v2,x
    interpolate1D = v1*(1-x) + x*v2
  end function interpolate1D

  function interpolate2D(v1,v2,v3,v4,x,y)
    real*8::x,y,s,t
    real*8::interpolate2D,v1,v2,v3,v4
    s = interpolate1D(v1,v2,x)
    t = interpolate1D(v3,v4,x)
    interpolate2D = interpolate1D(s,t,y)
  end function interpolate2D

  function interpolate3D(v1,v2,v3,v4,v5,v6,v7,v8,x,y,z)
    real*8::x,y,z,s,t
    real*8::interpolate3D,v1,v2,v3,v4,v5,v6,v7,v8
    s = interpolate2D(v1,v2,v3,v4,x,y)
    t = interpolate2D(v5,v6,v7,v8,x,y)
    interpolate3D = interpolate1D(s,t,z)
  end function interpolate3D

  function interpolate4D(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,x,y,z,w)
    real*8::x,y,z,w,s,t
    real*8::interpolate4D,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16
    s = interpolate3D(v1,v2,v3,v4,v5,v6,v7,v8,x,y,z)
    t = interpolate3D(v9,v10,v11,v12,v13,v14,v15,v16,x,y,z)
    interpolate4D = interpolate1D(s,t,w)
  end function interpolate4D

end module starting_densities

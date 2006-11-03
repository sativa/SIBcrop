!...module to hold values for constants that are not contained in the
!...BUGS module physical_parameters.F

module sib_const_module

    use kinds
    use physical_parameters

    implicit none
    save

    real(kind=dbl_kind) :: dtt       ! model time step (seconds) 
    real(kind=dbl_kind) :: dti       ! inverse time step

    !--variables that retain a constant value throughout the simulation
    integer(kind=int_kind) ::     &
        nsib,           & !  number of SiB points in datasets
        subcount,       & !  actual number of SiB point in simulation
        snowl,          & !  number of (actual) snow layers
        ihr,            & !  global points in x-direction
        jhr,            & !  global points in y-direction
        nhr               !  ihr*jhr (total global points)
        
    integer(kind=long_kind) :: &
	endtime,        & !  end time of integration -- units can vary
        starttime,      & !  start time of integration -- units can vary
	dtsib,          & !  timestep in seconds
        dtsibmetin,     & !  driver data input interval (seconds)
        dtsibres,       & !  restart interval (see namel_sibdrv for exp)
        dtsibout	  !  output interval (see namel_sibdrv for exp)

     integer(kind=int_kind) :: &
        dtsibbcin,      & !  sib boundary condition input interval
        nsecond,        & !  simulation time, in seconds
        numsib,         & !  check of # of sib points--init_sibdrv
        nstepsib          !  number of timesteps integrated

    integer(kind=int_kind) :: endyear
    integer(kind=int_kind) :: ndtsibpbp

    !itb...to nsib points in init_sibdrv.F
    real(kind=real_kind), dimension(:), allocatable    ::    &
        latsib,       & !  SiB point latitude
        lonsib,       & !  SiB point longitude
        latitude,     &
        longitude,    &
        lonpbp,       &
        latpbp

    integer (kind=int_kind), parameter ::          &
        nsoil = 10,                 &!  number of soil layers         
        nsnow = 5,                  &!  maximum number of snow layers       
        physmax = 5                  !  maximum number of physiology types
                                     !    (only C3 and C4 for now, but 
                                     !     the capability for more)

    !itb...some extra variables
    real(kind=real_kind) ::  &
        sin_dec,        & !  function of earth's orbit
        cos_dec,        & !  function of earth's orbit
        tau

    integer(kind=int_kind) ::     & ! 
        startyear = 1,   & !  time manager stuff
        eqnx   = 80     !  day of vernal equinox

    real(kind=real_kind) ::   &
        lonearth        ! (rad) Earth lon about Sun from vernal equinox

    integer (kind=int_kind), dimension(:), allocatable   ::    &
        subset,       & !  array of landpoint indices for subgrid
        latindex,     & !  latitude index array of all landpoints
        lonindex,     & !  longitude index array of all landpoints
        sublat,       & !  latitude index array of subset
        sublon          !  longitude index array of subset

    !itb...some SCALARS
    real(kind=dbl_kind) :: &
        c3day,        & !  timesteps per day
        ztemp,        & !  height of temperature measurement (m)
        zwind           !  height of wind measurement (m)

    !------------------------------------------------------------------------
    real (kind=dbl_kind), parameter  ::     &
        version = 3.0,                &!  code version identifier
        snomel  = 3.705185e8,         &!  latent heat of fusion of ice (J m^-3) 
        cv     = 1952.0,              &!  specific heat of water vapor at 
                                       !  constant pressure (J deg^-1 kg^-1)
        cpice  = 2117.27,             &!  specific heat of ice (J kg^-1 deg^-1)
        cpliq  = 4188.0,              &!  spec heat of water (J kg^-1 deg^-1)

!itb...playing with vegetation heat capacity

!itb...original.....
        clai   = 4.186*1000.0*0.2,    &!  leaf heat capacity  (J m^-2 deg^-1)
        cww    = 4.186*1000.0*1000.0, &!  water heat capacity (J m^-3 deg^-1)
        asnow  = 16.7,                &!  UNKNOWN
        rotper = 24.0,                &!  hours per day
        day    = rotper * 3600.0,     &!  seconds per day
        vkrmn  = 0.35,                &!  Von Karmann's constant (unitless)
        ribc   = 3.05,                &!  critical Richardson Number (unitless)
        pr0    = 0.74,                &!  turb Prandtl Number at neutral stblty
        tkemin = 0.01,                &!  minimum allowed value for tke
        rgfac  = 100.0/gas_const_r,   &!  
        cpdgrv = spec_heat_cp/grav,   &! 
        po2m   = 20900.0,             &!  mixed layer O2 concentration
        perhl  = 102.7,               &!  UNKNOWN     
        denh2o = 1000.0,              &!  density of water (kg/m^3)
        denice = 917.0,               &!  density of ice (kg/m^3) 
        tkair  = 0.023,               &!  thermal conductivity of air (W/m/K)
        tkwat  = 0.6,                 &!  thermal conductivity of water (W/m/K)
        tkice  = 2.29,                &!  thermal conductivity of ice (W/m/K)
        snofac = hltm/(hltm + snomel * 1.E-3), & !  ratio of hltm to hltm+ht 
                                       ! of fusion (see Sellers (1986) appendix B)
        wimp = 0.05,                  &!  water impermeable if porosity 
                                       !  below this value
        phmin = -1.e8,                &!  minimum value for soil potential (mm)
        !vwcmin = 0.1,                 &!  wilting point volumetric water content
        wtfact = 0.3,                 &!  fraction of area with high 
                                       !  (HARDWIRE PATCH) water table
        ssi = 0.033,                  &!  irreducible water fraction of snow
        zlnd = 0.01,                  &!  roughness length for land (m)
        eccn   = 0.016715,            &!  eccentricity
        daypyr = 365.0,               &!  days per  year
        decmax = 23.441,              &!  max declination
        cn_fact = 0.5,                &!  Crank-Nicholson factor
        cosz_min = -0.1045            !  minimum cosine of zenith angle value
                                       !  -0.1045 is 96 deg. which includes
                                       !  civil twilight





    !...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX
    real(kind=dbl_kind),parameter ::        &
        pdb = 0.0112372   
    ! 13C/12C ratio of Pee Dee 
    !   Belemnite (no units)

    !...Carbon isotopic fractionation constants (units = per mil)         
    !...KIEC refers to Kinetic Isotope Effect (KIE) for Carbon (C), 
    !...and can be converted to alpha notation by alpha = (1 - KIEC/1000). 
    !...For a chemical reaction, alpha = Rreactant/Rproduct.  
    !...KIEs are sometimes referred to as epsilon factors. 


    real(kind=dbl_kind),parameter ::  &
        kieclfbl  = - 2.9    
    ! canopy air space to leaf 
    !  boundary layer
    real(kind=dbl_kind),parameter ::  &
        kiecstom  = - 4.4    
    ! leaf boundary layer to 
    !  stomatal cavity
    real(kind=dbl_kind),parameter ::  &
        kieclphas =  -0.7    
    ! liquid phase fractionation 
    real(kind=dbl_kind),parameter ::  &
        kiecdis   =  -1.1    
    ! dissolution

    real(kind=dbl_kind),parameter ::  &
        kiecrbsco = -28.2    
    ! C3 C-fixation enzyme rubisco

    real(kind=dbl_kind),parameter ::  &
        tref = 298.16       
    ! standard temperature (K)

    real(kind=dbl_kind),parameter ::  &
        pref = 101325.0     
    ! standard pressure (Pa)


end module sib_const_module

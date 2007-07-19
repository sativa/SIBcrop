module sibtype

!----------------------------------------------------------------------
!
!   New SiB module: SiB returns to a single-point model, for 
!   reasons of compatability with BUGS and MPI offline code (SiBDRV)
!   as well as to adhere to the CCM4 coding standard
!
! Modifications
!  created Ian Baker, 06 Feb 2002
!  Kevin Schaefer changed tdew1/2 to sh1/2 (8/17/04)
!
!----------------------------------------------------------------------

use kinds
use sib_const_module, only: nsoil, nsnow

implicit none

!jlc...These need to be public
public  sib_t
public  sib_local_vars

!jlc...Not sure what these need to be
public param_vars
public prognostic_vars
public diagnostic_vars
public sib_status


!------------------------------------------------------------------
!                   BOUNDARY CONDITION VARIABLES
!------------------------------------------------------------------
type param_vars


    !...boundary conditions--TIME INVARIANT
    real(kind=dbl_kind) :: biome     ! biome type (see refs for description)
    real(kind=dbl_kind) :: chil      ! leaf angle distribution factor (-)
    real(kind=dbl_kind) :: phc       ! 1/2 crit leaf water pot limit (m)
    real(kind=dbl_kind) :: z1        ! canopy bottom (m)
    real(kind=dbl_kind) :: z2        ! canopy top (m)
    real(kind=dbl_kind) :: poros     ! soil porosity (zero to one)
    real(kind=dbl_kind) :: satco     ! hydraulic conductivity at 
    !               saturation (m/s)
    real(kind=dbl_kind) :: bee       ! Clapp & Hornberber 'b' exponent (-)
    real(kind=dbl_kind) :: phsat     ! soil tension at saturation (m)
    real(kind=dbl_kind) :: slope     ! cosine of mean slope (-)
    real(kind=dbl_kind) :: vcover    ! fraction of vegetation cover (0-1)
    real(kind=dbl_kind) :: vmax0(5)  ! Rubisco vel of sun leaf (mol/m^2/sec)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: trop(5)   ! temp coeff in GS-A model (K)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: trda(5)   ! temp coeff in GS-A model (K^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: trdm(5)   ! temp coeff in GS-A model (K)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: respcp(5) ! respiration fraction of vmax0 (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: slti(5)   ! slope of lo-temp inhibition (K^-1)
    !   function (K^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: shti(5)   ! slop of hi-temp inhibition (K^-1)
    !   function (K-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: hltii(5)  ! 1/2 point of lo-temp inhibition (K)
    !   function
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: hhti(5)   ! 1/2 point of hi-temp inhibition (K)
    !   function
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: soref(2)  ! soil reflectance (-)
    !   (1) shortwave
    !   (2) longwave
    real(kind=dbl_kind) :: effcon(5) ! quantum efficiency (mol/mol)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: binter(5) ! conductance-photosynthesis 
    !   intercept (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: gradm(5)  ! conductance-photosynthesis slope 
    !   parameter (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: atheta(5) ! WC WE coupling parameter (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: btheta(5) ! WC&WE, WS coupling parameter (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: zm        ! respiration parameter - exponent
    real(kind=dbl_kind) :: wopt      ! respiration parameter - optimum (-)
    !  soil moisture
    real(kind=dbl_kind) :: wsat      ! respiration parameter (-)
    real(kind=dbl_kind) :: sandfrac  ! soil sand fraction
    real(kind=dbl_kind) :: clayfrac  ! soil clay fraction
    real(kind=dbl_kind) :: physfrac(5)
    ! physiology fraction-5 elements must add
    !  up to 1.0
    !  (1) C3 vegetation
    !  (2) C4 vegetation
    !  (3) open
    !  (4) open
    !  (5) open
    real(kind=dbl_kind) :: physfrac1(5)
    ! physiology fraction-5 elements must add
    !  up to 1.0
    !  (1) C3 vegetation
    !  (2) C4 vegetation
    !  (3) open
    !  (4) open
    !  (5) open
    real(kind=dbl_kind) :: physfrac2(5)
    ! physiology fraction-5 elements must add
    !  up to 1.0
    !  (1) C3 vegetation
    !  (2) C4 vegetation
    !  (3) open
    !  (4) open
    !  (5) open

    integer(kind=int_kind) :: phystype(5)
    ! physiology type-will be either 3 or 4
    !  (1) C3 vegetation - always 3
    !  (2) C4 vegetation - always 4
    !  (3) open
    !  (4) open
    !  (5) open
    real(kind=dbl_kind) :: rootf(nsoil)
    ! root fraction
    real(kind=dbl_kind) :: rootr(nsoil)
    ! adjusted root fraction
    real(kind=dbl_kind) :: tran(2,2) ! leaf transmittance (-)
    !  (1,1) - shortwave, green plants
    !  (1,2) - longwave, green plants
    !  (2,1) - shortwave, brown plants
    !  (2,2) - longwave, brown plants
    real(kind=dbl_kind) :: ref(2,2)  ! leaf reflectance (-)
    !  (1,1) - shortwave, green plants
    !  (1,2) - longwave, green plants
    !  (2,1) - shortwave, brown plants
    !  (2,2) - longwave, brown plants
    real(kind=dbl_kind) :: respfactor(nsoil) 
    ! factor for balancing soil respiration
    !  with annual assimilation (-)

    real(kind=dbl_kind) :: tkmg(nsoil)   ! thermal conductivity, soil 
    ! minerals (W m^-1 K^-1)
    real(kind=dbl_kind) :: tksatu(nsoil) ! thermal conductivity, saturated
    !   soil (W m^-1 K^-1)
    real(kind=dbl_kind) :: tkdry(nsoil)
    ! thermal conductivity, dry soil 
    !                   (W m^-1 K^-1)      
    real(kind=dbl_kind) :: tksoil(-nsnow+1:nsoil) 
    ! ground and snow thermal 
    !  conductivity (W m^-1 K^-1)
    real(kind=dbl_kind) :: slamda(-nsnow+1:nsoil) 
    ! CLM heat flux term
    !   (see begtem.F) (m^2 K W^-1)
    real(kind=dbl_kind) :: csolid(nsoil)
    ! heat capacity, soil solids
    !                         (J/m^3/K)
    real(kind=dbl_kind) :: shcap(-nsnow+1:nsoil)
    ! soil total heat capacity
    !                         (J/m^2/K)
    real(kind=dbl_kind) :: satcap(2) ! saturation capacity depth for vegetation
    !   (1) and ground (2) (m)
    real(kind=dbl_kind) :: czc       ! canopy heat capacity (J m^-2 K-1)

    real(kind=dbl_kind) :: vwcmin    ! soil wilting point (volumetric)
    real(kind=dbl_kind) :: fieldcap  ! soil field capacity (volumetric)




    !...boundary conditions-TIME VARYING
    real(kind=dbl_kind) :: aparc      ! absorbed fraction of PAR (-)
    real(kind=dbl_kind) :: aparc1     ! absorbed fraction of PAR (-)
    real(kind=dbl_kind) :: aparc2     ! absorbed fraction of PAR (-)
    real(kind=dbl_kind) :: zlt        ! leaf area index (-)
    real(kind=dbl_kind) :: zlt1       ! leaf area index (-)
    real(kind=dbl_kind) :: zlt2       ! leaf area index (-)
    real(kind=dbl_kind) :: green      ! green fraction of LAI (-)
    real(kind=dbl_kind) :: green1     ! green fraction of LAI (-)
    real(kind=dbl_kind) :: green2     ! green fraction of LAI (-)
    real(kind=dbl_kind) :: z0d        ! roughness length (m)
    real(kind=dbl_kind) :: z0d1       ! roughness length (m)
    real(kind=dbl_kind) :: z0d2       ! roughness length (m)
    real(kind=dbl_kind) :: z0         ! roughness length adjust for 
    !   snow-covered canopy (m)
    real(kind=dbl_kind) :: z01        ! roughness length adjust for 
    !   snow-covered canopy (m)
    real(kind=dbl_kind) :: z02        ! roughness length adjust for 
    !   snow-covered canopy (m)
    real(kind=dbl_kind) :: zp_disp    ! zero-plane displacement (m)
    real(kind=dbl_kind) :: zp_disp1   ! zero-plane displacement (m)
    real(kind=dbl_kind) :: zp_disp2   ! zero-plane displacement (m)
    real(kind=dbl_kind) :: zpd_adj    ! zp_disp adjusted for snow on canopy (m)
    real(kind=dbl_kind) :: zpd_adj1   ! zp_disp adjusted for snow on canopy (m)
    real(kind=dbl_kind) :: zpd_adj2   ! zp_disp adjusted for snow on canopy (m)
    real(kind=dbl_kind) :: cc1        ! bulk pbl resistance coefficient (s/m)^0.5
    real(kind=dbl_kind) :: cc11       ! bulk pbl resistance coefficient (s/m)^0.5
    real(kind=dbl_kind) :: cc12       ! bulk pbl resistance coefficient (s/m)^0.5
    real(kind=dbl_kind) :: cc2        ! ground to CAS resistance (-)
    ! coefficient
    real(kind=dbl_kind) :: cc21       ! ground to CAS resistance (-)
    ! coefficient
    real(kind=dbl_kind) :: cc22       ! ground to CAS resistance (-)
    ! coefficient
    real(kind=dbl_kind) :: rbc        ! cc1 adjusted for snow (s/m)^0.5
    real(kind=dbl_kind) :: rbc1       ! cc1 adjusted for snow (s/m)^0.5
    real(kind=dbl_kind) :: rbc2       ! cc1 adjusted for snow (s/m)^0.5
    real(kind=dbl_kind) :: rdc        ! cc2 adjusted for snow (-)
    real(kind=dbl_kind) :: rdc1       ! cc2 adjusted for snow (-)
    real(kind=dbl_kind) :: rdc2       ! cc2 adjusted for snow (-)
    real(kind=dbl_kind) :: gmudmu     ! time-mean leaf projection (-)
    real(kind=dbl_kind) :: gmudmu1    ! time-mean leaf projection (-)
    real(kind=dbl_kind) :: gmudmu2    ! time-mean leaf projection (-)

    ! isotope stuff 
    real(kind=dbl_kind) :: d13cresp   ! del13C of respiration (per mil vs PDB)
    real(kind=dbl_kind) :: d13cresp1  ! del13C of respiration (per mil vs PDB)
    real(kind=dbl_kind) :: d13cresp2  ! del13C of respiration (per mil vs PDB)


end type param_vars



!------------------------------------------------------------------
!                   PROGNOSTIC VARIABLES
!------------------------------------------------------------------
type prognostic_vars


    real(kind=dbl_kind) :: td(-nsnow+1:nsoil)
    ! soil temperature (K)
    real(kind=dbl_kind) :: tg        ! ground surface temp (K)
    real(kind=dbl_kind) :: ta        ! CAS temperature (K)
    real(kind=dbl_kind) :: tc        ! vegetation temperature (K)
    real(kind=dbl_kind) :: tha       ! CAS potential temperature (K)
    real(kind=dbl_kind) :: sha       ! CAS water vapor mixing ratio (kg/kg)
    real(kind=dbl_kind) :: ea        ! CAS water vapor pressure (hPa or mb)

    real(kind=dbl_kind) :: tke       ! turbulent kinetic energy (UNITS??)

    !...soil/snow arrays

    real(kind=dbl_kind) :: www_liq(-nsnow+1:nsoil)
    ! soil liquid (kg/m^2)
    real(kind=dbl_kind) :: www_ice(-nsnow+1:nsoil)
    ! soil ice (kg/m^2)
    real(kind=dbl_kind) :: vol_liq(-nsnow+1:nsoil)
    ! soil liquid - unitless, fraction
    !   of total layer volume
    real(kind=dbl_kind) :: vol_ice(-nsnow+1:nsoil)
    ! soil liquid - unitless, fraction
    !   of total layer volume

    real(kind=dbl_kind) :: node_z(-nsnow+1:nsoil)
    ! layer node depth (m)
    real(kind=dbl_kind) :: layer_z(-nsnow:nsoil)
    ! layer interface depth (m)
    real(kind=dbl_kind) :: dz(-nsnow+1:nsoil)
    ! layer thickness (m)

    real(kind=dbl_kind) :: snow_veg  ! vegetation snow cover (kg/m^2)
    real(kind=dbl_kind) :: snow_mass ! mass of snow on ground (kg/m^2)
    real(kind=dbl_kind) :: snow_depth 
    ! depth of snow on ground (m)



    real(kind=dbl_kind) :: snow_age  ! non-dimensional snow age

    integer(kind=int_kind) :: nsl    ! number of (actual) snow layers
    !   THIS NUMBER WILL BE NEGATIVE WHEN
    !   SNOW IS PRESENT
    real(kind=dbl_kind) :: capac(2)  ! vegetation and ground surface liquid
    !   water interception storage (kg m^-2)
    !   (1) canopy   (2) ground

    real(kind=dbl_kind) :: rst(6)    ! stomatal resistance (sec/m)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    real(kind=dbl_kind) :: pco2ap     ! (Pa) CAS CO2 partial pressure
    real(kind=dbl_kind) :: pco2ap_old ! (Pa) previous timestep pco2ap
    real(kind=dbl_kind) :: cas        ! (mole m^-2) CO2 store in Canopy air space
    real(kind=dbl_kind) :: cas_old    ! (mole m^-2) previous timestep CO2 store in Canopy air space
    real(kind=dbl_kind) :: expand     ! (mole m^-2) CAS expansion loss from previous timestep CO2 store
    real(kind=dbl_kind) :: pco2m      ! mixed layer CO2 partial pressure (Pa)

    !...driver data/forcing
    real(kind=dbl_kind) :: sw_dwn    ! surface incident shortwave radiation (W/m^2)
    real(kind=dbl_kind) :: sw_dwn1   ! surface incident shortwave radiation (W/m^2)
    real(kind=dbl_kind) :: sw_dwn2   ! surface incident shortwave radiation (W/m^2)
    real(kind=dbl_kind) :: radvbc    ! visible beam radiation (W/m^2)
    real(kind=dbl_kind) :: radvdc    ! visible diffuse radiation (W/m^2)
    real(kind=dbl_kind) :: radnbc    ! nir beam radiation (W/m^2)
    real(kind=dbl_kind) :: radndc    ! nir diffuse radiation (W/m^2)
    real(kind=dbl_kind) :: dlwbot    ! surface incident longwave 
    real(kind=dbl_kind) :: dlwbot1   ! surface incident longwave 
    real(kind=dbl_kind) :: dlwbot2   ! surface incident longwave 
    real(kind=dbl_kind) :: vdcsav    ! jk used to save vdcsib2 for single point runs
    !                 radiation (W/m^2)
    real(kind=dbl_kind) :: tm        ! mixed layer temperature (K)
    real(kind=dbl_kind) :: tm1       ! mixed layer temperature (K)
    real(kind=dbl_kind) :: tm2       ! mixed layer temperature (K)
    real(kind=dbl_kind) :: thm       ! mixed layer potential temperature (K)
    real(kind=dbl_kind) :: sh        ! mixed layer water vapor mixing ratio (kg/kg)
    real(kind=dbl_kind) :: sh1       ! mixed layer water vapor mixing ratio (kg/kg)
    real(kind=dbl_kind) :: sh2       ! mixed layer water vapor mixing ratio (kg/kg)
    real(kind=dbl_kind) :: em        ! mixed layer water vapor 
    !  pressure (hPa or mb)
    real(kind=dbl_kind) :: ps        ! surface pressure (hPa or mb)
    real(kind=dbl_kind) :: ps1       ! surface pressure (hPa or mb)
    real(kind=dbl_kind) :: ps2       ! surface pressure (hPa or mb)
    real(kind=dbl_kind) :: bps(2)    ! (ps/1000)**kapa
    !   multiplying by bps turns a theta
    !   into a temperature
    real(kind=dbl_kind) :: psb       ! boundary layer mass depth (hPa or mb)
    real(kind=dbl_kind) :: zb        ! boundary layer thickness (m)
    real(kind=dbl_kind) :: ros       ! surface air density (kg/m^3)
    real(kind=dbl_kind) :: cupr      ! cumulus precipitation rate (mm/sec)
    real(kind=dbl_kind) :: cupr1     ! cumulus precipitation rate (mm/sec)
    real(kind=dbl_kind) :: cupr2     ! cumulus precipitation rate (mm/sec)
    real(kind=dbl_kind) :: lspr      ! stratiform precipitation rate (mm/sec) 
    real(kind=dbl_kind) :: lspr1     ! stratiform precipitation rate (mm/sec)  
    real(kind=dbl_kind) :: lspr2     ! stratiform precipitation rate (mm/sec) 

    real(kind=dbl_kind) :: spdm      ! wind speed (m/sec) 
    real(kind=dbl_kind) :: spdm1     ! wind speed (m/sec) 
    real(kind=dbl_kind) :: spdm2     ! wind speed (m/sec) 

    real(kind=dbl_kind) :: tcc1      ! cloud cover fraction
    real(kind=dbl_kind) :: tcc2      ! cloud cover fraction

    real(kind=dbl_kind) :: d13cca    ! del13C of canopy CO2 (per mil vs PDB)
    real(kind=dbl_kind) :: d13cm     ! del13C of ref level (per mil vs PDB)

end type prognostic_vars





!------------------------------------------------------------------
!                   DIAGNOSTIC VARIABLES
!------------------------------------------------------------------
type diagnostic_vars


!itb_crop...diagnostics for crop model development
    real(kind=dbl_kind) :: ta_bar    ! daily mean CAS air temp (K)
                                     ! used for calculating growing
                                     ! degree days
    real(kind=dbl_kind) :: gdd       ! growing degree days

    real(kind=dbl_kind) :: tb_temp(20000) 

                                     ! placeholder for accumulating
                                     ! temperature for GDD
    integer(kind=int_kind) :: tb_indx
                                     ! index for counting up timesteps

    integer(kind=int_kind) :: year

	integer(kind=int_kind) :: doy	!to calculate planting dates and growth stages of certain crops- EL
                                     

!itb...end crop variables


    real(kind=dbl_kind) :: eastar    ! CAS saturation vapor pressure (hPa or mb)
    real(kind=dbl_kind) :: rha       ! CAS relative humidity (-)
    real(kind=dbl_kind) :: psy       ! psycrometric constant (gamma) (hPa K^-1)
    real(kind=dbl_kind) :: salb(2,2) ! total albedo
    !   (1,1) - visible, beam
    !   (1,2) - visible, diffuse
    !   (2,1) - nir, beam
    !   (2,2) - nir, diffuse
    real(kind=dbl_kind) :: cas_cap_heat
    ! CAS heat capacity (J/m^2/K)
    real(kind=dbl_kind) :: cas_cap_vap
    ! CAS vapor capacity (J Pa^-1 m^-2)
    real(kind=dbl_kind) :: cas_cap_co2
    ! depth of 'canopy'  (m)

    real(kind=dbl_kind) :: cas_e_storage   ! CAS change in energy/timestep (W/m^-2)



    real(kind=dbl_kind) :: canex     ! snow depth on vegetation factor
    !     (unitless)
    real(kind=dbl_kind) :: wc        ! canopy wetness fraction (-)
    real(kind=dbl_kind) :: wg        ! ground wetness fraction (-)
    real(kind=dbl_kind) :: rstfac(4) ! stress factors (-)
    !  (1) leaf surface RH stress
    !  (2) rootzone water stress
    !  (3) temperature stress
    !  (4) product of factors 1-3
    real(kind=dbl_kind) :: wssp      ! water stress shape parameter. unitless.
                                     ! this variable controls the shape of the 
                                     ! water stress curve. value set/adjusted
                                     ! in begtem. 


    real(kind=dbl_kind) ::  paw_tot  ! paw, summed over all soil layers (kg m^-2)
    real(kind=dbl_kind) ::  paw_max  ! maximum paw, column (kg m^-2)
    real(kind=dbl_kind) ::  pawfrac  ! fraction of water available to plants
    real(kind=dbl_kind) ::  paw(nsoil) ! plant available water (volumetric)


    real(kind=dbl_kind) :: areas     ! snow cover fraction (zero to one)
    real(kind=dbl_kind) :: a_areas   ! 'apparent' areas, for computation
    !   purposes (will have value 0.0 or 1.0)
    real(kind=dbl_kind) :: tsnow     ! snow surface temp (K)
    real(kind=dbl_kind) :: snow_end(3) !julian day that snow ends
                                       ! array positions are for 3 criteria
                                       ! 1) snow depth =0
                                       ! 2) areas < threshold value (0.05)
                                       ! 3) snow_mass < threshold value (1.0)

    real(kind=dbl_kind) :: eff_poros(-nsnow+1:nsoil)
    ! effective porosity for liquid 
    !  (unitless)
    real(kind=dbl_kind) :: snowmelt  ! snow water converted to liquid
    !  and added to soil sfc (kg m^-2 sec^-1) 
    real(kind=dbl_kind) :: www_tot_soil  ! total soil water-all layers, water+ice (kg m^-2)
    real(kind=dbl_kind) :: roff      ! total subsurface runoff out of 
              !soil layers during whole sib timestep dtt (mm)
    real(kind=dbl_kind) :: roffo     ! overland runoff (mm)
    real(kind=dbl_kind) :: qqq       ! part of the subsurface runoff (mm)
    real(kind=dbl_kind) :: hr        ! soil surface relative humidity
    real(kind=dbl_kind) :: hrr       ! (copy) soil surface relative humidity
    real(kind=dbl_kind) :: soilscale(nsoil) 
    ! 'R-star' from Denning et al (1996) 
    !    (eqn 6) (UNITS?)
    real(kind=dbl_kind) :: tot_ss(13,nsoil) ! accumulated soilscale used for
                                     ! calculating respfactor
    real(kind=dbl_kind) :: soilq10(nsoil)
    ! soil temperature respiration dependence 
    !    function (UNITS?)
    real(kind=dbl_kind) :: respg     ! ground respiration (mol m^-2 sec^-1)
    real(kind=dbl_kind) :: www_inflow
    ! water inflow at ground surface 
    !    (kg m^-2 sec^-1)

    real(kind=dbl_kind) :: cu        ! momentum transfer coefficient (-)
    real(kind=dbl_kind) :: ct        ! thermal transfer coefficient (-)
    real(kind=dbl_kind) :: ustar     ! friction velocity (m sec^-1)
    real(kind=dbl_kind) :: drag(2)   ! drag (kg m^-2 sec^-1)
    real(kind=dbl_kind) :: ventmf    ! ventilation mass flux (kg m^-2 sec^-1)
    real(kind=dbl_kind) :: thvgm     ! sfc-reference height deficit of moisture
    !   UNITS ARE UNCLEAR

    real(kind=dbl_kind) :: ecmass    ! canopy evapotranspiration (kg m^-2 or 
    !                              mm water)
    real(kind=dbl_kind) :: egmass    ! ground evapotranspiration (kg m^-2 or 
    !                              mm water)
    real(kind=dbl_kind) :: chf       ! canopy heat storage flux (W m^-2)
    real(kind=dbl_kind) :: shf       ! soil heat storage flux (W m^-2)

    !...resistance
    real(kind=dbl_kind) :: ra        ! CAS-mixed layer resistance (sec/m) 
    real(kind=dbl_kind) :: rb        ! leaf-CAS resistance (sec/m)
    real(kind=dbl_kind) :: rc        ! canopy-CAS resistance (stomatal
    !   resistance + 2rb) (sec/m)
    real(kind=dbl_kind) :: rd        ! ground-CAS resistance (sec/m)
    real(kind=dbl_kind) :: rsoil     ! soil surface resistance (sec m^-1)
    real(kind=dbl_kind) :: rds       ! rd + rsoil (sec m^-1)

    real(kind=dbl_kind) :: ggl(6)    ! leaf conductance (sec/m)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    !...carbon

    real(kind=dbl_kind) :: pco2i(6)  ! leaf internal CO2 partial 
    !        pressure (Pa) 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes   
    real(kind=dbl_kind) :: pco2c(6)  ! chloroplast CO2 partial pressure (Pa)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: pco2s(6)  ! leaf surface CO2 partial
    !        pressure (Pa)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    real(kind=dbl_kind) :: thermk    ! canopy gap fraction for thermal IR 
    !                    radiation (-)
    real(kind=dbl_kind) :: tgeff     ! effective skin temp (K)
    ! (takes into effect vegetation, soil
    !       and snow)
    real(kind=dbl_kind) :: thgeff    ! effective skin potential temp (K)
    real(kind=dbl_kind) :: shgeff    ! saturation mixing ratio of effective
    !   skin temperature (kg/kg) 
    real(kind=dbl_kind) :: radt(3)   ! net radiation (W/m^2)
    !   1) canopy leaves
    !   2) ground
    !   3) snow 
    real(kind=dbl_kind) :: p0        ! ground surface precip (after canopy
    !  interception, before snow/rain
    !  partition) (m/sec)
    real(kind=dbl_kind) :: pcpg_rain ! rain fraction of precip reaching ground
    !  (m/sec)
    real(kind=dbl_kind) :: pcpg_snow ! snow fraction of precip reaching ground
    !  (m/sec)
    real(kind=dbl_kind) :: cuprt     ! copy of cupr (m/sec)
    real(kind=dbl_kind) :: lsprt     ! copy of lspr (m/sec)

    real(kind=dbl_kind) :: radc3(2)  ! absorbed radiation (W/m^2)
    ! (1) - radiation absorbed by canopy
    ! (2) - radiation absorbed by ground 

    real(kind=dbl_kind) :: radfac(2,2,2)
    ! radiation absorption factors
    !   (1,1,1)
    !   (1,1,2)
    !   (1,2,1)
    !   (1,2,2)
    !   (2,1,1)
    !   (2,1,2)
    !   (2,2,1)
    !   (2,2,2)

    !...diagnostics/output
    real(kind=dbl_kind) :: hg        ! ground sfc sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: hc        ! canopy sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: hs        ! snow sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: fss       ! CAS-BL sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: fws       ! CAS-BL latent heat flux (W m^-2)
    !...ec,eg, and es are intermediate values used in delef.F and sibslv.F
    real(kind=dbl_kind) :: ec        ! canopy latent heat flux (W m^-2)
    real(kind=dbl_kind) :: eg        ! ground latent heat flux (W m^-2)
    real(kind=dbl_kind) :: es        ! snow   latent heat flux (W m^-2)


    real(kind=dbl_kind) :: egi       ! latent heat flux, ground interception
    !   (puddles) (W m^-2) 
    real(kind=dbl_kind) :: eci       ! latent heat flux, canopy interception
    !   (puddles) (W m^-2) 
    real(kind=dbl_kind) :: egs       ! latent heat flux, ground evaporation
    !    (W m^-2)
    real(kind=dbl_kind) :: ess       ! snow latent heat flux (W m^-2)
    real(kind=dbl_kind) :: ect       ! latent heat flux, canopy transpiration
    !                 (W m^-2)





    real(kind=dbl_kind) :: aparkk    ! canopy PAR use factor (-)
    real(kind=dbl_kind) :: respc(6)  ! canopy respiration (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) sum of physiological types 1-5
    real(kind=dbl_kind) :: pfd       ! incident PAR flux density 
    !           (moles m^-2 sec^-1)
    real(kind=dbl_kind) :: assim(6)  ! gross assimilation (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: assimn(6) ! net assimilation (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: tot_an(13) ! accumulated net assimilation
                                      ! used for calculating respfactor
    real(kind=dbl_kind) :: cflux     ! carbon flux between CAS and reference 
    !   level (mol C  m^-2 sec^-1)
    real(kind=dbl_kind) :: assimnp(6)! assimn scaled by aparkk (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: antemp(6) ! assimn bottom stopped at zero 
    !                  (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: ansqr(6)  ! intermediate squared antemp 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: omepot(6) ! potential light limited photosynthesis
    !                  (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes 
    real(kind=dbl_kind) :: assimpot(6)  
    ! potential biochemical limited 
    !   photosynthesis (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: assimci(6)! potential stress limited 
    !   photosynthesis (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: wsfws(6)  ! intermediate assimpot weighted water 
    !   stress factor (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: wsfht(6)  ! intermediate assimpot weighted high
    !   temperature stress factor (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: wsflt(6)  ! intermediate assimpot weighted low
    !   temperature stress factor (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: wci(6)    ! intermediate antemp weighted 
    !   intercellular CO2 ()
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: whs(6)    ! intermediate antemp weighted 
    !   relative humidity (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: wags(6)   ! intermediate antemp weighted stomatal
    !   conductance ()
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: wegs(6)   ! intermediate evaporation weighted 
    !   stomatal conductance ()
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes


    !   ISOTOPE variables


    real(kind=dbl_kind) :: kiecps(6)  
    ! Kinetic isotope effect during  
    !  photosynthesis by phystype i 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: d13cassimn(6)  
    ! del13C of CO2 assimilated by phystype i 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: c13assimn(6)
    !total flux of C13 in CO2 assimilated 
    !  by phystype i 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: c12assimn(6) 
    !total flux of C12 in CO2 assimilated 
    !  by phystype i
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: rcassimn(6)
    ! isotope ratio (13C/12C) of CO2 
    !   assimilated phystype i (unitless)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    real(kind=dbl_kind) :: flux13c   !
    real(kind=dbl_kind) :: flux12c   !
    real(kind=dbl_kind) :: flux_turb !


end type diagnostic_vars



!------------------------------------------------------------------
!                   STATUS VARIABLES
!------------------------------------------------------------------
type sib_status

    real(kind=dbl_kind) :: coszbar
    real(kind=dbl_kind) :: cosz
    real(kind=dbl_kind) :: dayflag
    real(kind=dbl_kind) :: julday
    integer(kind=int_kind) :: pt_num  ! (-) point number in nSiB vector

end type sib_status


!------------------------------------------------------------------
!                   TOP LEVEL ENCAPSULATING TYPE 'SIB_T'
!------------------------------------------------------------------
type sib_t

    type(param_vars) :: param          ! parameter variables
    type(prognostic_vars) :: prog      ! prognostic variables
    type(diagnostic_vars) :: diag      ! diagnostic variables
    type(sib_status) :: stat            ! sib status variables

end type sib_t




!------------------------------------------------------------------
!                   LOCAL VARIABLES
!------------------------------------------------------------------
type sib_local_vars


    !...canopy air space (CAS) and canopy
    real(kind=dbl_kind) :: dtg       ! delta ground surface temp (K)
    real(kind=dbl_kind) :: dtd(-nsnow+1:nsoil)
    ! delta soil temperature (K)
    real(kind=dbl_kind) :: dtc       ! change in canopy temperature (K)
    real(kind=dbl_kind) :: dts       ! change in snow surface temperature (K)
    real(kind=dbl_kind) :: dth       ! change in ref level temperature (K)
    real(kind=dbl_kind) :: dqm       ! change in ref level moisture (Pa)
    real(kind=dbl_kind) :: dta       ! change in CAS temperature (K)
    real(kind=dbl_kind) :: dea       ! change in CAS moisture (Pa)

    real(kind=dbl_kind) :: etc       ! saturation vapor pressure at Tc (hPa)
    !   ('e-star' of Tc)
    real(kind=dbl_kind) :: getc      ! derivative of etc with respect to temp
    !   (d(etc)/dTc (hPa K^-1)
    real(kind=dbl_kind) :: etg       ! 'e-star' of ground surface (Pa)
    real(kind=dbl_kind) :: getg      ! d(etg)/dTg (hPa K^-1)
    real(kind=dbl_kind) :: ets       ! 'e-star' of snow surface (Pa)
    real(kind=dbl_kind) :: gets      ! d(ets)/dTs (hPa K^-1)

    real(kind=dbl_kind) :: fc        ! direction of vapor flux, canopy to CAS
    !  =1 when flux is canopy to CAS (etc>ea)
    !  =0 when flux is CAS to canopy (ea>etc)
    real(kind=dbl_kind) :: fg        ! direction of vapor flux, ground to CAS
    !  =1 when flux is ground to CAS (etg>ea)
    !  =0 when flux is CAS to ground (ea>etg)
    real(kind=dbl_kind) :: gect      ! dry fraction of veg / rc
    real(kind=dbl_kind) :: geci      ! wet fraction of veg / 2rb
    real(kind=dbl_kind) :: gegs      ! dry fraction of ground / rds
    real(kind=dbl_kind) :: gegi      ! wet fraction of ground /rd
    real(kind=dbl_kind) :: coc       ! gect + geci
    real(kind=dbl_kind) :: cog1      ! gegi + gegs*hrr
    real(kind=dbl_kind) :: cog2      ! gegi + gegs


    !...radiation


    integer(kind=int_kind) :: imelt(-nsnow+1:nsoil)
    ! flag for melting/freezing
    !   1=> melting, 2=> freezing
    real(kind=dbl_kind) :: frac_iceold(-nsnow+1:nsoil)
    ! start-of-timestep value for
    !   ice fraction of total liquid
    real(kind=dbl_kind) :: dtc4      ! d(canopy thermal em)/dT (W m^-2 K^-1)
    real(kind=dbl_kind) :: dtg4      ! d(ground thermal em)/dT (W m^-2 K^-1)
    real(kind=dbl_kind) :: dts4      ! d(snow thermal em)/dT   (W m^-2 K^-1)
    real(kind=dbl_kind) :: lcdtc     ! d(canopy thermal em)/dtc (W m^-2 K^-1)
    real(kind=dbl_kind) :: lcdtg     ! d(canopy thermal em)/dtg (W m^-2 K^-1)
    real(kind=dbl_kind) :: lcdts     ! d(canopy thermal em)/dts (W m^-2 K^-1)
    real(kind=dbl_kind) :: lgdtc     ! d(ground thermal em)/dtc (W m^-2 K^-1)
    real(kind=dbl_kind) :: lgdtg     ! d(ground thermal em)/dtg (W m^-2 K^-1)
    real(kind=dbl_kind) :: lsdts     ! d(snow thermal em)/dts   (W m^-2 K^-1)
    real(kind=dbl_kind) :: lsdtc     ! d(snow thermal em)/dtc   (W m^-2 K^-1)
    real(kind=dbl_kind) :: hcdtc     ! d(canopy H)/dtc  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hcdta     ! d(canopy H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hgdta     ! d(ground H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hgdtg     ! d(ground H)/dtg  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hsdta     ! d(snow H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hsdts     ! d(snow H)/dtsnow  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hadta     ! d(CAS H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hadth     ! d(CAS H)/dtheta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: ecdtc     ! d(canopy LE)/dtc (W m^-2 K^-1)
    real(kind=dbl_kind) :: ecdea     ! d(canopy LE)/dea (W m^-2 K^-1)
    real(kind=dbl_kind) :: egdtg     ! d(ground LE)/dtg (W m^-2 K^-1)
    real(kind=dbl_kind) :: egdea     ! d(ground LE)/dea (W m^-2 K^-1)
    real(kind=dbl_kind) :: esdts     ! d(snow LE)/dtsnow (W m^-2 Pa^-1)
    real(kind=dbl_kind) :: esdea     ! d(snow LE)/dea (W m^-2 K^-1)
    real(kind=dbl_kind) :: eadea     ! d(CAS LE)/dea (W m^-2 Pa^-1)
    real(kind=dbl_kind) :: eadem     ! d(CAS LE)/dem (W m^-2 Pa^-1)
    real(kind=dbl_kind) :: closs     ! canopy thermal loss     (W m^-2)
    real(kind=dbl_kind) :: gloss     ! ground thermal loss     (W m^-2)
    real(kind=dbl_kind) :: sloss     ! snow thermal loss       (W m^-2)
    real(kind=dbl_kind) :: fac1      ! effective ground cover for 
    !   thermal radiation (-)

    real(kind=dbl_kind) :: td_old(-nsnow+1:nsoil)
    ! prev timestep soil temperature (K)

end type sib_local_vars


end module sibtype

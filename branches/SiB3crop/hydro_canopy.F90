
!=====================SUBROUTINE HYDRO_CANOPY===========================

subroutine hydro_canopy(sib,sib_loc)

use kinds
use sibtype

use physical_parameters, only: &
    tice

use sib_const_module, only: &
    cww,    &
    clai,   &
    snomel, &
    denice, &
    denh2o, &
    dtt,    &
    dti,    &
    cpliq,  &
    cpice,  &
    tkice


implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib
type(sib_local_vars)     ,intent(inout) :: sib_loc
! variables local to SiB

!----------------------------------------------------------------------  


!=======================================================================
!
!     CALCULATION OF  INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW
!     INCORPORATING EFFECTS OF PATCHY SNOW COVER AND TEMPERATURE
!     ADJUSTMENTS.
!
!----------------------------------------------------------------------
!
!     (1) NON-UNIFORM PRECIPITATION
!         CONVECTIVE PPN. IS DESCRIBED BY AREA-INTENSITY
!         RELATIONSHIP :-
!
!                   F(X) = A*EXP(-B*X)+C
!
!         THROUGHFALL, INTERCEPTION AND INFILTRATION
!         EXCESS ARE FUNCTIONAL ON THIS RELATIONSHIP
!         AND PROPORTION OF LARGE-SCALE PPN.
!         REFERENCE: SA-89B, APPENDIX.
!
!
!
!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       TC             CANOPY TEMPERATURE (K)
!       CAPAC(1)       CANOPY LIQUID INTERCEPTION STORE (kg m^-2)
!       SNOW_VEG       CANOPY SNOW INTERCEPTION STORE (Kg M^-2) or (mm water)
!       P0             THROUGHFALL PRECIPITATION (M/SEC)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     local variables
real(kind=dbl_kind) :: pcoefs(2,2)
real(kind=dbl_kind) :: ap
real(kind=dbl_kind) :: bp
real(kind=dbl_kind) :: cp
real(kind=dbl_kind) :: totalp    ! total precip (meters)
real(kind=dbl_kind) :: pinf
real(kind=dbl_kind) :: fpi
real(kind=dbl_kind) :: xsc   ! canopy excess water (kg/m^2)
real(kind=dbl_kind) :: xss   ! canopy excess snow (kg/m^2)
real(kind=dbl_kind) :: xs
real(kind=dbl_kind) :: capacp    ! copy of veg liquid store
real(kind=dbl_kind) :: snowwp    ! copy of veg snow
real(kind=dbl_kind) :: spechc    ! specific heat of canopy (intercepted
!  water and vegatation )(J m^-2 deg^-1)
real(kind=dbl_kind) :: chiv      ! leaf angle dist factor (unitless)
real(kind=dbl_kind) :: aa
real(kind=dbl_kind) :: bb
real(kind=dbl_kind) :: exrain
real(kind=dbl_kind) :: zload
real(kind=dbl_kind) :: tti       ! direct throughfall in meters
real(kind=dbl_kind) :: tex       ! canopy drainage in meters
real(kind=dbl_kind) :: thru      ! total throughfall (tti + tex)
real(kind=dbl_kind) :: freeze
real(kind=dbl_kind) :: diff
real(kind=dbl_kind) :: ccp
real(kind=dbl_kind) :: cct
real(kind=dbl_kind) :: tsd       ! temperature modified for canopy 
! interception
real(kind=dbl_kind) :: tta
real(kind=dbl_kind) :: ttb
real(kind=dbl_kind) :: cca
real(kind=dbl_kind) :: ccb
real(kind=dbl_kind) :: ccc
real(kind=dbl_kind) :: arg
real(kind=dbl_kind) :: fliq      ! fraction of liquid in precip (0-1)
real(kind=dbl_kind) :: bifall    ! density of falling snow (kg/m^-3)
real(kind=dbl_kind) :: dz_snowfall ! new snow depth (m)
real(kind=dbl_kind) :: capac1m   ! capac(1) ( canopy water in meters)

integer(kind=int_kind) :: pcptype    ! precip type; 1=rain 2=snow
integer(kind=int_kind) :: newnode    ! new snow layer indicator
integer(kind=int_kind) :: i,j        ! loop index

data pcoefs(1,1)/ 20. /, pcoefs(1,2)/ .206E-8 /, &
    pcoefs(2,1)/ 0.0001 /, pcoefs(2,2)/ 0.9999 /, bp /20. /

    !-----------------------------------------------------------------------
    !
    !     PREC ( PI-X )   : EQUATION (C.3), SA-89B
    !
    !-----------------------------------------------------------------------
    dz_snowfall = 0.0

    ap = pcoefs(2,1)
    cp = pcoefs(2,2)
    totalp  = (sib%diag%cuprt + sib%diag%lsprt) *  dtt

    !itb...check against ridiculously low precip amounts...
    if(totalp < 1.0E-10) then
        totalp = 0.0
        sib%diag%cuprt = 0.0
        sib%diag%lsprt = 0.0
    endif

    if( sib%prog%snow_veg >  0.0 .or. sib%prog%snow_mass >  0.0 &
        .or. sib%prog%tm < tice ) sib%diag%cuprt = 0.0
    sib%diag%lsprt = totalp/dtt - sib%diag%cuprt

    if(totalp > 1.e-8) then
        ap = sib%diag%cuprt * dtt / totalp * pcoefs(1,1) + &
            sib%diag%lsprt * dtt / totalp * pcoefs(2,1)
        cp = sib%diag%cuprt * dtt / totalp * pcoefs(1,2) + &
            sib%diag%lsprt * dtt / totalp * pcoefs(2,2)
    endif

    thru = 0.
    fpi  = 0.

    !----------------------------------------------------------------------
    !     PRECIP INPUT INTO HYDRO_CANOPY IN M/SEC; TOTALP IS IN METERS
    !----------------------------------------------------------------------

    sib%diag%p0 = totalp  !sib%diag%p0 now in meters

    !itb...calculate capac1m
    capac1m = sib%prog%capac(1)/denh2o

    xsc = max(0.0_dbl_kind, capac1m - sib%param%satcap(1) )
    capac1m = capac1m - xsc
    xss = max(0.0_dbl_kind, (sib%prog%snow_veg/denice - sib%param%satcap(1)) ) 
    sib%prog%snow_veg = (sib%prog%snow_veg/denice - xss)*denice

    !itb...add excess to total throughfall
    thru = thru + xsc + xss

    capacp = capac1m                !m of liquid on veg
    snowwp = sib%prog%snow_veg/denice    !m snow on veg

    !...capac + snow_veg needs to be in meters for units to work here.
    !...units of spechc will be J m^-2 deg^-1

    spechc = &
        min( 0.05_dbl_kind, ( capac1m + sib%prog%snow_veg/denice ) ) &
        * cww + sib%param%zlt * clai

    !----------------------------------------------------------------------
    !    PROPORTIONAL SATURATED AREA (XS) AND LEAF DRAINAGE(TEX)
    !
    !     TTI ( D-D )     : EQUATION (C.4), SA-89B
    !     XS  ( X-S )     : EQUATION (C.7), SA-89B
    !     TEX ( D-C )     : EQUATION (C.8), SA-89B
    !
    !    SA-89B is Sato et al, Implementing the Simple Biosphere Model 
    !    (SiB) in a General Circulation Model: Methodologies and Results
    !    NASA Contractor Report 185509 (1989)
    !
    !-----------------------------------------------------------------------

    chiv = sib%param%chil
    if ( abs(chiv) <= 0.01 ) chiv = 0.01

    aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv   ! unitless
    bb = 0.877 * ( 1. - 2. * aa )                  ! unitless
    exrain = aa + bb                               ! unitless

    zload = capac1m + sib%prog%snow_veg/denice            ! meters
    fpi   = ( 1.-exp( - exrain*sib%param%zlt/sib%param%vcover ) )* sib%param%vcover

    !...tti is direct throughfall (meters)

    tti   = sib%diag%p0 * ( 1.-fpi )

    !...xs is fraction of canopy where intercepted plus existing storage
    !...exceeds capacity (unitless) 

    xs    = 1.

    if ( sib%diag%p0 >= 1.e-9 ) then  

        arg = ( sib%param%satcap(1)-zload )/ ( sib%diag%p0*fpi*ap ) -cp/ap 

        if ( arg >= 1.e-9 ) then                                 
            xs = -1./bp * log( arg )                                      
            xs = min( xs, 1.0_dbl_kind ) 
            xs = max( xs, 0.0_dbl_kind ) 
        endif   
    endif 

    !...tex is canopy drainage (meters) 

    tex = sib%diag%p0*fpi * &
        ( ap/bp*( 1.- exp( -bp*xs )) + cp*xs ) - &
        ( sib%param%satcap(1) - zload ) * xs                        
    tex = max( tex, 0.0_dbl_kind ) 

    !-------------------------------------------------------------
    !    TOTAL THROUGHFALL (THRU) AND STORE AUGMENTATION         
    !-----------------------------------------------------------

    !...thru is throughfall (meters) 

    thru = thru + tti + tex 

    !...pinf is rainfall intercepted by canopy (meters)
    pinf = sib%diag%p0 - thru

    !...make sure interception is not negative
    pinf = MAX(pinf,0.0_dbl_kind)

    if( sib%prog%tm > tice ) then
        capac1m = capac1m + pinf
    else
        sib%prog%snow_veg = (sib%prog%snow_veg/denice + pinf)*denice 
    endif       

    !itb...this is the old sib routine, 'ADJUST'
    !=======================================================================
    !
    !     TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION
    !
    !=======================================================================

    freeze = 0.
    diff = ( capac1m+sib%prog%snow_veg/denice - capacp - snowwp )*cww
 
    ccp = spechc
    cct = spechc + diff

    tsd = ( sib%prog%tc * ccp + sib%prog%tm * diff ) / cct

    if ( ( sib%prog%tc >  tice .AND. sib%prog%tm <= tice ) .or. &
        ( sib%prog%tc <= tice .AND. sib%prog%tm >  tice ) )then

        tta = sib%prog%tc
        ttb = sib%prog%tm
        cca = ccp
        ccb = diff
        if ( tsd <= tice ) then

            !----------------------------------------------------------------
            !    FREEZING OF WATER ON CANOPY 
            !----------------------------------------------------------------

            ccc = capacp * 0.001 * snomel
            if ( sib%prog%tc < sib%prog%tm ) ccc = diff * snomel / cww
            tsd = ( tta * cca + ttb * ccb + ccc ) / cct

            freeze = ( tice * cct - ( tta * cca + ttb * ccb ) )
            freeze = (min ( ccc, freeze )) / snomel
            if(tsd > tice) tsd = tice - 0.01


        else

            !----------------------------------------------------------------
            !    MELTING OF SNOW ON CANOPY 
            !----------------------------------------------------------------

            ccc = - sib%prog%snow_veg/denice * snomel 
            if ( sib%prog%tc > sib%prog%tm ) ccc = - diff * snomel / cww

            tsd = ( tta * cca + ttb * ccb + ccc ) / cct

            freeze = ( tice * cct - ( tta * cca + ttb * ccb ) )
            freeze = (max( ccc, freeze )) / snomel
            if(tsd <= tice)tsd = tice - 0.01
           

        endif
    endif

    sib%prog%snow_veg = (sib%prog%snow_veg/denice + freeze)*denice
    capac1m = capac1m - freeze
    sib%prog%snow_veg = max(sib%prog%snow_veg,0.0_dbl_kind)
    capac1m = max(capac1m,0.0_dbl_kind)

    xs = max( 0.0_dbl_kind, ( capac1m - sib%param%satcap(1) ) )
    if ( sib%prog%snow_veg/denice >= 1.0e-7_dbl_kind ) xs = capac1m
    capac1m = capac1m - xs
    sib%prog%tc = tsd

    !...end of 'ADJUST' code


    !...now assign to p0 (precip) that precipitation that was
    !...NOT intercepted and stored on the canopy

    sib%diag%p0 = thru + xs 

    !itb...change sib%diag%p0 back to rate (mm/sec)
    sib%diag%p0 = sib%diag%p0 * dti * 1000.0  ! conversion from meters to mm/sec 


    !...precipitation onto ground (follows clm_hydro_canopy)

    !...determine precip type

    if(sib%prog%tm >= tice + 2.5) then
        pcptype = 1
    else
        pcptype = 2
    endif

    !...percentage of liquid water by mass is arbitrarily set to vary 
    !...linearly with air temp, from 0% at tice to 40% max at (tice + 2.0)

    if (pcptype == 1 ) then  ! RAIN
        fliq = 1.0
        sib%diag%pcpg_snow = 0.0
        sib%diag%pcpg_rain = sib%diag%p0

    else
        if(sib%prog%tm <= tice) then
            fliq = 0.0
        elseif(sib%prog%tm < tice + 2.0) then
            fliq = -54.61 + 0.2*sib%prog%tm
        else
            fliq = 0.40
        endif
        fliq = max(0.0_dbl_kind,fliq)


        !...use Alta Relationship, Anderson (1976); LaChapelle (1961),
        !...U.S. Department of Agriculture Forest Service, Project F,
        !...Progress Report 1, Alta Avalanche Study Center: Snow
        !...Layer Densification

        if(sib%prog%tm > tice + 2.0) then
            bifall = 189.0
        elseif(sib%prog%tm > tice - 15.0) then
            bifall = 50.0 + 1.7 * (sib%prog%tm - tice + 15.0)**1.5
        else
            bifall = 50.0
        endif

        sib%diag%pcpg_snow = sib%diag%p0 * (1.0 - fliq)
        sib%diag%pcpg_rain = sib%diag%p0 * fliq

        !...snowfall rate will be in m/sec
        !...BE CAREFUL HERE; convervsion units of kg/m^3 for bifall with 
        !...                 kg/m^2/sec for precip (mm/sec) results in
        !...                 output of m/sec snow accumulation.

        dz_snowfall = sib%diag%pcpg_snow/bifall 

        !...snow_depth will be in meters...
        sib%prog%snow_depth = sib%prog%snow_depth + dz_snowfall * dtt

        !...add snow depth to top snow layer
        if(dz_snowfall > 0.0 .and. sib%prog%nsl < 0) then
            sib%prog%dz(sib%prog%nsl+1) = sib%prog%dz(sib%prog%nsl+1) + &
                dz_snowfall * dtt
        endif

        sib%prog%snow_mass  = sib%prog%snow_mass + sib%diag%pcpg_snow * dtt

    endif

    !...initialize snow layer if accumulation exceeds 10mm..


    newnode = 0
    if(   sib%prog%nsl == 0 .and. sib%diag%pcpg_snow > 0.0 ) then 
        newnode                = 1
        sib%prog%nsl           = -1
        sib%prog%dz(0)         = sib%prog%snow_depth
        sib%prog%node_z(0)     = -0.5 * sib%prog%dz(0)
        sib%prog%layer_z(-1)   = -sib%prog%dz(0)
        sib%prog%layer_z(0)    = 0.0
        sib%prog%snow_age      = 0.0
        sib%prog%td(0)         = MIN(tice, sib%prog%ta)
        sib%prog%www_ice(0)    = sib%prog%snow_mass
        sib%prog%www_liq(0)    = 0.0
        sib_loc%frac_iceold(0) = 1.0
        sib%param%shcap(0) = cpliq*sib%prog%www_liq(0) + &
            cpice*sib%prog%www_ice(0)
        !--itb this is a patch
        sib%param%tksoil(0)       = tkice
    endif

    !itb...new snowfall is added here
    if(sib%prog%nsl < 0 .and. newnode == 0) then
        sib%prog%www_ice(sib%prog%nsl+1) = sib%prog%www_ice(sib%prog%nsl+1) + &
            dtt * sib%diag%pcpg_snow
        sib%prog%dz(sib%prog%nsl+1) = sib%prog%dz(sib%prog%nsl+1) + &
            dz_snowfall * dtt
    endif

    sib%prog%capac(1) = capac1m * denh2o


end subroutine hydro_canopy


!=================SUBROUTINE UPDATE=====================

subroutine update(sib,sib_loc)

!itb===============================================
!itb
!itb  MOVING STORAGE TERM UPDATES (SNOW, CAPAC) HERE FROM ENDTEM, WHICH
!itb     NO LONGER EXISTS. FLUXES PREVIOUSLY CALCULATED IN ENDTEM ARE
!itb     TAKEN CARE OF IN THE PROGNOSTIC C.A.S. CALCULATIONS, SO WE
!itb     MERELY NEED TO TAKE CARE OF STORAGE TERMS NOW.
!itb
!itb      November 2000
!
!=================================================
!
!     UPDATING OF ALL HYDROLOGICAL PROGNOSTIC VARIABLES.  SNOW AND
!        RUNOFF CALCULATIONS (SEE ALSO INTER2).  SUBROUTINES SNOW2 AND
!        RUN2 OF 1D MODEL ARE INLINED IN THIS CODE.
!=================================================


!+++++++++++++++++++++++OUTPUT++++++++++++++++++++++++
!
!       DTC            CANOPY TEMPERATURE INCREMENT (K)
!       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
!       WWW(3)         GROUND WETNESS 
!       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (kg m^-2)
!       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
!       ROFF           RUNOFF (MM/HR)
!
!++++++++++++++++++++++++++DIAGNOSTICS+++++++++++++++++
!
!       ECMASS         CANOPY EVAPOTRANSPIRATION (MM)
!       EGMASS         GROUND EVAPOTRANSPIRATION (MM)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++


use kinds
use sibtype
use sib_const_module, only:  &
    snofac,  &
    nsoil,   &
    snomel,  &
    denh2o,  &
    denice,  &
    phmin,   &
    dtt,     &
    dti,     &
    tkair,   &
    tkice,   &
    cpliq,   &
    cpice

use physical_parameters, only:  &
    hltm,                &
    pi,                  &
    tice,                &
    cp => spec_heat_cp

use eau_params, only : &
    lfus


implicit none

!-----------------------------------------------------------------

type(sib_t), intent(inout) :: sib
type(sib_local_vars)     ,intent(inout) :: sib_loc
  ! variables local to SiB

!------------------------------------------------------------------  



!...LOCAL VARIABLES
real(kind=dbl_kind) :: flux_imb(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: water_imb(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: wice0(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: wliq0(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: wmass(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: sflux(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: sflx1(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: htstor(sib%prog%nsl+1:nsoil)
real(kind=dbl_kind) :: rsnow  ! fraction of snow that is ice
real(kind=dbl_kind) :: ecpot
real(kind=dbl_kind) :: egpot
real(kind=dbl_kind) :: espot
real(kind=dbl_kind) :: hrr
real(kind=dbl_kind) :: cogs1
real(kind=dbl_kind) :: cogs2
real(kind=dbl_kind) :: ectmax(nsoil)  ! layer max transpiration (J m^-2)
real(kind=dbl_kind) :: ectmax_tot !column total max transp (J m^-2)
real(kind=dbl_kind) :: egsmax
real(kind=dbl_kind) :: facks
real(kind=dbl_kind) :: cpdpsy
real(kind=dbl_kind) :: timcon
real(kind=dbl_kind) :: hltmi
real(kind=dbl_kind) :: hfus
real(kind=dbl_kind) :: ecidif
real(kind=dbl_kind) :: egidif
real(kind=dbl_kind) :: egit
real(kind=dbl_kind) :: temp
real(kind=dbl_kind) :: propor
real(kind=dbl_kind) :: heatr
real(kind=dbl_kind) :: le_phase

real(kind=dbl_kind) :: btran
real(kind=dbl_kind) :: tinc(sib%prog%nsl+1:nsoil) ! temperature increment
real(kind=dbl_kind) :: shcap_temp ! new heat capacity of updated top snow level 

real(kind=dbl_kind), dimension(-nsnow+1:nsoil) :: &
        thk ! thermal conductivity of layer (W/m/K)

real(kind=dbl_kind) :: bw   ! partial density of water (ice+liquid)

integer(kind=int_kind) :: i,j

!------------------------------------------------------------------
!
!------------------------------------------------------------------


    timcon = pi / 86400.0

    hltmi = 1. / hltm      ! inverse heat of vap, units kg/J

    hfus = snomel/1000.0   ! heat of fusion, units J/kg

    cpdpsy = cp/sib%diag%psy

    if(sib%prog%nsl < 0) then
        wice0(1) = 0.0
        wliq0(1) = 0.0
        do i=sib%prog%nsl+1,0,-1
            wice0(1) = wice0(1) + sib%prog%www_ice(i)
            wliq0(1) = wliq0(1) + sib%prog%www_liq(i)
        enddo
        !itb...patch
        if (wliq0(1) > 0.0 ) then
            rsnow = wice0(1) / (wliq0(1) + sib%prog%capac(2) )
        else
            rsnow = 1.0
        endif
        wice0(1) = 0.0
        wliq0(1) = 0.0
    else
        rsnow = 0.0
    endif

    !itb...effective porosity
    do i = sib%prog%nsl+1,nsoil
        sib%prog%vol_ice(i) = min(sib%param%poros, sib%prog%www_ice(i)/(sib%prog%dz(i)*denice))
        sib%diag%eff_poros(i) = sib%param%poros - sib%prog%vol_ice(i)
        sib%prog%vol_liq(i) = min(sib%diag%eff_poros(i),  &
            sib%prog%www_liq(i)/(sib%prog%dz(i)*denh2o))
        !if no effective porosity, lose liquid to runoff
        !itb...this is to maintain water balance...
        if (sib%prog%vol_liq(i) == 0.0 .and. sib%prog%www_liq(i) > 0.0 ) then
            sib%diag%roff = sib%diag%roff + sib%prog%www_liq(i)
            sib%prog%www_liq(i) = 0.0
        endif
    enddo

!itb...total soil water
    sib%diag%www_tot_soil = 0.0
    do i=1,nsoil
     sib%diag%www_tot_soil = sib%diag%www_tot_soil +      &
               sib%prog%www_liq(i) + sib%prog%www_ice(i)
    enddo



    !pl this is the potential gradient in Pa
    !pl this WAS realized in sibslv

    ecpot = (sib_loc%etc + sib_loc%getc*sib_loc%dtc) - sib%prog%ea 

    !pl and this is the  INTERCEPTION flux in J/m2

    sib%diag%eci = ecpot * sib_loc%geci * sib%prog%ros * cpdpsy * dtt

    !pl and this is the TRANSPIRATION flux in J/m2

    sib%diag%ect = ecpot * sib_loc%gect * sib%prog%ros * cpdpsy * dtt


    !itb...ground interception terms will be split up further in
    !itb...the new code; previously, mean values were used.

    !pl this is the potential gradient in Pa (ground to CAS and snow to CAS)
    !pl this WAS realized in sibslv

    if(sib%prog%nsl == 0 ) then
        egpot = (sib_loc%etg + sib_loc%getg*sib_loc%dtg) - sib%prog%ea 
        espot = 0.0

    else

        !itb...i'm going to gamble that ESS will generally be low, and not
        !itb...larger than the available amount of snow. This is also because
        !itb...small snow amounts are folded into top soil layer by subroutine
        !itb...combine_snow

        egpot = 0.0
        espot = (sib_loc%ets  + sib_loc%gets* sib_loc%dts) - sib%prog%ea  


        !itb...snow surface evaporation

        sib%diag%ess = espot * sib%prog%ros * cpdpsy / sib%diag%rd * snofac 

        !itb...HARDWIRE-MAKE SNOW EVAP ZERO
        sib%diag%ess = 0.0


        !itb...??? when is sib%diag%ess subtracted from snow amount?
        !itb...is it taken from snow water or snow ice?
!            sib%prog%www_liq(sib%prog%nsl+1)*hltm*dti,              &
!            sib%prog%www_ice(sib%prog%nsl+1)/snofac*hltm*dti

        !itb...can't evaporate more snow than you have...
        sib%diag%ess = min(sib%diag%ess,sib%prog%www_ice(sib%prog%nsl+1)/  &
            snofac*hltm*dti)

        !itb...currently keeping ESS positive
        sib%diag%ess = max(sib%diag%ess,0.0_dbl_kind)

        !itb...now subtract it out
        sib%prog%www_ice(sib%prog%nsl+1) = sib%prog%www_ice(sib%prog%nsl+1)  &
            /snofac*hltm*dti - sib%diag%ess


        sib%prog%www_ice(sib%prog%nsl+1) = sib%prog%www_ice(sib%prog%nsl+1)  &
            * snofac * dtt /hltm

        if(sib%prog%www_ice(sib%prog%nsl+1) == 0.0 ) then
            if(sib%prog%www_liq(sib%prog%nsl+1) > 0.0) then
                sib%prog%capac(2) = sib%prog%capac(2) +  &
                    sib%prog%www_liq(sib%prog%nsl+1)
                sib%prog%www_liq(sib%prog%nsl+1) = 0.0
                sib%prog%vol_liq(sib%prog%nsl+1) = 0.0
                sib%prog%td(sib%prog%nsl+1) = 0.0
            endif                   ! get rid of snow layer
            sib%prog%dz(sib%prog%nsl+1)     = 0.0
            sib%prog%node_z(sib%prog%nsl+1) = 0.0
            sib%prog%layer_z(sib%prog%nsl)  = 0.0
            sib%prog%nsl = sib%prog%nsl + 1  
        endif

    endif



    !pl and this is the  INTERCEPTION flux in J/m2
!      sib%diag%egi = egpot * sib_loc%gegi * (1.-sib%diag%areas) * sib%prog%ros           &
!            * cpdpsy * dtt + espot * sib%diag%areas/sib%diag%rd             &
!            * sib%prog%ros * cpdpsy * dtt
    !itb...keeping things simple-not dealing with areas for now...
!      sib%diag%egi = egpot * sib_loc%gegi * (1.-sib%diag%areas) * sib%prog%ros           &
!            * cpdpsy * dtt  

    sib%diag%egi = egpot * sib_loc%gegi * sib%prog%ros * cpdpsy * dtt

    hrr = sib%diag%hr
    if ( sib_loc%fg < 0.5_dbl_kind ) hrr = 1.

    cogs1 =  sib_loc%gegs * hrr * (1.-sib%diag%areas)                
    cogs2 =  sib_loc%gegs       * (1.-sib%diag%areas)

    !pl and this is the EVAPORATIVE flux in J/m2
    sib%diag%egs =  (sib_loc%etg + sib_loc%getg * sib_loc%dtg) * cogs1  & 
        -sib%prog%ea * cogs2

    sib%diag%egs = sib%diag%egs * sib%prog%ros * cpdpsy *dtt


    !itb...make sure you don't evap more than you have (units are J m^-2)
    !itb...SiB2 restricted evaporation by limiting egs to 1/2 the water in 
    !itb...soil moisture layer 1. 
    !itb...In SiB3, we will restrict by wilting point; WP is a soil  
    !itb...parameter, currently must be set by user. We're working on it...

    egsmax = 0.5_dbl_kind * sib%prog%www_liq(1)
    egsmax = egsmax * hltm
    egsmax = max(egsmax,0.0_dbl_kind)

    sib%diag%egs = min( sib%diag%egs, egsmax )

    !itb...Make sure you don't transpire more water than is in the soil.
    !itb...Have to play around with this; both evaporation and transpiration
    !itb...can remove water from top soil layer. Plan is to give evaporation
    !itb...priority in top layer, then push any transpiration deficit downwards
    !itb...in the soil. Idea is that lower roots will supply water (if
    !itb...available) to make up deficit from upper soil. 

    ectmax(:)  = 0.0
    ectmax_tot = 0.0

    do i=1,nsoil
        ectmax(i) =  sib%diag%paw(i) * denh2o * sib%prog%dz(i) * hltm
        ectmax_tot = ectmax_tot + ectmax(i)
    enddo


    !itb...cannot transpire more than PAW

    sib%diag%ect = min ( ectmax_tot, sib%diag%ect )

    !itb...do not allow negative transpiration...

    sib%diag%ect = max (sib%diag%ect, 0.0_dbl_kind)


    !itb...need to do a check for positive values...
    !itb...FOR NOW, not allowing negative egs (dew or frost) accumulate
    !itb...on ground-have had problems with this in the past.

    btran       = 1.0E-10


    if(sib%diag%egs > 0.0 )then
       if(sib%prog%www_liq(1) > 1.0E-10) then
        if( (sib%diag%egs + sib%diag%ect*sib%param%rootf(1)) > egsmax ) then
             sib%param%rootr(1) = (0.5_dbl_kind * sib%prog%www_liq(1) * hltm -      &
                                                    sib%diag%egs)/sib%diag%ect
            sib%param%rootr(1) = MIN(sib%param%rootr(1), 1.0_dbl_kind)
            sib%param%rootr(1) = MAX(sib%param%rootr(1), 0.0_dbl_kind)
            sib%param%rootr(1) = sib%param%rootr(1) * sib%param%rootf(1)
            if(sib%prog%www_liq(1) <= sib%param%vwcmin) sib%param%rootr(1) = 0.0_dbl_kind
            
        else
            sib%param%rootr(1) = (1.0 - sib%param%vwcmin/sib%prog%vol_liq(1)) /  &
                (1.0 - sib%param%vwcmin/sib%param%fieldcap)
            sib%param%rootr(1) = MIN(sib%param%rootr(1), 1.0_dbl_kind)
            sib%param%rootr(1) = MAX(sib%param%rootr(1), 0.0_dbl_kind)
            sib%param%rootr(1) = sib%param%rootr(1) * sib%param%rootf(1)
        endif
      else
        sib%param%rootr(1) = 0.0
      endif
    else
        if(sib%prog%www_liq(1) > 1.0E-10) then
            sib%param%rootr(1) = (1.0 - sib%param%vwcmin/sib%prog%vol_liq(1)) /  &
                (1.0 - sib%param%vwcmin/sib%param%fieldcap)
            sib%param%rootr(1) = MIN(sib%param%rootr(1), 1.0_dbl_kind)
            sib%param%rootr(1) = MAX(sib%param%rootr(1), 0.0_dbl_kind)
            sib%param%rootr(1) = sib%param%rootr(1) * sib%param%rootf(1)
        else
            sib%param%rootr(1) = 0.0
        endif
    endif

    

    !itb...no transpiration loss from frozen soil...
    if(sib%prog%td(1) < tice) sib%param%rootr(1) = 0.0
    btran = btran + sib%param%rootr(1)


    do i=2,nsoil
        if(sib%prog%td(i) >= tice .and. sib%prog%vol_liq(i) > 0.0 ) then
            sib%param%rootr(i) = (1.0 - sib%param%vwcmin/sib%prog%vol_liq(i)) /  &
                (1.0 - sib%param%vwcmin/sib%param%fieldcap)

            sib%param%rootr(i) = MIN(sib%param%rootr(i), 1.0_dbl_kind)
            sib%param%rootr(i) = MAX(sib%param%rootr(i), 0.0_dbl_kind)

            sib%param%rootr(i) = sib%param%rootr(i) * sib%param%rootf(i)
            btran = btran + sib%param%rootr(i)           
        else
            sib%param%rootr(i) = 0.0
        endif
    enddo 

    !itb...normalize to get layer contribution
    do i=1,nsoil
        sib%param%rootr(i) = sib%param%rootr(i) / btran
        sib%param%rootr(i) = max(sib%param%rootr(i), 0.0_dbl_kind)
    enddo
!
    !itb...any deficit left over will have to be partitioned into sensible heat


    !itb...these fluxes were all realized in sibslv. If positive, they
    !itb...imply transfer of water vapor INTO the CAS. If negative,
    !itb...they imply transfer OUT OF the CAS. We need to adjust
    !itb...the various reserviors as well as the CAS vapor capacity, 
    !itb...making sure that none go to negative values.

    !itb...the actual movement of vapor is taken care of in the
    !itb...equations in sibslv. All we do now is adjust the surface
    !itb...and vegetation storage reservoirs to reflect what we've
    !itb...already added or taken out.

    !pl this is the limitation to the ECI flux in J/m2

    ecidif=MAX(0.0_dbl_kind,(sib%diag%eci - sib%prog%snow_veg * hltm -      &
        sib%prog%capac(1)  * hltm))

    sib%diag%eci   =MIN(sib%diag%eci,                                       &
        ( sib%prog%snow_veg*hltm + sib%prog%capac(1) * hltm))

    !pl this is the EGI flux in J/m2

    egidif=                                                       &
        MAX(0.0_dbl_kind,sib%diag%egi-(sib%prog%snow_mass +                  &
        sib%prog%capac(2)) * hltm) * (1.-rsnow)

    egit  =                                                       &
        MIN(sib%diag%egi, ((sib%prog%snow_mass +                             &
        sib%prog%capac(2)) * hltm + sib%diag%ess) ) * (1.-rsnow)

    !itb...this is a little less complicated than the old technique...
    sib%diag%egi = egit

    !...need to reduce ground surface store (ponding) by sib%diag%egi...

    sib%prog%capac(2) = sib%prog%capac(2) - sib%diag%egi * (1.0-rsnow) * hltmi
!------------------------------------------------------------------
!     CALCULATION OF SENSIBLE HEAT FLUXES FOR THE END OF THE TIMESTEP.
!        SEE FIGURE (2) OF SE-86.  NOTE THAT INTERCEPTION LOSS EXCESS
!        ENERGIES (ECIDIF, EGIDIF) ARE ADDED.
!
!      HC          (HC)    : EQUATION (63) , SE-86
!      HG          (HGS)   : EQUATION (65) , SE-86
!------------------------------------------------------------------

    !itb...i've left the leaf one-sided, for now...

    sib%diag%hc = ( sib%prog%tc  - sib%prog%ta)  /sib%diag%rb       &
        * sib%prog%ros * cp * dtt + ecidif 

    sib%diag%chf = sib%param%czc * dti * sib_loc%dtc

    !itb...keep in mind that currently we are not worrying about 
    !itb...partial snow cover issues. with that in mind, there will
    !itb...only exist ground H when the is NO snow, snow H when
    !itb...there IS snow.

    !itb...ground sensible heat flux 

    !itb...snow fluxes
    if(sib%prog%nsl < 0 ) then  !SNOW case

        !...remember that td was updated in addinc
        sib%diag%hs = (sib%prog%td(sib%prog%nsl+1)  - sib%prog%ta )    &
            / sib%diag%rd * sib%prog%ros * cp * dtt
        sib%diag%hg = 0.0
    else
        sib%diag%hg = ( sib%prog%tg - sib%prog%ta ) /sib%diag%rd       &
            * sib%prog%ros * cp * dtt + egidif
        sib%diag%hs = 0.0
    endif

!-------------------------------------------------------------------
!     CALCULATION OF STORAGE HEAT FLUXES
!------------------------------------------------------------------ 


    !itb...this ugly beast represents G, soil heat flux. G is calculated
    !itb...as normalized (snow and non-snow) sum of change in storage in 
    !itb...top layer and flux from first layer downwards.

    shcap_temp = cpliq * sib%prog%www_liq(sib%prog%nsl+1) +                    &
                 cpice * sib%prog%www_ice(sib%prog%nsl+1)

    !...thermal conductivity of snow
    if (sib%prog%nsl < 0) then
        do j= sib%prog%nsl+1,0
            bw = (sib%prog%www_liq(j) + sib%prog%www_ice(j))/sib%prog%dz(j)
            thk(j) = tkair + (7.75E-5*bw + 1.105E-6*bw*bw)*(tkice-tkair) 
        enddo
    endif

    do j = sib%prog%nsl+1,0

        sib%param%tksoil(j) = thk(j)*thk(j+1)* &
            (sib%prog%node_z(j+1)-sib%prog%node_z(j))  &
            /(thk(j)*(sib%prog%node_z(j+1)-sib%prog%layer_z(j)) + &
            thk(j+1)*(sib%prog%layer_z(j)-sib%prog%node_z(j)))

    enddo

    do j=sib%prog%nsl+1,0
        sib%param%slamda(j) = sib%param%tksoil(j) / (sib%prog%node_z(j+1) - &
            sib%prog%node_z(j))
    enddo



    sib%diag%shf =  &
        (1.0 - sib%diag%areas) * ((sib_loc%dtd(1) * sib%param%shcap(1) * dti *  &
        sib%prog%dz(1))   +                                                  &
        
        (sib%param%slamda(1) * (sib%prog%td(1) - sib%prog%td(2) )) )   +        &

        sib%diag%areas * ((sib_loc%dtd(sib%prog%nsl+1) * shcap_temp *        &
        dti * sib%prog%dz(sib%prog%nsl+1)) +                                 &

        (sib%param%tksoil(sib%prog%nsl+1) * (sib%prog%td(sib%prog%nsl+2) -      &
        sib%prog%td(sib%prog%nsl+1))/(sib%prog%node_z(sib%prog%nsl+2) -      &
        sib%prog%node_z(sib%prog%nsl+1)))  )



!-----------------------------------------------------------------
!    INTERCEPTION LOSSES APPLIED TO SURFACE WATER STORES.                      
!    EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY
!    ( HLTM*1000.) LOSS IS IN M M-2. MASS TERMS ARE IN KG M-2 DT-1
!    INTERCEPTION AND DRAINAGE TREATED IN INTER2.
!
!      CAPAC/SNOWW(1) (M-C)   : EQUATION (3)  , SE-86
!      CAPAC/SNOWW(2) (M-G)   : EQUATION (4)  , SE-86
!------------------------------------------------------------------

    !PL HERE WE DO A CHECK FOR CONDENSATION AND MAKE SURE THAT IT ONLY
    !PL HAPPENS TRHOUGH ECI AND EGI

    rsnow = (sib%prog%snow_veg/denice)/   &
        (sib%prog%capac(1)/denh2o + sib%prog%snow_veg/denice + 1.0E-10)

    facks = 1. + rsnow * ( snofac-1. )

    if ( (sib%diag%ect+sib%diag%eci) <= 0.0) then
        sib%diag%eci = sib%diag%ect+sib%diag%eci
        sib%diag%ect = 0.
        facks = 1. / facks
    endif

    sib%prog%capac(1) = sib%prog%capac(1)-(1.-rsnow)*sib%diag%eci*facks*hltmi
    sib%prog%snow_veg = sib%prog%snow_veg - rsnow*sib%diag%eci * facks * hltmi
    sib%prog%snow_veg = max(sib%prog%snow_veg,0.0_dbl_kind)
    sib%prog%capac(1) = max(sib%prog%capac(1),0.0_dbl_kind)
    sib%diag%ecmass   = sib%diag%eci*facks * hltmi
    sib%diag%egmass   = sib%diag%egi*facks * hltmi

!!-----------------------------------------------------------------
!
!   Calculation of phase change within snow and soil layers.
!   This code based on routine CLM_MELTFREEZE from the Common
!   Land Model (CLM)  (Dai et al, submitted)
!   CLM web info:  http://clm.gsfc.nasa.gov
!------------------------------------------------------------------

    !...initialize some values
    sib%diag%snowmelt = 0.0

    do j = sib%prog%nsl+1,nsoil
        sib_loc%imelt(j)     = 0
        flux_imb(j)  = 0.0
        water_imb(j) = 0.0
        wice0(j)     = sib%prog%www_ice(j)
        wliq0(j)     = sib%prog%www_liq(j)
        wmass(j)     = sib%prog%www_ice(j) + sib%prog%www_liq(j)
    enddo

    !...calculate some fluxes. 

    do j = sib%prog%nsl+1,nsoil-1

        !itb...remember: temperature increment was added in addinc.

        !itb...sflux is the 'previous timestep' flux value
        sflux(j) = sib%param%slamda(j)*(sib_loc%td_old(j+1)-sib_loc%td_old(j))

        !itb...sflx1 is the 'n+1 timestep' flux value
        sflx1(j) = sib%param%slamda(j)*(sib%prog%td(j+1) -sib%prog%td(j))

    enddo


    sflux(nsoil) = 0.0
    sflx1(nsoil) = 0.0

    !...melting check
    do j = sib%prog%nsl+1,nsoil

        !...if ice exists in warm conditions, melt it
        if(sib%prog%www_ice(j) > 0.0 .and. sib%prog%td(j) > tice) then
            sib_loc%imelt(j)  = 1
            tinc(j) = sib%prog%td(j) - tice
            sib%prog%td(j) = tice
        endif

        !...if water exists in cold conditions, freeze it
        if(sib%prog%www_liq(j) > 0.0 .and. sib%prog%td(j) < tice) then
            sib_loc%imelt(j)  = 2
            tinc(j) = sib%prog%td(j) - tice
            sib%prog%td(j) = tice
        endif

    enddo


    !...check for existence of snow, less depth than 0.01 m
    if(sib%prog%nsl == 0  .and. sib%prog%snow_mass > 0.0_dbl_kind) then
        if(sib%prog%td(1) > tice) then
            sib_loc%imelt(1) = 1
            tinc(1) = sib%prog%td(1) - tice
            sib%prog%td(1)    = tice
        endif
    endif

    !...change in storage
    do j = sib%prog%nsl + 1, nsoil
        htstor(j) = sib%param%shcap(j)*dti*   &
            (sib%prog%td(j) - sib_loc%td_old(j))
    enddo

    !itb...calculate energy surplus or loss if there was melting or freezing

    do j = sib%prog%nsl+1,nsoil
        if(sib_loc%imelt(j) > 0 ) then    ! did melting/freezing occur?
            flux_imb(j) = tinc(j) * sib%param%shcap(j) ! (J/m^2)
        endif
    enddo



    !...the CLM boys say these have been checked carefully-they are the 
    !...result of computed error in the tridiagonal matrix solver.

    do j=sib%prog%nsl+1,nsoil
        if(sib_loc%imelt(j) == 1 .and. flux_imb(j) < 0.0 ) then
            sib_loc%imelt(j)    = 0
            flux_imb(j) = 0.0
        endif

        if(sib_loc%imelt(j) == 2 .and. flux_imb(j) > 0.0 ) then
            sib_loc%imelt(j)    = 0
            flux_imb(j) = 0.0
        endif
    enddo 


    do j=sib%prog%nsl+1,nsoil

        !...melting or freezing occurring, some error present
        if(sib_loc%imelt(j) > 0 .and. abs(flux_imb(j))>0.0) then

            !itb...water_imb > 0 ==> melting ice
            !itb...water_imb < 0 ==> freezing water

            water_imb(j) = flux_imb(j) / lfus  !(kg of water)

            !...snow exists, but less than critical depth. CLM boys say this
            !...needs work.


            !...LEVEL 1 - top soil layer, small amt of snow present, 
            !...no snow layers
            if(j == 1) then  
                if(sib%prog%nsl == 0 .and. sib%prog%snow_depth >  &
                    0.01_dbl_kind .and. water_imb(j) > 0.0_dbl_kind) then 

                    temp = sib%prog%snow_mass
                    sib%prog%snow_mass = max(0.0_dbl_kind,temp-water_imb(j))
                    propor = sib%prog%snow_mass/temp
                    sib%prog%snow_depth = sib%prog%snow_depth * propor

                    heatr = flux_imb(j)*dti -       &
                        lfus*(temp-sib%prog%snow_mass)*dti

                    if(heatr > 0.0) then
                        water_imb(j) = heatr*dtt/lfus 
                        flux_imb(j) = heatr
                    else
                        water_imb(j) = 0.0
                        flux_imb(j)  = 0.0
                    endif

                    sib%diag%snowmelt = max(0.0_dbl_kind,  &
                        (temp-sib%prog%snow_mass))*dti
                    le_phase = lfus*sib%diag%snowmelt  ! w/m^2
                endif  
            endif ! j==1

            heatr = 0.0
            if(water_imb(j) > 0.0) then
                sib%prog%www_ice(j) = max(0.0_dbl_kind,wice0(j)-water_imb(j))
                heatr = flux_imb(j)*dti  &
                    - lfus*(wice0(j)-sib%prog%www_ice(j))*dti
            elseif(water_imb(j) < 0.0_dbl_kind) then
                sib%prog%www_ice(j) = min(wmass(j),wice0(j)-water_imb(j))
                heatr = flux_imb(j)*dti  &
                    - lfus*(wice0(j)-sib%prog%www_ice(j))*dti
            endif

            sib%prog%www_liq(j) = max(0.0_dbl_kind,wmass(j)-sib%prog%www_ice(j))

            if(abs(heatr) > 0.0) then
                sib%prog%td(j) = sib%prog%td(j) + dtt*heatr/sib%param%shcap(j)

                if(sib%prog%www_liq(j)*sib%prog%www_ice(j) > 0.0 )  &
                    sib%prog%td(j) = tice
            endif 

            le_phase = le_phase + lfus*(wice0(j) - sib%prog%www_ice(j))*dti

            if(sib_loc%imelt(j) == 1 .and. j < 1 ) then
                sib%diag%snowmelt = sib%diag%snowmelt +  &
                    max(0.0_dbl_kind, (wice0(j) - sib%prog%www_ice(j)))*dti
            endif



        endif   ! imelt \= 0 and flux_imb \= 0
    enddo 


    !...convert to (joules) - what is this needed for?
    !        energy_of_snowmelt = sib%diag%snowmelt * lfus

!-----------------------------------------------------------------
!
!      END OF CLM_MELTFREEZE CODE
!
!------------------------------------------------------------------

!itb...calculate change in energy of CAS

   sib%diag%cas_e_storage = sib_loc%dta * sib%diag%cas_cap_heat * dti


end subroutine update

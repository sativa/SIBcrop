!-----------------------------------------------------------------

subroutine sib_main ( sib )

    use kinds
    use sibtype
    use sib_const_module
    use physical_parameters

    implicit none

!---PARAMETERS----------------------------------------------------

    type(sib_t), intent(inout) :: sib
    type(sib_local_vars) :: sib_loc


!-----------------------------------------------------------------

!     REFERENCES: Sato, N., P. J. Sellers, D. A. Randall, E. K. Schneider, 
!          J. Shukla, J. L Kinter III, Y-T, Hou, and Albertazzi (1989) 
!          "Effects of implementing the simple biosphere model in a general
!          circulation model. J. Atmos. Sci., 46, 2767-2782.

!                 Sellers, P. J., D. A. Randall, C. J. Collatz, J. A. Berry,
!          C. B. Field, D. A. Dazlich, C. Zhang, G. Collelo (1996) A revised 
!          land-surface parameterization (SiB2) for atmospheric GCMs. Part 1:
!          Model formulation.


!     SUBROUTINES CALLED: 


!     FUNCTIONS CALLED:


!------------------------------------------------------------------
!
!   local variables
!
!------------------------------------------------------------------


    real(kind=dbl_kind) :: zzwind    ! adjusted wind measurement height (m) 
    real(kind=dbl_kind) :: zztemp    ! adjusted temp measurement height (m)


    !...water/energy balance variables...
    real(kind=dbl_kind), dimension(sib%prog%nsl+1:nsoil) :: &
        wwwliq_old   ! beginning-of-timetep soil/snow liquid
                                    !    (kg/m^2)
    real(kind=dbl_kind), dimension(sib%prog%nsl+1:nsoil) :: &
        wwwice_old   ! end-of-timestep soil/snow ice
                                    !    (kg/m^2)
    real(kind=dbl_kind)    :: capac_old(2) ! beginning-of-timestep
                                          ! interception storage (kg/m^2)
    integer(kind=int_kind) :: nsl_old  
                                    ! beginning-of-timestep # of snow layers


    real(kind=dbl_kind) :: tsum   ! temp variable
    real(kind=dbl_kind) :: cas_q  ! beginning-of-timestep CAS moisture
    real(kind=dbl_kind) :: tliq
    real(kind=dbl_kind) :: icet,sbeg,send

    real(kind=dbl_kind),dimension(1) :: ppl, ttl,tess
                                               ! temp vars for call to eau_sat
                                               ! dimension(1) to keep SGI compiler happy

    integer(kind=int_kind) :: i, j, k, n, ksoil, l


!------------------------------------------------------------------
!
!   end local variables
!
!------------------------------------------------------------------

    !...first guesses for ta and ea

    !...load previous timestep soil temps
    !...also load soil ice/water for water balance check
    !...at end of timestep

    nsl_old = sib%prog%nsl

    do i = sib%prog%nsl+1,nsoil
        sib_loc%td_old(i) = sib%prog%td(i)
        wwwliq_old(i) = sib%prog%www_liq(i)
        wwwice_old(i) = sib%prog%www_ice(i)
    enddo

    do i = sib%prog%nsl+1, 0
        sib_loc%frac_iceold(i) = sib%prog%www_ice(i)     &
            /(sib%prog%www_ice(i) + sib%prog%www_liq(i))
    enddo

    capac_old(1) = sib%prog%capac(1)
    capac_old(2) = sib%prog%capac(2)

    sib%prog%tha = sib%prog%ta  / sib%prog%bps(1)
    sib%prog%ea  = sib%prog%sha * sib%prog%ps / (0.622 + sib%prog%sha)
    sib%prog%em  = sib%prog%sh  * sib%prog%ps / (0.622 + sib%prog%sh)

    !...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX
    !
    !...temporarily hardwiring the carbon isotopic ratios of the mixed layer
    !...and respireation into SiBDRIVE
    sib%prog%d13cm    = -7.8              ! del13C of mixed layer
    !sib%param%d13cresp   = -28.0             ! del13C of respiration
    !
    !...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX

    sib%diag%cuprt = sib%prog%cupr * 0.001
    sib%diag%lsprt = sib%prog%lspr * 0.001  ! converting units to m/sec
    sib%diag%roff       = 0.0
    sib%diag%roffo      = 0.0
    sib%prog%pco2ap_old = sib%prog%pco2ap
    sib%prog%cas_old = sib%prog%cas

    !...calculate albedo/reflectance/transmissivity
    call rada2(sib,sib_loc)

    !...distribute incident radiation between canopy and surface
    call rnload(sib)

    sib%param%zpd_adj   = sib%param%z2 - ( sib%param%z2-sib%param%zp_disp )   &
                                           * sib%diag%canex

    sib%param%z0        = sib%param%z0d/( sib%param%z2-sib%param%zp_disp )    &
                                        * ( sib%param%z2-sib%param%zpd_adj )

    !print*,'SiB, zpd:',sib%param%zp_disp,sib%param%zpd_adj
    !print*,'SiB, z0:',sib%param%z0d,sib%param%z0
    !print*,'SiB,z2:',sib%param%z2,sib%param%z1

    sib%param%rbc       = sib%param%cc1/sib%diag%canex
    sib%param%rdc       = sib%param%cc2*sib%diag%canex


!     initialize energy and water budgets
!     Initialize heat capacities, soil properties
    call begtem(sib,sib_loc)

!   after begtem, get beginning-of-timestep CAS water in kg/m^2
    cas_q = sib%prog%sha * sib%prog%ros * sib%diag%cas_cap_co2

!     CALCULATE RADT USING RADIATION FROM PHYSICS AND CURRENT
!     LOSSES FROM CANOPY AND GROUND!

    call netrad(sib,sib_loc)

!     GET RESISTANCES FOR SIB, update stomatal resistance 
    call vntlat(sib,sib_loc)

!   this call for ustar, cu for oceanic value of z0 
!   this is only used for near-coastal gridpoints...
!    zzwind = sib%param%z2 - sib%param%zpd_adj + zwind
!    zztemp = sib%param%z2 - sib%param%zpd_adj + ztemp

!    call vmfcalzo(sib,zzwind,zztemp)

    !itb...calculate partial derivatives of the various heat fluxes
    !itb...with respect to ground/canopy/snow temp, as well as
    !itb...some other derivatives.

    call dellwf(sib,sib_loc)

    call delef(sib,sib_loc)

    call delhf(sib,sib_loc)

!   solve matrix of prognostic variables
    call sibslv(sib,sib_loc)

!    update prognostic variables, get total latent and sensible fluxes
    call addinc(sib,sib_loc)

    !itb...call eau_sat instead of vnqsat...
    ppl(1) = sib%prog%ps*100.0
    ttl(1) = sib%prog%ta

    call ess_eau(1,ppl,ttl,tess)

    sib%diag%eastar = tess(1)/100.0

    !...CAS relative humidity
    sib%diag%rha = sib%prog%ea / sib%diag%eastar

    !...update canopy and ground surface water stores; check that fluxes 
    !...are correct, no evaporating more water than is available.
    call update(sib,sib_loc)

    !...update water storage/interception on canopy and
    !...calculate throughfall

    call hydro_canopy(sib,sib_loc)

    !itb...precip in mm/sec after hydro_canopy: convert to mm/hour
    sib%diag%p0 = sib%diag%p0 * 3600.0

    !...update precipitation onto snowcover
    call hydro_snow(sib)

    !...update infiltration and soil water
    call hydro_soil(sib)

!itb...the next three routines need only be called if there is snow
!itb...snow layers are initialized in hydro_canopy.

!    if(sib%prog%nsl < 0) then              !snow layers present
    !...compact snow levels 

       call compact_snow(sib,sib_loc)

    !...combine snow levels where necessary
       call combine_snow(sib)

    !...subdivide snow levels where necessary
       call subdivide_snow(sib)

!    endif                                  !snow layers present

    !     for perpetual conditions, do not update soil moisture
!         if(.not.fixday) then
!            www(1) = wwwtem(1)
!            www(2) = wwwtem(2)
!            www(3) = wwwtem(3)
!         endif

!         do i = 1,len
!            xgpp(i) = assim(i) * 1.e6
!         enddo
!     Calculate the surface flux of CO2 due to SiB (tracer T18)

!     The net flux to the atmosphere is given by the release of CO2
!     by soil respiration minus the uptake of CO2 by photosynthesis

!     ASSIMN is the net assimilation of CO2 by the plants (from SiB2)
!     respFactor*soilScale is the rate of release of CO2 by the soil

!     soilScale is a diagnostic of the instantaneous rate of 
!        soil respiration (derived by Jim Collatz, similar to TEM)
!     respFactor is the annual total accumulation of carbon in the
!        previous year at each grid cell (annual total ASSIMN) 
!        divided by the annual total of soilScale at the same grid pt.

!     Surface flux of CO2 used to be merely Assimn-Respg. With the
!     prognostic CAS, the calculation becomes
!
!     co2flux =  (CO2A - CO2M)/ra
!
!     with a temperature correction thrown in. This calculation is 
!     performed in phosib.


    !...check that energy and water balance is maintained
    call balan(sib,nsl_old,wwwliq_old,wwwice_old,capac_old,cas_q)

end subroutine sib_main

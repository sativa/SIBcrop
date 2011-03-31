!
!================SUBROUTINE BEGTEM=======================================
!
subroutine begtem(sib,sib_loc)


    use kinds
    use sibtype

    use physical_parameters, only: &
        cp => spec_heat_cp, &
        grav,               &
        hltm,               &
        tice

    use sib_const_module, only: &
        nsoil,  &
        denh2o, &
        denice, &
        tkwat,  &
        tkice,  &
        tkair,  &
        cww,    & ! water heat capacity
        clai,   &
        cv,     &
        cpice,  &
        cpliq

    implicit none

    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    type(sib_local_vars)     ,intent(inout) :: sib_loc
    ! variables local to SiB

    !----------------------------------------------------------------------  


    !      References

    !      Sellers, P.J., M.D. Heiser, F.G. Hall, S.J. Goetz, D.E. Strebel,
    !                     S.B. Verma, R.L. Desjardins, P.M. Schuepp,
    !                     J.I. MacPherson, 1995: Effects of Spatial Variability
    !                     in Topography, Vegetation Cover and Soil Moisture on 
    !                     Area-aAveraged Surface Fluxes: A Case Study Using the
    !                     FIFE 1989 Data. JGR, 100(D12), 25607-25629.

    !      Sellers, P.J. and Mintz, Y., Y.C. Sud, A. Dalcher, 1986: A Simple 
    !                     Biospher Model (SiB) for use Within General 
    !                     Circulation Models. JAS, 43(6),505-531.

    !========================================================================
    !
    !     Calculation of flux potentials and constants prior to heat 
    !         flux calculations.  
    !
    !======================================================================== 


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       CZC            CANOPY HEAT CAPACITY (J M-2 K-1)
    !       PSY            PSYCHROMETRIC CONSTANT (hPa K-1)
    !       TKSOIL         SOIL/SNOW THERMAL CONDUCTIVITY (W M-1 K-1)
    !       SHCAP          SOIL/SNOW HEAT CAPACITY (J M-2 K-1)
    !       ETC            E(star) TC-vapor pressure at leaf sfc
    !       GETC           dE(star)/dTC
    !       ETG            E(star) TG-vapor pressure at ground sfc
    !       GETG           dE(star)/dTG
    !       ETS            E(star) TS-vapor pressure at snow sfc
    !       GETS           dE(star)/dTS
    !       WC             CANOPY WETNESS FRACTION
    !       WG             GROUND WETNESS FRACTION
    !       RSTFAC(2)      SOIL MOISTURE STRESS FACTOR 
    !       RSOIL          SOIL SURFACE RESISTANCE (S M-1)
    !       HR             SOIL SURFACE RELATIVE HUMIDITY


    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !......LOCAL VARIABLES
    integer :: i,j
    real(kind=dbl_kind) :: & 
        one,   & ! must be for numerical purposes...
        satw,  & ! fraction of saturation, soil level
        fl,    & ! fraction liquid in soil layer
        dke,   & ! Kersten number
        bw,    & ! partial density of water (ice+liquid)
        fac,   & ! fraction of saturation of top soil layer
        psit,  & ! moisture potential of top soil layer
        argg,  & ! RH of air at soil surface
        phroot   ! normalized value of soil moisture potential
    ! over all soil levels


    real(kind=dbl_kind), dimension(nsoil) :: &
        dksat,        & ! thermal conductivity, saturated soil (W/m/K)
        phroot_layer, & ! layer values of phi-soil moisture potential
        vwc             ! volumetric water content (theta) (-)
   

    real(kind=dbl_kind), dimension(-nsnow+1:nsoil) :: &
        thk ! thermal conductivity of layer (W/m/K)

    real(kind=dbl_kind),dimension(1) :: ppl,ttl,esst,dtesst
    ! holder variables to make SGI 
    ! compiler happy for calls to
    ! eau_sat


    one = 1.0          

    !----------------------------------------------------------------------
    !     CALCULATION OF CANOPY HEAT CAPACITIES.
    !----------------------------------------------------------------------
    sib%param%czc = sib%param%zlt*clai+(0.5*sib%prog%snow_veg/denice + &
        sib%prog%capac(1)/denh2o)*cww
    sib%diag%psy = cp / hltm * sib%prog%ps / .622    


    !-------------------------------------------------------------------
    !   calculation of soil/snow thermal conductivity 
    !   and heat capacity
    !   Based on Farouki (1981), Jordan (1991), de Vires (1963),
    !   and on CLM coding
    !-------------------------------------------------------------------

    do j=1,nsoil

        !...fractional expression of saturation
        satw = ((sib%prog%www_liq(j)/denh2o) + (sib%prog%www_ice(j)/denice))/ &
            (sib%prog%dz(j) * sib%param%poros)


        satw = min(1.0_dbl_kind,satw)

        if(satw > 1.0E-6) then                 ! water present in soil
            fl = sib%prog%www_liq(j) / (sib%prog%www_liq(j) + &
                sib%prog%www_ice(j)) ! frac liq

            if(sib%prog%td(j) >= tice) then ! unfrozen soil
                ! kersten number
                dke = max(0.0_dbl_kind, log10(satw) + 1.0_dbl_kind)  
                dksat(j) = sib%param%tksatu(j)
            else ! frozen soil
                dke = satw
                dksat(j) = sib%param%tkmg(j) * 0.249**(fl*sib%param%poros) &
                    *2.29**sib%param%poros
            endif

            thk(j) = dke*dksat(j) + (1.0-dke)*sib%param%tkdry(j)

        else ! soil is very dry
            thk(j) = sib%param%tkdry(j)
        endif

    enddo  !nsoil loop

    !...thermal conductivity of snow
    if (sib%prog%nsl < 0) then
        do j= sib%prog%nsl+1,0
            bw = (sib%prog%www_liq(j) + sib%prog%www_ice(j))/sib%prog%dz(j)
            thk(j) = tkair + (7.75E-5*bw + 1.105E-6*bw*bw)*(tkice-tkair) 
        enddo
    endif

    !...thermal conductivity at the layer interface
    do j = sib%prog%nsl+1,nsoil-1

        sib%param%tksoil(j) = thk(j)*thk(j+1)* &
            (sib%prog%node_z(j+1)-sib%prog%node_z(j))  &
            /(thk(j)*(sib%prog%node_z(j+1)-sib%prog%layer_z(j)) + &
            thk(j+1)*(sib%prog%layer_z(j)-sib%prog%node_z(j)))

    enddo


    !itb...THIS IS A PATCH
    !      sib%param%tksoil(nsoil) = 0.0
    sib%param%tksoil(nsoil) = sib%param%tksoil(nsoil-1)

    !......heat capacity, soil
    do j=1,nsoil
        sib%param%shcap(j) = sib%param%csolid(j)*(1.0-sib%param%poros)*sib%prog%dz(j) &
            + sib%prog%www_ice(j)*cpice + sib%prog%www_liq(j)*cpliq
    enddo


    !......heat capacity, snow
    if(sib%prog%nsl < 0) then
        do j=sib%prog%nsl+1,0
            sib%param%shcap(j) = cpliq*sib%prog%www_liq(j) + &
                cpice*sib%prog%www_ice(j)
        enddo
    endif

    !...slamda is the lambda(z_sub h,j)/(zj+1 - zj) term from the CLM
    !...document. slamda*(Tj+1 - Tj) is the heat flux from layer j to 
    !...layer j+1. This is a different numerical scheme than Bonan.

    do j=sib%prog%nsl+1,nsoil-1
        sib%param%slamda(j) = sib%param%tksoil(j) / (sib%prog%node_z(j+1) - &
            sib%prog%node_z(j))
    enddo
    

    !...THIS IS A PATCH
    sib%param%slamda(nsoil) = sib%param%slamda(nsoil-1)

    !
    !----------------------------------------------------------------------
    !      Calculation of ground surface temperature and wetness fractions
    !        
    !----------------------------------------------------------------------
    !

    !...get saturation vapor pressure ('e-star') values for canopy,
    !...soil, and snow

    sib%prog%tg = sib%prog%td(sib%prog%nsl+1)

    ppl(1) = sib%prog%ps*100.0
    ttl(1) = sib%prog%tc

    call dtess_eau(1,ppl,ttl,esst,dtesst)

    sib_loc%etc  = esst(1)
    sib_loc%getc = dtesst(1)

    ttl(1) = sib%prog%tg

    call dtess_eau(1,ppl,ttl,esst,dtesst)

    sib_loc%etg = esst(1) 
    sib_loc%getg = dtesst(1)

    !...vapor pressures come out of eau_sat in Pa, so convert to mb
    sib_loc%etc  = sib_loc%etc/100.0
    sib_loc%getc = sib_loc%getc/100.0
    sib_loc%etg  = sib_loc%etg/100.0
    sib_loc%getg = sib_loc%getg/100.0


    if(sib%prog%nsl < 0 ) then
        sib%diag%tsnow = sib%prog%td(sib%prog%nsl+1)

        ttl(1) = sib%diag%tsnow

        call dtess_eau(1,ppl,ttl,esst,dtesst)

        sib_loc%ets = esst(1)/100.0
        sib_loc%gets = dtesst(1)/100.0

    else
        sib%diag%tsnow = sib%prog%tg
        sib_loc%ets   =  sib_loc%etg
        sib_loc%gets  =  sib_loc%getg
    endif

    !...canopy and ground fractional wetness...

    sib%diag%wc = MIN( one,( sib%prog%capac(1)/denh2o + &
        sib%prog%snow_veg)/sib%param%satcap(1) )

    sib%diag%wg = MAX( 0.*one, &
        (sib%prog%capac(2)/(sib%param%satcap(2)*denh2o)) )*0.25

    !-----------------------------------------------------------------------
    !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
    !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
    !     SOURCE FOR TRANSPIRATION.
    !
    !      PHROOT      (PSI-R) : EQUATION (48) , SE-86
    !      RSTFAC(2)  F(PSI-L) : MODIFICATION OF EQUATION (12), SE-89
    !
    !     CALCULATION OF WATER STRESS AND PLANT AVAILABLE WATER HAS BEEN
    !     CHANGED IN SiB3. 
    !-----------------------------------------------------------------------
    !

!bio...wssp is the 'water stress shape parameter'. reasonable values will be
!bio...between 0.1 and 1.0. This variable controls the shape of the water stress 
!bio...curve. A value of 1.0 yields a linear stress between wilt point and field 
!bio...capacity, and a small value gives a steep drop-off when WP is neard-not as
!bio...much stress at higher vwc amounts.

    sib%diag%wssp = 0.2
    
    sib%diag%paw_tot    = 0.0
    sib%diag%paw_max    = 0.0
    sib%diag%paw(:)     = 0.0

!bio...Assumption is that there will be no stress if volumetric water content
!bio...is above field capacity. Entire soil column is considered. Calculate
!bio...amount that soil water is below FC.

    do i=1,nsoil
        sib%prog%vol_ice(i) = min(sib%param%poros,sib%prog%www_ice(i)/   &
            (sib%prog%dz(i)*denice))
        sib%diag%eff_poros(i) = sib%param%poros - sib%prog%vol_ice(i)
        sib%prog%vol_liq(i) = min(sib%diag%eff_poros(i),    &
            sib%prog%www_liq(i)/(sib%prog%dz(i)*denh2o))
	    
        sib%diag%paw(i) = sib%prog%vol_liq(i) - sib%param%vwcmin
        sib%diag%paw(i) = max(sib%diag%paw(i),0.0_dbl_kind)
        sib%diag%paw_tot = sib%diag%paw_tot +      &
                                sib%diag%paw(i) * sib%prog%dz(i)
    	sib%diag%paw_max = sib%diag%paw_max +           &
                         ((sib%param%fieldcap - sib%param%vwcmin) * sib%prog%dz(i))
	enddo
	
!bio...calculate stress factor. If total column soil moisture is at or below wilt 
!bio...point, stress=0 (total stress). No stress if column-mean vwc is above FC.

	 sib%diag%pawfrac = sib%diag%paw_tot/sib%diag%paw_max
     sib%diag%pawfrac = MAX(0.0_dbl_kind, sib%diag%pawfrac)
     sib%diag%pawfrac = MIN(sib%diag%pawfrac, 1.0_dbl_kind)

	 sib%diag%rstfac(2) = ((1+sib%diag%wssp) * sib%diag%pawfrac)/    &
	                      (sib%diag%wssp + sib%diag%pawfrac)
	
    !itb...maintain rstfac2 at or above 0.1
    sib%diag%rstfac(2) = MAX(sib%diag%rstfac(2), 0.1_dbl_kind)

    !EL..under irrigated conditions, changed the rstfac(2) to 1.0 
    !to avoid the moisture stress.
    !kdcorbin, 03/11 - added test for crops
    !if (sib%param%biome .ge. 20) then
    !    sib%diag%rstfac(2) = MAX(sib%diag%rstfac(2), 1.0_dbl_kind)
    !endif

    !----------------------------------------------------------------------
    !
    !      RSOIL FUNCTION FROM FIT TO FIFE-87 DATA.  Soil surface layer
    !         relative humidity.
    !
    !      RSOIL      (RSOIL) : EQUATION (3), Sellers et al, JGR, 100(D12)
    !                                          25607-25629
    !      HR         (Fh)    : EQUATION (66) , SE-86
    !----------------------------------------------------------------------
    !
    fac = MIN( (sib%prog%www_liq(1)+sib%prog%www_ice(1))/ &
        (sib%param%poros*sib%prog%dz(1)*1000.0), one )
    fac = MAX( fac, 0.02*one  )

    sib%diag%rsoil =  exp(8.206 - 4.255 * fac)    

    psit = sib%param%phsat * fac ** (- sib%param%bee )
    argg = max(-10.*one,(psit*grav/ (461.5*sib%prog%tg) ))
    sib%diag%hr = exp(argg)

    !...storage inertia terms...
    sib%diag%cas_cap_heat = sib%prog%ros * cp * max(4.0_dbl_kind,sib%param%z2)
    sib%diag%cas_cap_vap  = sib%prog%ros * cv * max(4.0_dbl_kind,sib%param%z2) &
        / sib%diag%psy
    sib%diag%cas_cap_co2  = max(4.0_dbl_kind,sib%param%z2)          
    ! this goes 
    ! out to phosib

end subroutine begtem

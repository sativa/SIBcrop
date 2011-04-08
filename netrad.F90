!==================SUBROUTINE NETRAD====================================
subroutine netrad (sib,sib_loc)

    use kinds
    use sibtype
    use physical_parameters, only: &
        stefan, &
        tice

    implicit none

    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    type(sib_local_vars)     ,intent(inout) :: sib_loc
    ! variables local to SiB

    !----------------------------------------------------------------------  
    !
    !=======================================================================
    !
    !                                                                       
    !        CALCULATE RADT USING RADIATION FROM PHYSICS AND CURRENT        
    !        LOSSES FROM CANOPY AND GROUND                                  
    !
    !
    !=======================================================================
    !
    !... bands in sib: 0.2 to 0.7 microns are VIS, then 
    !...               0.7 to 4.0 is NIR, above 4.0 it is thermal

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       RADt (2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2) 
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !Bio...LOCAL VARIABLES
    real(kind=dbl_kind) :: &
        tc4,     & ! canopy temp **4
        tg4,     & ! ground temp **4
        ts4,     & ! snow temp **4
        radtbar, & ! bulk weighted net rad
        feedfac, & ! feedback factor
        zlwup      ! total thermal rad up from sfc (W/m^2)

    real(kind=dbl_kind),dimension(1) :: ppl,ttl,qsatst
    ! holder arrays to make SGI compiler happy
    ! for calls 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    data feedfac/1.0/

    if(sib%prog%nsl < 0) then
        sib%diag%tsnow = sib%prog%td(sib%prog%nsl+1)
    else
        sib%diag%tsnow = min(sib%prog%td(1),tice)
    endif

    tc4 = sib%prog%tc*sib%prog%tc*sib%prog%tc*sib%prog%tc
    tg4 = sib%prog%td(1)*sib%prog%td(1)*sib%prog%td(1)*sib%prog%td(1)
    ts4 = sib%diag%tsnow*sib%diag%tsnow*sib%diag%tsnow*sib%diag%tsnow

    !...effective ground cove3r for thermal radiation
    sib_loc%fac1 = sib%param%vcover * ( 1.-sib%diag%thermk )

    !...derivatives
    sib_loc%dtc4 = 4*stefan * sib%prog%tc**3
    sib_loc%dtg4 = 4*stefan * sib%prog%tg**3
    sib_loc%dts4 = 4*stefan * sib%diag%tsnow**3

    !...canopy leaves thermal radiation loss
    sib_loc%closs =  2. * sib_loc%fac1 * stefan * tc4
    sib_loc%closs =  sib_loc%closs - sib_loc%fac1 * stefan * &
        ( (1.-sib%diag%areas)*tg4+sib%diag%areas*ts4)

    !...ground thermal radiation loss 
    sib_loc%gloss =  stefan * tg4 - sib_loc%fac1 * stefan * tc4

    !...snow thermal radiation loss 
    sib_loc%sloss =  stefan * ts4 - sib_loc%fac1 * stefan * tc4

    !...canopy leaves net radiation
    sib%diag%radt(1) = sib%diag%radc3(1) - sib_loc%closs

    !...ground net radiation
    sib%diag%radt(2) = sib%diag%radc3(2) - sib_loc%gloss

    !...snow net radiation 
    sib%diag%radt(3) = sib%diag%radc3(2) - sib_loc%sloss

    !...bulk, weighted net radiation from combined ground and snow
    radtbar = sib%diag%areas*sib%diag%radt(3) + (1.-sib%diag%areas)* &
        sib%diag%radt(2)

    !...this is the exchange meant to help out exchanges between
    !...ground and snow
    sib%diag%radt(2) = radtbar + (1.+feedfac)*(sib%diag%radt(2)-radtbar)
    sib%diag%radt(3) = radtbar + (1.+feedfac)*(sib%diag%radt(3)-radtbar)

    !...total thermal radiation up from surface
    zlwup = sib_loc%fac1 * tc4 + &
        (1.-sib_loc%fac1) * (sib%diag%areas*ts4+(1.-sib%diag%areas)*tg4)
    !...effective (combined) skin temperature from surface thermal radiation
    sib%diag%tgeff =  zlwup ** 0.25

    !...potential skin temp 
    sib%diag%thgeff = sib%diag%tgeff / sib%prog%bps(1)

    !itb...use a call to eau_sat instead of vnqsat
    !...have to make SGI compiler happy with array arguments, not scalars...

    ppl(1) = sib%prog%ps*100.0
    ttl(1) = sib%diag%thgeff

    call qsat_eau(1,ppl,ttl,qsatst)

    sib%diag%shgeff = qsatst(1)


end subroutine netrad

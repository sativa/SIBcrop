
!===================SUBROUTINE DELLWF=====================================


subroutine dellwf(sib,sib_loc)

    use kinds
    use sibtype


    !========================================================================
    !
    !     Calculation of partial derivatives of canopy and ground radiative
    !        heat fluxes with respect to Tc, Tg
    !     Here we are doing only the long wave radiative loss, which is the
    !     only radiative quantity we are trying to bring to the next time step.
    !
    !======================================================================== 

    !------------------------------INPUT is coming from Netrad-------------
    !
    !       dtc4, dtg4, dts4, which are the derivatives of the LW loss
    !
    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       LCDTC          dLC/dTC 
    !       LCDTG          dLC/dTG
    !       LCDTS          dLC/dTS
    !       LGDTG          dLG/dTG
    !       LGDTC          dLG/dTC
    !       LSDTS          dLS/dTS
    !       LSDTC          dLS/dTC
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT none

    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    type(sib_local_vars)     ,intent(inout) :: sib_loc
    ! variables local to SiB

    !----------------------------------------------------------------------  
    !
    !   canopy leaves:
    sib_loc%lcdtc =   2 * sib_loc%dtc4 * sib_loc%fac1
    sib_loc%lcdtg =     - sib_loc%dtg4 * sib_loc%fac1 * (1.-sib%diag%areas)
    sib_loc%lcdts =     - sib_loc%dts4 * sib_loc%fac1 * (   sib%diag%areas)
    !
    !   ground:
    !
    sib_loc%lgdtg =   sib_loc%dtg4
    sib_loc%lgdtc = - sib_loc%dtc4 * sib_loc%fac1
    !
    !   snow:
    !
    sib_loc%lsdts =   sib_loc%dts4
    sib_loc%lsdtc = - sib_loc%dtc4 * sib_loc%fac1


end subroutine dellwf                                             

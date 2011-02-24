!==============SUBROUTINE RNLOAD========================
subroutine rnload(sib)

    use kinds
    use sibtype

    implicit none

    !---------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    !---------------------------------------------------------------  

    !
    !====================================================
    !
    !    calculation of absorption of radiation by surface.  Note that
    !       output from this calculation (radc3) only accounts for the 
    !       absorption of incident longwave and shortwave fluxes.  The
    !       total net radiation calculation is performed in subroutine
    !       netrad.
    !
    !====================================================
    !

    !+++++++++++++++++++++++++OUTPUT+++++++++++++++++++++
    !
    !       RADN(2,3)      INCIDENT RADIATION FLUXES (W M-2)
    !       RADC3(2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2) 
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++



    integer(kind=int_kind) :: i, iveg, iwave, irad
    real(kind=dbl_kind)    :: radn(2,2)

    !-------------------------------------------------------------
    !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
    !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
    !     SOURCE FOR TRANSPIRATION.
    !
    !      RADN        (F(IW,IMU,O)) : EQUATION (19-22) , SE-86
    !      RADC3       (FC,FGS)      : EQUATION (21,22) , SE-86
    !--------------------------------------------------------------


    sib%diag%radc3(1) = 0.
    sib%diag%radc3(2) = 0.
    radn(1,1) = sib%prog%radvbc
    radn(1,2) = sib%prog%radvdc
    radn(2,1) = sib%prog%radnbc
    radn(2,2) = sib%prog%radndc

    do iveg=1,2
        do iwave=1,2
            do irad=1,2
                sib%diag%radc3(iveg) = sib%diag%radc3(iveg) +  &
                    sib%diag%radfac(iveg,iwave,irad) *         &
                    radn(iwave,irad)
            enddo
        enddo
    enddo

    !...absorb downwelling radiation 

    sib%diag%radc3(1) = sib%diag%radc3(1) + sib%prog%dlwbot *   &
        sib%param%vcover * (1.- sib%diag%thermk)
    sib%diag%radc3(2) = sib%diag%radc3(2) + sib%prog%dlwbot *   &
        (1.-sib%param%vcover * (1.-sib%diag%thermk))

end subroutine rnload

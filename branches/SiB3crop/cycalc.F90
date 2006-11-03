
!=====================SUBROUTINE CYCALC=================================


subroutine cycalc( aparkk, vm, atheta, btheta, par, &
        gammas, respc, rrkk, omss, c3, c4, &
        pco2i, assimn, assim)

    use kinds

    implicit none

    !=======================================================================
    !
    !     CALCULATION EQUIVALENT TO STEPS IN FIGURE 4 OF SE-92A
    !     C4 CALCULATION BASED ON CO-92.
    ! 
    !=======================================================================


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !     
    !       ASSIM         
    !       ASSIMN
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       OMC            RUBISCO LIMITED ASSIMILATION (MOL M-2 S-1)
    !       OME            LIGHT LIMITED ASSIMILATION (MOL M-2 S-1)
    !       OMS            SINK LIMITED ASSIMILATION (MOL M-2 S-1)
    !       CO2S           CANOPY SURFACE CO2 CONCENTRATION (MOL MOL-1)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !Bio...INPUT VARIABLES
    real(kind=dbl_kind),intent(in) :: &
        aparkk, & !
        vm,     & !
        atheta, & !
        btheta, & !
        gammas, & !
        par,    & !
        respc,  & !
        rrkk,   & !
        omss,   & !
        c3,     & !
        c4,     & !
        pco2i     !

    !Bio...OUTPUT VARIABLES
    real(kind=dbl_kind),intent(out) :: &
        assimn, & !
        assim     !

    !Bio...LOCAL VARIABLES
    real(kind=dbl_kind) :: &
        ome,   & !
        omc,   & !
        omp,   & !
        oms,   & !
        sqrtin   !


    !-----------------------------------------------------------------------
    !     CALCULATE ASSIMILATION RATE
    !
    !      OMC         (OMEGA-C): EQUATION (11) , SE-92A
    !      OME         (OMEGA-E): EQUATION (12) , SE-92A
    !      OMS         (OMEGA-S): EQUATION (13) , SE-92A
    !      ASSIMN      (A-N)    : EQUATION (14,15), SE-92A
    !-----------------------------------------------------------------------

    omc = vm *(pco2i-gammas)/(pco2i + rrkk)*c3    &
        + vm * c4
    ome = par*(pco2i-gammas)/(pco2i+2.*gammas)*c3 &
        + par * c4
    sqrtin= MAX( 0.0_dbl_kind, ( (ome+omc)**2 - 4.*atheta*ome*omc ) )
    omp  = ( ( ome+omc ) - SQRT( sqrtin ) ) / ( 2.*atheta )
    oms  = omss * c3 + omss*pco2i * c4
    sqrtin= MAX( 0.0_dbl_kind, ( (omp+oms)**2 - 4.*btheta*omp*oms ) ) 
    assim = ( ( oms+omp ) - SQRT( sqrtin ) ) / ( 2.*btheta )
    assimn= ( assim - respc) * aparkk

    return
end subroutine cycalc   

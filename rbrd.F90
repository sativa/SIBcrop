!==================SUBROUTINE RBRD=======================================
subroutine rbrd(sib,u2)

use kinds
use sibtype
use physical_parameters, only: grav    

implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  

!      Reference

!      Sellers, P.J. and Mintz, Y., Y.C. Sud, A. Dalcher, 1986: A Simple 
!                     Biospher Model (SiB) for use Within General 
!                     Circulation Models. JAS, 43(6),505-531.

!========================================================================
!
!      CALCULATION OF RB AND RD AS FUNCTIONS OF U2 AND TEMPERATURES
!
!======================================================================== 



!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       RB (GRB)       CANOPY TO CAS AERODYNAMIC RESISTANCE (Sec M-1)
!       RD (GRD)       GROUND TO CAS AERODYNAMIC RESISTANCE (Sec M-1)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!...input variable

real(kind=dbl_kind) :: u2

!Bio...LOCAL VARIABLES
real(kind=dbl_kind) :: temdif    
! vegetation-CAS temperature difference (K)

real(kind=dbl_kind) :: fac    ! factor for vegetation-CAS resistance calc
real(kind=dbl_kind) :: fih    ! factor for ground-CAS resistance calc
real(kind=dbl_kind) :: tgs    ! composite soil/snow sfc temperature

    tgs = (1.0_dbl_kind - sib%diag%areas) * sib%prog%td(1) +  &
        sib%diag%areas * sib%prog%td(sib%prog%nsl+1)

    !-----------------------------------------------------------------------
    !      RB       (RB)       : EQUATION (A9), SE-86
    !-----------------------------------------------------------------------

    temdif  = MAX( 0.1_dbl_kind,  sib%prog%tc - sib%prog%ta)
    fac     = sib%param%zlt / 890.* (temdif * 20.0)**0.25

    sib%diag%rb  = 1.0 / (SQRT(u2) / sib%param%rbc+fac)

    !-----------------------------------------------------------------------
    !      RD       (RD)       : EQUATION (A15), SE-86
    !-----------------------------------------------------------------------

    !itb...for rd, we use the composite soil/snow temperature
    tgs = sib%diag%areas*sib%prog%td(sib%prog%nsl+1) +   &
        (1.0 - sib%diag%areas)*sib%prog%td(1)
    temdif = MAX( 0.1_dbl_kind, tgs-sib%prog%ta )
!print*,sib%param%z2

    fih = SQRT( 1.+9.* grav * temdif * sib%param%z2 / (tgs*u2*u2) )
    
    sib%diag%rd  = sib%param%rdc / (u2 * fih) 

end subroutine rbrd

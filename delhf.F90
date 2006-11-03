 
!===================SUBROUTINE DELHF25=====================================
 
subroutine delhf(sib,sib_loc)

use kinds
use sibtype

use physical_parameters, only :                                 &
    cp => spec_heat_cp

use sib_const_module, only    :                                 &
    dtt

!========================================================================
!
!     Calculation of partial derivatives of canopy and ground sensible
!        heat fluxes with respect to Tc, Tg, and Theta-m.
!     Calculation of initial sensible heat fluxes.
!
!======================================================================== 


!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       HC             CANOPY SENSIBLE HEAT FLUX (J M-2)
!       HG             GROUND SENSIBLE HEAT FLUX (J M-2)
!       HS             SNOW   SENSIBLE HEAT FLUX (J M-2)
!       HA             CAS    SENSIBLE HEAT FLUX (J M-2)
!       HCDTC          dHC/dTC 
!       HCDTA          dHC/dTA
!       HGDTG          dHG/dTG
!       HGDTA          dHG/dTA
!       HSDTS          dHS/dTS
!       HSDTA          dHS/dTA
!       HADTA          dHA/dTA
!       HADTH          dHA/dTH
!       AAC            dH/dTC
!       AAG            dH/dTG
!       AAM            dH/dTH
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

type(sib_local_vars)     ,intent(inout) :: sib_loc
                           ! variables local to SiB

!----------------------------------------------------------------------  


!     local variables
!     REALX8 D1, d1i,    
real(kind=dbl_kind) :: rai    ! 1/ra
real(kind=dbl_kind) :: rbi    ! 1/rb
real(kind=dbl_kind) :: rdi    ! 1/rd
integer i                                      

!-----------------------------------------------------------------------
!                                                                       
!     FLUXES EXPRESSED IN JOULES M-2, although in SIBSLV WE THEN WANT W/m2
!    WHY ????
!                                                                       
!    if we were to keep things simple, there is no need to separate
!    HG and HS, but it helps the derivatives keep clean.
!
!      HC          (HC)    : EQUATION (63) , SE-86
!      HG          (HG)    : EQUATION (65) , SE-86
!      HS          (HS)    : EQUATION (65) , SE-86
!      HA          (HA)    : EQUATION ???
!-----------------------------------------------------------------------

    rai = 1.0 / sib%diag%ra
    rbi = 1.0 / sib%diag%rb  
    rdi = 1.0 / sib%diag%rd              

    !    these are the current time step fluxes in J/m2
    !    can we change this to W/m2 ???
    sib%diag%hc   = cp * sib%prog%ros * (sib%prog%tc    - sib%prog%ta)   &
        * rbi * dtt

    if(sib%prog%nsl == 0 ) then   !no snow case
        sib%diag%hg   = cp * sib%prog%ros *               &
            (sib%prog%td(1) - sib%prog%ta) * rdi * dtt 
        sib%diag%hs   = 0.0
    else                          ! snow case
        sib%diag%hg   = cp * sib%prog%ros *               &
            (sib%prog%td(1) - sib%prog%ta) * rdi * dtt 

        !         sib%diag%hg   = sib%diag%hg * (1.0 - sib%diag%areas)
        sib%diag%hg = 0.0

        sib%diag%hs   = cp * sib%prog%ros *               &
            (sib%diag%tsnow - sib%prog%ta) * rdi * dtt 

        !         sib%diag%hs   = sib%diag%hs * sib%diag%areas
    endif


    sib%diag%fss  = cp * sib%prog%ros * (sib%prog%ta    - sib%prog%tm)    &
        * rai * dtt

    !    now we do the partial derivatives
    !    these are done assuming the fluxes in W/m2        

    !    for canopy leaves sensible heat flux: W/(m2 * K)
    ! 
    sib_loc%hcdtc =   cp * sib%prog%ros * rbi
    sib_loc%hcdta = - sib_loc%hcdtc
    !
    !    for ground and snow sensible heat fluxes: W/(m2 * K)
    !
    sib_loc%hgdtg =   cp * sib%prog%ros * rdi  
    sib_loc%hsdts =   sib_loc%hgdtg  
    sib_loc%hgdta = - sib_loc%hgdtg  
    sib_loc%hsdta = - sib_loc%hgdtg  
    !
    !    for the canopy air space (CAS) sensible heat flux: W/(m2 * K)
    !
    sib_loc%hadta =   cp * sib%prog%ros * rai
    sib_loc%hadth = - sib_loc%hadta/sib%prog%bps(1)

    !    ATTENTION !!!! DANGER !!!!! THIS WILL NOT WORK WITHOUT sibdrv = true
    !    for mixed layer (ref temp if not sibdrv): YET TO BE DONE
    !itb...LOOK AT SATO ET AL...
    !        AAG(I) = rdi * d1i                                    
    !        AAC(I) = rbi * d1i                                    
    !        AAM(I) = rai * d1i * bps(i)                          


end subroutine delhf                                              

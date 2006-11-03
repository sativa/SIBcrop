
!-SUBROUTINE: addinc-------------------------------------------------

subroutine addinc ( sib, sib_loc )

use kinds
use sibtype
use physical_parameters, only: cp => spec_heat_cp

implicit none

!-Parameters-----------------------------------------------------
type(sib_t), intent(inout) :: sib
type(sib_local_vars),      intent(inout) :: sib_loc

!-Local Variables------------------------------------------------
integer(kind=int_kind) :: j

!----------------------------------------------------------------
!
! Add prognostic variable increments to prognostic variables and
! diagnose heat fluxes and mixing ratios.
!
!----------------------------------------------------------------
!
!                           OUTPUT
!
!   tc      Canopy Temperature (K)
!   tg      Ground Surface Temperature (K)
!   ta      CAS Temperature (K)
!   ea      CAS Vapro Preassure (hPa or mb)
!   td      Deep Siol Temperature (K)
!   fss     CAS to Mixed Layer Sensible Heast Flux (W m^-2)
!   fws     CAS to Mixed Layer Latent Heat Flux (W m^-2)
!   sh      Mixed Layer Mixing Ratio (kg/kg)
!   sha     CAS Mixing Ratio (kg/kg)
!
!----------------------------------------------------------------

    sib%prog%tc = sib%prog%tc + sib_loc%dtc
    sib%prog%tg = sib%prog%tg + sib_loc%dtg
    sib%prog%ta = sib%prog%ta + sib_loc%dta
    sib%prog%ea = sib%prog%ea + sib_loc%dea

    if ( sib%prog%tc < 0.0 .or. sib%prog%tg < 0.0 .or.  &
        sib%prog%ta < 0.0 .or. sib%prog%ea < 0.0 ) then

        print *, 'point',sib%stat%pt_num
        print *, 'BAD ta OR ea VALUE:'
        print *, 'tc: ', sib%prog%tc-sib_loc%dtc, sib%prog%tc
        print *, 'tg: ', sib%prog%tg-sib_loc%dtg, sib%prog%tg
        print *, 'ta: ', sib%prog%ta-sib_loc%dta, sib%prog%ta
        print *, 'ea: ', sib%prog%ea-sib_loc%dea, sib%prog%ea
        print *, ''
    endif

    ! Calucluate latent and sunsible fluxes between CAS and
    ! mixed (boundry) layer
    sib%diag%fss = sib%prog%ros * cp * (sib%prog%ta - sib%prog%tm) / &
        sib%diag%ra
    sib%diag%fws = (sib%prog%ea - sib%prog%em) / sib%diag%ra * cp *  &
        sib%prog%ros / sib%diag%psy
    ! Recalculate mixing ratios
    sib%prog%sh  = 0.622 / (sib%prog%ps / sib%prog%em - 1.)
    sib%prog%sha = 0.622 / (sib%prog%ps / sib%prog%ea - 1.)
    do j = sib%prog%nsl+1, nsoil
        sib%prog%td(j) = sib%prog%td(j) + sib_loc%dtd(j)
    enddo


end subroutine addinc


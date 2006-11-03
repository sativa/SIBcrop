
subroutine hydro_snow(sib)

!----------------------------------------------------------------------
!
!     Evaluate and update snow mass/snow water.
!                 based on code in CLM_HYDRO_SNOW.F90
!
!
!     the following is taken verbatim from the comments in CLM_HYDRO_SNOW
!
!     Description:
!       Evaluate the change of snow mass and the snow water onto soil. 
!     Water flow within snow is computed by an expicit and non-physical
!     capacity (a tentative value is used, i.e. equal to 0.33*porosity)
!     to percolate into the underlying layer.  Except for cases where 
!     the porosity of one of the two neighboring layers is less than 0.05,
!     zero flow is assumed. The water flow out of the bottom of the snow
!     pack will participate as the input of the soil water and runoff.
!
!     Revision History:
!      15 September 1999: Yonjiu Dai; original code
!      15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!      15 November  2000: Mariana Vertenstein
!      22 January   2002: Ian Baker - integration into SiB
!----------------------------------------------------------------------

use kinds
use sibtype

use physical_parameters, only: &
    tice, &
    hltm

use sib_const_module, only: &
    cww,    &
    denice, &
    denh2o, &
    ssi,    &
    snomel, &
    dtt,    &
    dti


implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  


!...LOCAL VARIABLES...
integer(kind=int_kind) :: j,i
real(kind=dbl_kind)    :: wgdif
real(kind=dbl_kind)    :: qin
real(kind=dbl_kind)    :: qout
real(kind=dbl_kind)    :: qout_snow
real(kind=dbl_kind)    :: hfus

!----------------------------------------------------------------------


    hfus = snomel/1000.0   ! heat of fusion, units J/kg


    !...porosity and partial volume - recalculated
    if(sib%prog%nsl < 0) then
        do i = sib%prog%nsl + 1, 0
            sib%prog%vol_ice(i) = min(sib%param%poros, &
                sib%prog%www_ice(i)/(sib%prog%dz(i)*denice))
            sib%diag%eff_poros(i) = 1.0 - sib%prog%vol_ice(i)
            sib%prog%vol_liq(i) = min(sib%diag%eff_poros(i), &
                sib%prog%www_liq(i)/(sib%prog%dz(i)*denh2o))
        enddo
    endif

    if (sib%prog%nsl == 0) then  !no snow case

        sib%diag%www_inflow =  sib%diag%pcpg_rain
        ! units kg/m^2/sec 

    else  !snow present

        j = sib%prog%nsl + 1


        !itb...wgdif is the top snow layer ice amount, plus deposition (dew/frost) 
        !itb...minus sublimation. Keeping it zero for now, as a PATCH
        wgdif = 0.0

        if(wgdif < 0.0) then
            sib%prog%www_ice(j) = 0.0
            sib%prog%www_liq(j) = sib%prog%www_liq(j) + wgdif
        endif

        !itb...add rain
        if(sib%prog%tm >= tice) sib%prog%www_liq(j) = sib%prog%www_liq(j) &
            + sib%diag%pcpg_rain*dtt

        !itb...for the moment, not adjusting top snow layer liquid by fluxes.
        !        sib%prog%www_liq(j) =  sib%prog%www_liq(j) + precip + dew - evap/sublimation


        !itb---this directly from the comments in CLM_HYDRO_SNOW---
        !
        !...Capillary forces within snow are usually two or more orders of 
        !...magnitude less than those of gravity. Only gravity terms are
        !...considered.  The general expression for water flow is "K" * ss**3",
        !...however, no effective parameterization for "K".  Thus, a very
        !...simple consideration (not physically based) is introduced: when
        !...the liquid water of layer exceeds the layer's holding capacity,
        !...the excess meltwater adds to the underlying neighbor water.
        !

        !itb...CLM sets 0.033 as irreducible water saturation for snow (ssi)
        qin = 0.0
        do j = sib%prog%nsl+1,0
            sib%prog%www_liq(j) = sib%prog%www_liq(j) + qin
            if(j <= -1 ) then
                !...no runoff over snow surface, just ponding on surface
                !itb...USE SATCAP(2) HERE?
                if(sib%diag%eff_poros(j)<0.05 .or. &
                    sib%diag%eff_poros(j+1)< 0.05) then
                    qout = 0.0
                else
                    qout = max(0.0_dbl_kind,(sib%prog%vol_liq(j)- &
                        ssi * sib%diag%eff_poros(j)) * sib%prog%dz(j))
                    qout = min(qout,(1.0-sib%prog%vol_ice(j+1)- &
                        sib%prog%vol_liq(j+1))*sib%prog%dz(j+1))
                endif

            else
                qout = max(0.0_dbl_kind,(sib%prog%vol_liq(j)-ssi* &
                    sib%diag%eff_poros(j))*sib%prog%dz(j))
            endif

            qout           = qout * 1000.0
            sib%prog%www_liq(j) = sib%prog%www_liq(j) - qout
            qin            = qout
        enddo

        !itb...liquid out of the bottom of the snow into the top soil layer...
        qout_snow = qout*dti
        sib%diag%www_inflow = qout_snow
    endif   ! snow present/not present condition

    !itb...put sib%diag%www_inflow into sib%prog%capac(2) - ground interception
    !itb...before it is infiltrated.

    sib%prog%capac(2) = sib%prog%capac(2) + sib%diag%www_inflow * dtt


end subroutine hydro_snow

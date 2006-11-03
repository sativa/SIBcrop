!----------------------------------------------------------------------

subroutine compact_snow(sib,sib_loc)

!----------------------------------------------------------------------
!
!   Based on CLM subroutine clm_compact
!
!   CLM web info : http://clm.gsfc.nasa.gov
!
!   Description:
!   Three metamorphisms of changing snow characteristics are 
!   implemented, i.e., destructive, overburden and melt. The
!   treatments of the former two are from SNTHERM.89 and SNTHERM.99 
!   (1991, 1999). The contribution due to melt metamorphism is 
!   simply taken as a ratio of snow ice fraction after the melting 
!   versus before the melting.
!
!   Revision History:
!   15 September 1999: Yongjiu Dai; initial code
!   15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!   30 January   2002: Ian Baker, SiB integration

use kinds
use sibtype

use sib_const_module, only: &
    denice, &
    denh2o, &
    dtt,    &
    dti

use physical_parameters, only: &
    tice

implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

type(sib_local_vars)     ,intent(inout) :: sib_loc
! variables local to SiB

!----------------------------------------------------------------------  


!...local variables

integer(kind=int_kind)        :: j
real(kind=dbl_kind),parameter :: c2 = 23.0e-3   ! (m^3 kg^-1)
real(kind=dbl_kind),parameter :: c3 = 2.77e-6   ! (sec^-1)
real(kind=dbl_kind),parameter :: c4 = 0.04      ! (K^-1)
real(kind=dbl_kind),parameter :: c5 = 2.0       ! 
real(kind=dbl_kind),parameter :: dm = 100.0     ! upper limit on 
! destructive metamorphism compaction (kg m^-3)
real(kind=dbl_kind),parameter :: eta0 = 9.0e5   ! viscosity coefficient
!  (kg m^-3)

real(kind=dbl_kind) :: burden    ! pressure of overlying snow (kg m^-2)
real(kind=dbl_kind) :: wx        ! water mass (ice + liquid) (kg m^-2)
real(kind=dbl_kind) :: void      ! void = 1 - vol_ice - vol_liquid
real(kind=dbl_kind) :: bi        ! partial density of ice (kg m^-3)
real(kind=dbl_kind) :: fi        ! fraction of ice relative to total
!  water content
real(kind=dbl_kind) :: delt      ! snow sib%prog%td - tice (K)
real(kind=dbl_kind) :: dexpf     ! expf = exp(-c4*(tice-sib%prog%td))
real(kind=dbl_kind) :: ddz1      ! rate of settling snowpack due to 
!  destructive metamorphism
real(kind=dbl_kind) :: ddz2      ! rate of compaction of snowpack due
!  to overburden
real(kind=dbl_kind) :: ddz3      ! rate of compaction of snowpack due
!  to melt
real(kind=dbl_kind) :: pdzdtc    ! nodal rate of change in fractonal
!  thickness due to compaction 
!  (fraction sec^-1)


!----------------------------------------------------------------------

!    if ( sib%prog%nsl == 0 ) return

    burden = 0.0

    do j=sib%prog%nsl+1,0

        wx = sib%prog%www_ice(j) + sib%prog%www_liq(j)
        void = 1.0 - (sib%prog%www_ice(j)/denice + sib%prog%www_liq(j)/ &
            denh2o)/sib%prog%dz(j)

        !...disallow compaction for water saturated node and lower ice lens node
        if(void <= 0.001  .or.  sib%prog%www_ice(j) <= 0.1) then
            burden = burden + wx
            cycle
        endif

        bi = sib%prog%www_ice(j)/sib%prog%dz(j)
        fi = sib%prog%www_ice(j)/wx
        delt = tice - sib%prog%td(j)
        dexpf = exp(-c4*delt)

        !...settling as a result of desctructive metamorphism

        ddz1 = -c3*dexpf
        if(bi > dm) ddz1 = ddz1*exp(-46.e-3*(bi-dm))

        !...liquid water term

        if(sib%prog%www_liq(j) > 0.01*sib%prog%dz(j)) ddz1 = ddz1*c5

        !...compaction due to overburden

        ddz2 = -burden*exp(-0.08*delt-c2*bi)/eta0

        !...compaction occurring durin melt

        if(sib_loc%imelt(j) == 1 )then
            ddz3 = -1.*dti * max(0.0_dbl_kind, &
                (sib_loc%frac_iceold(j) - fi)/sib_loc%frac_iceold(j))
        else
            ddz3 = 0.0
        endif

        !...time rate of fractional change in dz (units of sec-1)

        pdzdtc = ddz1 + ddz2 + ddz3

        !...change in dz due to compaction

        sib%prog%dz(j) = sib%prog%dz(j) * (1.0+pdzdtc*dtt)

        !...pressure of overlying snow

        burden = burden + wx

    enddo ! nsnow loop


end subroutine compact_snow

!=======================================================================
subroutine phen_mapper(                 &
                   lat,                 &
                   DOY,                 &
                   fVCover,             &
                   ChiL,                &
                   LTran,               &
                   LRef,                &
                   MorphTab,            &
                   AeroVar,             &
                   LAIgrid,             &
                   fVCovergrid,         &
                   TimeVar )
!=======================================================================
! calculates time dependant boundary condition variables for SiB.

use kinds

!----------------------------------------------------------------------
!type(sib_t), intent(inout) :: sib
!----------------------------------------------------------------------  

IMPLICIT NONE

! begin input variables
!real(kind=dbl_kind) :: lai_in   ! input LAI from phenology model
real lat         ! center latitude of grid cell
real(kind=dbl_kind) :: fVCover     !
real(kind=dbl_kind) :: ChiL	 !
real LTran(2,2)  !
real LRef(2,2)   !

! begin input biome dependant, physical morphology variables
type biome_morph_var
   real (kind=real_kind) :: zc        ! Canopy inflection height (m)
   real (kind=real_kind) :: LWidth    ! Leaf width
   real (kind=real_kind) :: LLength   ! Leaf length	
   real (kind=real_kind) :: LAImax    ! Maximum LAI0.07

   real (kind=real_kind) :: stems     ! Stem area index
!   real (kind=real_kind) :: NDVImax   ! Maximum NDVI
!   real (kind=real_kind) :: NDVImin   ! Minimum NDVI
!   real (kind=real_kind) :: SRmax     ! Maximum mple ratio
!   real (kind=real_kind) :: SRmin     ! Minimum simple ratio
end type biome_morph_var

type(biome_morph_var) MorphTab

! begin input aerodynamic parameters
type aero_var
   real (kind=real_kind) :: zo	      ! Canopy roughness coeff 
   real (kind=real_kind) :: zp_disp  ! Zero plane displacement
   real (kind=real_kind) :: RbC      ! RB Coefficient
   real (kind=real_kind) :: RdC      ! RC Coefficient
end type aero_var

type(aero_var) AeroVar(50,50) ! aerodynamic interpolation tables

real (kind=real_kind) :: LAIgrid(50)       ! grid of LAI values for lookup table
real (kind=real_kind) :: fVCovergrid(50)   ! grid of fVCover values for 
                                         !  interpolation table

! begin time dependant, output variables
type time_dep_var
   real (kind=real_kind) :: fPAR    ! Canopy absorbed fraction of PAR
   real (kind=real_kind) :: LAI     ! Leaf-area index
   real (kind=real_kind) :: Green   ! Canopy greeness fraction of LAI
   real (kind=real_kind) :: zo      ! Canopy roughness coeff 
   real (kind=real_kind) :: zp_disp ! Zero plane displacement
   real (kind=real_kind) :: RbC     ! RB Coefficient (c1)
   real (kind=real_kind) :: RdC     ! RC Coefficient (c2)
   real (kind=real_kind) :: gmudmu  ! Time-mean leaf projection
end type time_dep_var

type(time_dep_var) TimeVar

! begin internal variables
real (kind=real_kind) ::DOY         ! Day of Year (DOY) of ndvi input map
real (kind=real_kind) ::prevfPAR    ! previous month's fPAR value
real(kind=real_kind), parameter :: fPARmax=0.95
                 ! Maximum possible FPAR corresponding to 98th percentile
real(kind=real_kind), parameter :: fPARmin=0.01
                 ! Minimum possible FPAR corresponding to 2nd percentile
!     For more information on fPARmin and fPARmax, see
!     Sellers et al. (1994a, pg. 3532); Los (1998, pg. 29, 37-39)

   !-----------------------------------------------------------------------
   ! Calculate time dependant variables
   !-----------------------------------------------------------------------


   call laigrn_phen (TimeVar%fPAR, prevfPAR, fPARmax, fVCover,         &
        	    MorphTab%stems, MorphTab%LAImax, TimeVar%Green,   &
        	    TimeVar%LAI)


   ! Interpolate to calculate aerodynamic, time varying variables
   call AeroInterpolate (TimeVar%LAI, fVCover, LAIgrid,fVCovergrid,   &
                         AeroVar, TimeVar%zo, TimeVar%zp_disp,        &
                         TimeVar%RbC, TimeVar%RdC)

   ! Calculate mean leaf orientation to par flux (gmudmu)
   call gmuder (lat, DOY, ChiL, TimeVar%gmudmu)


   ! recalculate fPAR adjusting for Sun angle, vegetation cover fraction,
   ! and greeness fraction, and LAI
   call aparnew (TimeVar%LAI, TimeVar%Green, LTran, LRef,   &
                 TimeVar%gmudmu, fVCover, TimeVar%fPAR,     &
                 fPARmax, fPARmin)


   return
end subroutine phen_mapper

!subrountine averageapar deleted- EL

!=======================================================================
subroutine laigrn_phen (fPAR,fPARm,fPARmax,fVCover,stems, LAImax,Green,LAI)
!=======================================================================
! calculate leaf area index (LAI) and greenness fraction (Green) from fPAR. 
! LAI is linear with vegetation fraction and exponential with fPAR.
! See Sellers et al (1994), Equations 7 through 13.

use kinds
                                                                      
implicit none

! begin input variables
real LAI      ! area average total leaf area index
real fPAR     ! fraction of PAR absorbed by plants at current time
real fPARm    ! fraction of PAR absorbed by plants at previous time
real fPARmax  ! maximum possible FPAR corresponding to 98th percentile
real(kind=dbl_kind) :: fVCover  ! vegetation cover fraction
real stems    ! stem area index for the specific biome type
real LAImax   ! maximum total leaf area index for specific biome type

! begin output variables
real Green    ! greeness fraction of the total leaf area index


! begin internal variables
real LAIg     ! green leaf area index at current time
real LAIgm    ! green leaf area index at previous time
real LAId     ! dead leaf area index at current time

! Calculate current and previous green leaf area index (LAIg and LAIgm):
! LAIg is log-linear with fPAR.  Since measured fPAR is an area average, 
! divide by fVCover to get portion due to vegetation.  Since fVCover can
! be specified, check to assure that calculated fPAR does not exceed fPARMax.

!EL-using phen_LAI to estimate fPAR
LAI=LAI
LAIg=LAI
fPAR=(1-exp((LAIg*alog(1-fPARmax))/LAImax))*real(fvcover)

   if(fPAR/fVCover.ge.fPARmax) then
      LAIg=LAImax
   else
     LAIg=alog(1.-fPAR/real(fVCover))*LAImax/alog(1-fPARmax)
   endif

   if(fPARm/fVCover.ge.fPARmax) then
      LAIgm=LAImax
   else
      LAIgm=alog(1.-fPARm/real(fVCover))*LAImax/alog(1-fPARmax)
   endif

print*,fpar,LAIg
   ! Calculate dead leaf area index (LAId):
   ! If LAIg is increasing or unchanged, the vegetation is in growth mode.
   ! LAId is then very small (very little dead matter).
   ! If LAIg is decreasing, the peak in vegetation growth has passed and
   ! leaves have begun to die off.  LAId is then half the change in LAIg,
   ! assuming half the dead leaves fall off.

   !     Growth mode dead leaf area index:
   if (LAIg.ge.LAIgm) LAId=0.0001

   !     die-off (post peak growth) dead leaf area index:
   if (LAIg.lt.LAIgm) LAId=0.5*(LAIgm-LAIg)


   ! Calculate greeness fraction (Green):
   ! Greeness fraction=(green leaf area index)/(total leaf area index)
   Green=LAIg/(LAIg+LAId+stems)

   return                                                                    
end subroutine laigrn_phen




!--------------------------------------------------------------
subroutine sibdrv_interp(sib, time)
!--------------------------------------------------------------
! This subroutine interpolates the sibdrv forcing surface meteorological
! variables between their read times
!
! Modifications:
!  Kevin Schaefer changed dayflag to cosz_min when checking if sw_rad in scaled (8/12/04)
!  Kevin Schaefer removed checks on sw_dwn between 1-5 (created light at night) (8/12/04)
!  Kevin Schaefer moved pressure conversion to read_drvr routines (8/16/04)
!  Kevin Schaefer move conversion from dew point to humidity to read_drvr routines (8/17/04)
!  Kevin Schaefer synched facsibdrv with 1 step earlier driver data read (8/18/04)
!  Kevin Schaefer deleted unused variables & commented code (11/15/04)
!--------------------------------------------------------------

use kinds
use sibtype
use timetype
use sib_const_module, only:  &
    rgfac,         &
    subcount,      &
    cosz_min
use physical_parameters, only:  kappa
!
implicit none
!
! Input parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time
!
! local variables
real(kind=dbl_kind):: facsibdrv  ! scaling factor between driver data points
real(kind=dbl_kind):: temp       ! temporary variable for testing
integer(kind=int_kind) :: i      ! index
!
! initialize CO2 partial pressure in boundary layer
    do i = 1,subcount
        sib(i)%prog%pco2m = 35.0
        sib(i)%prog%psb = 50.0
    enddo
!
! calculate cosine zenith angle
    call zenith_angle( time%hour, sib(:)%stat%cosz )
!
! scale downwelling driver SW radiation
    do i = 1,subcount
!
! scale SW radiation to cosine of zenith angle
	if((sib(i)%stat%cosz>cosz_min).and.(sib(i)%stat%coszbar.ne.0.0_dbl_kind))then
            	sib(i)%prog%sw_dwn = sib(i)%prog%sw_dwn2 *  max(sib(i)%stat%cosz-cosz_min,0.0_dbl_kind )/sib(i)%stat%coszbar
        else
            sib(i)%prog%sw_dwn=0.
        endif

!
! make sure downwelling SW radiation is positive
        if(sib(i)%prog%sw_dwn<0.) sib(i)%prog%sw_dwn = 0.
    enddo
!
! split downwelling SW radiation to visible and NIR (direct and diffuse)
    call raddrv(subcount,sib(:)%prog%sw_dwn,sib(:)%stat%cosz,sib%prog%radvbc,  &
        sib%prog%radvdc,sib%prog%radnbc,sib%prog%radndc)
!
! calculate driver data interpolation factor
    facsibdrv = 1. -float(mod(time%sec_day,time%driver_step)) / float(time%driver_step)
!
! loop through all land points
    do i = 1,subcount
!
! interpolate driver temperature
        sib(i)%prog%tm = facsibdrv*sib(i)%prog%tm1 + (1. - facsibdrv) *  &
            sib(i)%prog%tm2
!
! interpolate driver humidity
        sib(i)%prog%sh = facsibdrv*sib(i)%prog%sh1 + (1. - facsibdrv) *  &
            sib(i)%prog%sh2
!
! interpolate driver surface pressure
        sib(i)%prog%ps = facsibdrv*sib(i)%prog%ps1 +   &
            (1. - facsibdrv) * sib(i)%prog%ps2
!
! interpolate driver wind speed
        sib(i)%prog%spdm = facsibdrv*sib(i)%prog%spdm1 +   &
            (1. - facsibdrv) * sib(i)%prog%spdm2
        sib(i)%prog%spdm = MAX(sib(i)%prog%spdm,1.0_dbl_kind)
!
! interpolate driver large scale precipitation
        sib(i)%prog%lspr =  sib(i)%prog%lspr1 / time%driver_step
!
! interpolate driver cumulus or convective precipitation
        sib(i)%prog%cupr =  sib(i)%prog%cupr1 / time%driver_step
!
! interpolate driver downwelling longwave radiation
        sib(i)%prog%dlwbot = facsibdrv*sib(i)%prog%dlwbot1 +   &
            (1. - facsibdrv) * sib(i)%prog%dlwbot2
!
! calculate reference level potential temperature and some related parameters
        sib(i)%prog%bps(1) = (0.001*sib(i)%prog%ps)**kappa
        sib(i)%prog%bps(2) = (0.001*(sib(i)%prog%ps-sib(i)%prog%psb))**kappa
        sib(i)%prog%thm = sib(i)%prog%tm / sib(i)%prog%bps(1)
!
! calculate reference level air density
        sib(i)%prog%ros = rgfac * sib(i)%prog%ps / sib(i)%prog%tm

    enddo  ! land point loop
!
! initialize Canopy Air Space humidity
   !ogl...changed to an explicit do loop
    if( time%sec_tot == time%init_second ) then
        print*, 'sibdrv_interp: init CAS humidity'
        do i=1,subcount
!	print*, sib(i)%prog%sh,sib(i)%prog%ps,sib(i)%prog%spdm,sib(i)%prog%tm
            sib(i)%prog%sha = sib(i)%prog%sh
        enddo
    endif
!
end subroutine sibdrv_interp
!
!---------------------------------------------------------------------
subroutine raddrv(nsib,swdown,sunang,radvbc,radvdc,radnbc,radndc)
!---------------------------------------------------------------------
!               radiation radive code to use the downward sw at bottom 
!               and the formulation to estimate radvbc,radvdc, radndc,
!               radndc
!---------------------------------------------------------------------

use kinds
implicit none

integer(kind=int_kind) :: nsib, i
real(kind=real_kind) ::  &
    cloud,          &
    difrat,         &
    vnrat
real(kind=dbl_kind) ::  &
    swdown(nsib),   &
    sunang(nsib),   &
    stemp,           &
    localcosz
real(kind=dbl_kind) ::  &
    radvbc(nsib),   &
    radvdc(nsib),   &
    radnbc(nsib),   &
    radndc(nsib)

real(kind=dbl_kind),parameter :: c1 = 580.
real(kind=dbl_kind),parameter :: c2 = 464.
real(kind=dbl_kind),parameter :: c3 = 499.
real(kind=dbl_kind),parameter :: c4 = 963.
real(kind=dbl_kind),parameter :: c5 = 1160.

    do i=1,nsib
        localcosz = max( 0.001_dbl_kind, sunang(i) )
        stemp = swdown(i)
        stemp = MAX(stemp,0.01_dbl_kind)

        cloud = (c5 * localcosz - stemp) / (c4 * localcosz)                   
        cloud = max(cloud,0.)                                                
        cloud = min(cloud,1.)                                                  

        difrat = 0.0604 / ( sunang(i)-0.0223 + 1.0e-10 ) + 0.0683
        if ( difrat < 0. ) difrat = 0.
        if ( difrat > 1. ) difrat = 1.
        difrat = difrat + ( 1. - difrat ) * cloud

        vnrat = ( c1 - cloud*c2 ) / ( ( c1 - cloud*c3 ) + ( c1 - cloud*c2 ) )

        radvbc(i) = (1.-difrat)*vnrat*stemp
        radvdc(i) = difrat*vnrat*stemp
        radnbc(i) = (1.-difrat)*(1.-vnrat)*stemp
        radndc(i) = difrat*(1.-vnrat)*stemp
    enddo

end subroutine raddrv

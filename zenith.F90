!-------------------------------------------------------------------------------
subroutine init_solar_dec( time )
!-------------------------------------------------------------------------------
! initializes the declination of the Sun
!
! FUNCTIONS CALLED:
!  rt_asc
!
! MODIFICATIONS:
!  Created by Kevin Schaefer (3/4/03)
!  Kevin Schaefer moved time indep calculations to sib_main_mod (3/5/03)
!-----------------------------------------------------------------------

use kinds
use timetype
use sib_const_module
use physical_parameters
implicit none

! parameters
type(time_struct), intent(in) :: time

! local variables
real(kind=real_kind) :: t1       ! 1st factor to determine longitude of earth (lonearth)
real(kind=real_kind) :: t2       ! 2nd factor to determine longitude of earth (lonearth)
real(kind=real_kind) :: t3       ! 3rd factor to determine longitude of earth (lonearth)
real(kind=real_kind) :: t4       ! 4th factor to determine longitude of earth (lonearth)
integer(kind=int_kind) :: iday    ! day of year since vernal equinox variable
real(kind=real_kind) :: rt_asc

    ! lon of Earth from equinox at start of simulation
    iday = time%doy - eqnx
    if ( iday < 0 ) iday = iday + time%days_per_year
    lonearth=0.0
    if ( iday /= 0 ) then
      do while ( iday > 0 )
        iday = iday - 1
        t1 = rt_asc( lonearth ) * pidaypy
        t2 = rt_asc( lonearth+t1*.5 ) * pidaypy
        t3 = rt_asc( lonearth+t2*.5 ) * pidaypy
        t4 = rt_asc( lonearth+t3 ) * pidaypy
        lonearth = lonearth + (t1+2.*(t2+t3)+t4) / 6.
      enddo
    endif

end subroutine init_solar_dec


!-------------------------------------------------------------------------------
subroutine solar_dec( time )
!-------------------------------------------------------------------------------
! Calculates the declination of the Sun
!
! FUNCTIONS CALLED:
!  rt_asc
!
! MODIFICATIONS:
!  Created by Kevin Schaefer (3/4/03)
!  Kevin Schaefer moved time indep calculations to sib_main_mod (3/5/03)
!-----------------------------------------------------------------------

use kinds
use timetype
use sib_const_module
use physical_parameters
implicit none

! parameters
type(time_struct), intent(in) :: time

! local variables
real(kind=real_kind) :: t1  ! 1st factor to determine longitude of earth (lonearth)
real(kind=real_kind) :: t2  ! 2nd factor to determine longitude of earth (lonearth)
real(kind=real_kind) :: t3  ! 3rd factor to determine longitude of earth (lonearth)
real(kind=real_kind) :: t4  ! 4th factor to determine longitude of earth (lonearth)
real(kind=real_kind) :: rt_asc

    ! reset lon of Earth from equinox
    if( time%doy == eqnx ) lonearth = 0.0

    ! Increment Longitude of Earth
    t1 = rt_asc( lonearth ) * pidaypy
    t2 = rt_asc( lonearth+t1*.5 ) * pidaypy
    t3 = rt_asc( lonearth+t2*.5 ) * pidaypy
    t4 = rt_asc( lonearth+t3 ) * pidaypy
    lonearth = lonearth + (t1+2.*(t2+t3)+t4) / 6.

    ! Calculate the sine and cosine of Solar declination
    sin_dec = sin( decmax*(pi/180.) ) * sin( lonearth )
    cos_dec = sqrt( 1. - sin_dec * sin_dec )


end subroutine solar_dec


!-------------------------------------------------------------------------------
function rt_asc( ang )
!-------------------------------------------------------------------------------
! calculates correction for longitude (right ascension) of the Earth
! from vernal equinox based on the angle around Sun traversed
! since beginning of the year

    use kinds
    use sib_const_module
    use physical_parameters
    implicit none

real(kind=real_kind) :: rt_asc ! right ascension correction
real(kind=real_kind) :: ang    ! angle from beginning of year (radians)

! local variables
real(kind=real_kind) :: reccn
real(kind=real_kind) :: perhlr

    ! rt ascension correction based on Earth's orbit
    reccn  = 1. / (1. - eccn * eccn ) **1.5
    perhlr = perhl * (pi/180.)
    rt_asc = reccn * (1. - eccn * cos( ang - perhlr ) ) **2

end function rt_asc


!=======================================================================
subroutine mean_zenith_angle( sib, time )
!=======================================================================      
! Calculates a mean cosine zenith angle for civil twilight and daylight  
! hours between driver data points to scale driver radiation data.
! Daylight is defined as zenith angle<=90. (cosz>=0)
! Civil twilight is defined as 96<zenith angle<90 (-.1045<cosz<0).
! Night is defined as zenith angle>=96. (cosz<=-.1045)
!
! Modifications:
!  Kevin Schaefer created mean_zenith_angle subroutine (5/31/01)
!  Kevin Schaefer changed coszbar calculation to daylight only (5/31/01)
!  Kevin Schaefer included civil twilight bias in mean zenith angle (7/21/01)
!  Kevin Schaefer deleted use sib_main_module; only variable list (3/4/03)
!----------------------------------------------------------------------

use kinds
use sibtype
use timetype
use sib_const_module
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time

! local variables
integer i         ! sib point index
integer n         ! time of day index
integer nsteps    ! number SiB time steps per driver data time step
real loctofday    ! local version of time of day (GMT, in hours)

    ! Clear out arrays
    sib(:)%stat%coszbar = 0.
    sib(:)%stat%dayflag = 0.

    ! calculate number of sib time steps per driver data time step
    nsteps = time%driver_step / time%dtsib

    ! Integrate cosz over driver data time step (tofday to tofday+dtsibmetin)
    loctofday = time%hour
    do n = 1, nsteps

        ! calculate zenith angle at local time of day
        call zenith_angle ( loctofday, sib(:)%stat%cosz )

        ! check for civil twilight (-.1045<loccosz<0.) or daylight (loccosz>0)
        ! add loccosz to sum, include bias to account for civil twilight
        do i = 1, subcount
            if( sib(i)%stat%cosz > cosz_min ) then
                sib(i)%stat%dayflag = sib(i)%stat%dayflag + 1.
                sib(i)%stat%coszbar = sib(i)%stat%coszbar +  &
                    sib(i)%stat%cosz - cosz_min
            endif
        enddo

        ! increment local time of day
        loctofday = loctofday + real(time%dtsib)/3600.
    enddo

    ! Calculate mean cosz during civil twilight and daylight hours
    do i = 1, subcount
        if( sib(i)%stat%coszbar > 0. )  &
            sib(i)%stat%coszbar = sib(i)%stat%coszbar / nsteps
    enddo

end subroutine mean_zenith_angle


!=======================================================================
subroutine zenith_angle ( hour, cosz )
!=======================================================================      
! calculates the zenith angle for an array of lat/lon points.
! The cos_hour variable accounts for the difference between GMT and local 
!  time
!
! Modifications:
!  Kevin Schaefer created zenith_angle subroutine (5/31/01)
!  Kevin Schaefer corrected hour angle calculation (5/31/01)
!  Kevin Schaefer removed minimum cosz limit to include twilight (7/21/01)
!----------------------------------------------------------------------

use kinds
use sib_const_module
implicit none

! parameters
real(kind=real_kind), intent(in) :: hour
real(kind=dbl_kind), dimension(subcount), intent(out) :: cosz

! internal variables
real pid180      ! conversion from degrees to radians
real sinlat      ! sine of latitude
real coslat      ! cosine of latitude
real hrang       ! hour angle; longitude of Sun from Greenwhich meridian
real cos_hour    ! cosine delta longitude between Sun & SiB point
integer i        ! sib point index

    ! calculate conversion from degrees to radians
    pid180 = 3.14159/180.

    ! Calculate hour angle (longitude of Sun from Greenwhich meridian)
    hrang = (12. - hour) * 360. / 24.

    ! Calculate cosine of solar zenith angle for each SiB point
    do i = 1,subcount
        cos_hour = cos( pid180 * ( hrang - lonsib(subset(i)) ) )
        sinlat = sin( pid180 * latsib(subset(i)) )
        coslat = cos( pid180 * latsib(subset(i)) )
        cosz(i) = coslat * cos_dec * cos_hour + sinlat * sin_dec

    enddo

end subroutine zenith_angle


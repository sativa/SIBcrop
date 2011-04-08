subroutine sibdrv_read_ncep2( sib, time )

!****--------------------------------------------------------------------
!    This subroutines reads the forcing surface meteorological data 
!    for the next time step.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.

!    precip, snowfall and radiation data are aggregated data over 6 hours
!    They have been stored such that the data at time x represent the
!    aggregation over the preceeding 6 hours.
!    Thus data at 
!     0 hours are the sum of the values 18-24 hours of the previous day
!     6 hours are the sum of the values  0- 6 hours of the current day
!    12 hours are the sum of the values  6-12 hours of the current day
!    18 hours are the sum of the values 12-18 hours of the current day
!
!    Consequently new data are read for all variables at the same time.
!    For all point data variables, e.g. temp., the data are for
!    now + 6 hours. For precip and radiation data the new data
!    are the aggregation over the next six hours.
!
! Modifications:
!  Kevin Schaefer moved conversion pascals to millibars from sibdrv_interp to here (8/16/04)
!  Kevin Schaefer changed from tdew1/2 to sh1/2 because thats whats there (8/17/04)
!****--------------------------------------------------------------------


use sib_const_module, only:   &
    nsib,               &
    latsib,             &
    lonsib,             &
    subset,             &
    subcount

use sib_io_module, only:   &
    dr_format,           &
    driver_id

use physical_parameters, only:              &
    kapa => kappa,   &
    pi

use kinds
use sibtype
use timetype

use netcdf
use typeSizes

#include "nc_util.h"

type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time

integer(kind=int_kind) :: i, iyear, imon, iday, idoy, ihour, imin

character*80 filename
integer(kind=int_kind) :: ncyid,ncmid,nctdid,ncdoyid,nchid
integer(kind=int_kind), dimension(2) :: mstart,mcount

integer(kind=int_kind) :: nct2mid ! Total Cloud Cover
integer(kind=int_kind) :: nctccid ! Total Cloud Cover
integer(kind=int_kind) :: ncswdid ! Surface solar radiation downwards           
integer(kind=int_kind) :: ncldwid ! Surface thermal radiation downwards
integer(kind=int_kind) :: ncuwdid ! U-wind at 10 m
integer(kind=int_kind) :: ncvwdid ! V-wind at 10 m
integer(kind=int_kind) :: ncshid  ! humidity at 2 m
integer(kind=int_kind) :: ncsfpid ! Log Surface Pressure
integer(kind=int_kind) :: nclspid ! Large Scale Precipitation
integer(kind=int_kind) :: nccvpid ! Convective Precipitation
integer(kind=int_kind) :: ncsflid ! Snow Fall

real(kind=real_kind), dimension(nsib) :: t2m ! Total Cloud Cover
real(kind=real_kind), dimension(nsib) :: tcc ! Total Cloud Cover
real(kind=real_kind), dimension(nsib) :: swd ! Surface solar radiation downwards           
real(kind=real_kind), dimension(nsib) :: ldw ! Surface thermal radiation downwards
real(kind=real_kind), dimension(nsib) :: sh  ! humidity at 2 m
real(kind=real_kind), dimension(nsib) :: sfp ! Log Surface Pressure
real(kind=real_kind), dimension(nsib) :: lsp ! Large Scale Precipitation
real(kind=real_kind), dimension(nsib) :: cvp ! Convective Precipitation

integer(kind=int_kind) :: status
real(kind=real_kind) :: xyear,xmonth,xdoy,xday,xhour
real(kind=real_kind), dimension(nsib) :: xx,uwd,vwd

character(len=13) :: subname
data subname/'sibdrv_read '/


   !*** Storing previous time steps data
    do i=1,subcount
        sib(i)%prog%ps1       = sib(i)%prog%ps2
        sib(i)%prog%tm1       = sib(i)%prog%tm2
        sib(i)%prog%tcc1      = sib(i)%prog%tcc2
        sib(i)%prog%sh1       = sib(i)%prog%sh2
        sib(i)%prog%spdm1     = sib(i)%prog%spdm2
        sib(i)%prog%lspr1     = sib(i)%prog%lspr2
        sib(i)%prog%cupr1     = sib(i)%prog%cupr2
        sib(i)%prog%dlwbot1   = sib(i)%prog%dlwbot2
        sib(i)%prog%sw_dwn1 = sib(i)%prog%sw_dwn2
    enddo

    ! switch files if needed
    if ( time%switch_driver ) then
        status = nf90_close( driver_id )

        write( filename, dr_format ) time%driver_year, time%driver_month
        CHECK( nf90_open( trim(filename), nf90_nowrite, driver_id ) )
    endif

    ! Read new driver data
    !print*, subname,tau,nextsecond,nextday,nextdoy,       &
    !    nmonth,nyear,nextmonth,nextyear

    ! check time values in driver data file
    ENSURE_VAR( driver_id,   'year', ncyid )
    ENSURE_VAR( driver_id,   'month',ncmid )
    ENSURE_VAR( driver_id,   'doy',  ncdoyid )
    ENSURE_VAR( driver_id,   'day',  nctdid )
    ENSURE_VAR( driver_id,   'hour', nchid )

    ! read time
    mstart(1) = time%driver_recnum
    CHECK( nf90_get_var( driver_id, ncyid,    xyear, mstart(1:1) ) )
    CHECK( nf90_get_var( driver_id, ncmid,   xmonth, mstart(1:1) ) )
    CHECK( nf90_get_var( driver_id, ncdoyid,   xdoy, mstart(1:1) ) )
    CHECK( nf90_get_var( driver_id, nctdid,    xday, mstart(1:1) ) )
    CHECK( nf90_get_var( driver_id, nchid,    xhour, mstart(1:1) ) )

    ihour=xhour
    iday =xday
    idoy =xdoy
    imon =xmonth
    iyear=xyear
    imin=0

    print*,subname,'Time level in file: ',ihour,iday,idoy,imon,iyear

    !* Get variable id's
    ENSURE_VAR( driver_id, 't2m', nct2mid ) ! Temperature at 2 m
    ENSURE_VAR( driver_id, 'tcc', nctccid ) ! Total Cloud Cover
    ENSURE_VAR( driver_id, 'swd', ncswdid ) ! Surface solar rad downwards
    ENSURE_VAR( driver_id, 'lwd', ncldwid ) ! Surface thermal rad down
    ENSURE_VAR( driver_id, 'uwd', ncuwdid ) ! U-wind at 10 m
    ENSURE_VAR( driver_id, 'vwd', ncvwdid ) ! V-wind at 10 m
    ENSURE_VAR( driver_id, 'shum', ncshid ) ! humidity at 2 m
    ENSURE_VAR( driver_id, 'sfp', ncsfpid ) ! Log Surface Pressure
    ENSURE_VAR( driver_id, 'lsp', nclspid ) ! Large Scale Precipitation
    ENSURE_VAR( driver_id, 'cvp', nccvpid ) ! Convective Precipitation
    ENSURE_VAR( driver_id, 'sfl', ncsflid ) ! Snow Fall

    !* get data
    mstart=(/1,time%driver_recnum/); mcount=(/nsib,1/)
    CHECK( nf90_get_var( driver_id, nct2mid, t2m, mstart, mcount ) ) !Temperature at 2 m
    CHECK( nf90_get_var( driver_id, nctccid, tcc, mstart, mcount ) ) !Total Cloud Cover
    CHECK( nf90_get_var( driver_id, ncswdid, swd, mstart, mcount ) ) !Surface solar rad downwards
    CHECK( nf90_get_var( driver_id, ncldwid, ldw, mstart, mcount ) ) !Surface thermal rad downwards
    CHECK( nf90_get_var( driver_id, ncuwdid, uwd, mstart, mcount ) ) ! U-wind at 10 m
    CHECK( nf90_get_var( driver_id, ncvwdid, vwd, mstart, mcount ) ) ! V-wind at 10 m
    CHECK( nf90_get_var( driver_id, ncshid, sh, mstart, mcount ) )   ! humidity at 2 m
    CHECK( nf90_get_var( driver_id, ncsfpid, sfp, mstart, mcount ) ) ! Surface Pressure
    CHECK( nf90_get_var( driver_id, nclspid, lsp, mstart, mcount ) ) ! Large Scale Precipitation
    CHECK( nf90_get_var( driver_id, nccvpid, cvp, mstart, mcount ) ) ! Convective Precipitation
    CHECK( nf90_get_var( driver_id, ncsflid, xx, mstart, mcount ) )  ! Snow Fall


    do i=1,subcount
        ! pull out landpoints in subdomain
        sib(i)%prog%tm2 = t2m(subset(i))
        sib(i)%prog%tcc2 = tcc(subset(i))
        sib(i)%prog%sw_dwn2 = swd(subset(i))
        sib(i)%prog%dlwbot2 = ldw(subset(i))
        sib(i)%prog%sh2 = sh(subset(i))
        sib(i)%prog%ps2 = sfp(subset(i))
        sib(i)%prog%lspr2 = lsp(subset(i))
        sib(i)%prog%cupr2 = cvp(subset(i))
    
        ! scale radiation to w/m2
        !sib(i)%prog%sw_dwn2 = sib(i)%prog%sw_dwn2/time%driver_step
        !sib(i)%prog%dlwbot2 = sib(i)%prog%dlwbot2/time%driver_step
        if ( sib(i)%prog%sw_dwn2 < 0 ) sib(i)%prog%sw_dwn2 = 0.0
        if ( sib(i)%prog%dlwbot2 < 0 ) sib(i)%prog%dlwbot2 = 0.0
        
        ! convert total cloud cover to fraction
        sib(i)%prog%tcc2 = sib(i)%prog%tcc2 * 0.01

        ! 10 m wind
        sib(i)%prog%spdm2=SQRT(uwd(subset(i))*uwd(subset(i))+vwd(subset(i))*vwd(subset(i)))

        ! add snowfall to large scale precip and let SiB decide about snow.
        sib(i)%prog%lspr2 = sib(i)%prog%lspr2+xx(subset(i))
        ! convert to mm
        sib(i)%prog%lspr2 = (sib(i)%prog%lspr2-sib(i)%prog%cupr2)*time%driver_step
        sib(i)%prog%cupr2 = sib(i)%prog%cupr2*time%driver_step
        ! make sure precip > 0
        if ( sib(i)%prog%lspr2 < 0.0 ) sib(i)%prog%lspr2 = 0.0
        if ( sib(i)%prog%cupr2 < 0.0 ) sib(i)%prog%cupr2 = 0.0

        ! Conversion Pa -> hPa (pascals to millibars)
        sib(i)%prog%ps2 = sib(i)%prog%ps2 * 0.01
    enddo

    print*,subname,'New driver data read ',ihour,iday,imon,iyear
    print*,'------------------------------------------------------'
    print*,'Extrema of new input data'
    print*, minval(sib%prog%tm2      ),maxval(sib%prog%tm2  ),' Temperature'
    print*, minval(sib%prog%tcc2     ),maxval(sib%prog%tcc2),' Total cloudiness'
    print*, minval(sib%prog%sh2    ),maxval(sib%prog%sh2 ),' dew point'
    print*, minval(sib%prog%spdm2    ),maxval(sib%prog%spdm2),' Surface wind'
    print*, minval(sib%prog%ps2      ),maxval(sib%prog%ps2 ),' Pressure'
    print*, minval(sib%prog%dlwbot2  ),maxval(sib%prog%dlwbot2),  &
        ' Long wave down'
    print*, minval(sib%prog%lspr2    ),maxval(sib%prog%lspr2),' Large sc precip'
    print*, minval(sib%prog%cupr2    ),maxval(sib%prog%cupr2),' Convective '
    print*, minval(sib%prog%sw_dwn2),maxval(sib%prog%sw_dwn2),  &
        ' Short wave down'
    print*,'-----------------------------------------------------'


end subroutine sibdrv_read_ncep2
!
!-------------------------------------------------------
subroutine sibdrv_read_ncep1(sib, time)
!-------------------------------------------------------
! This subroutines reads the forcing surface meteorological
! data from the NCEP1 2x2 reanalysis for the next driver
! data time step.  If required, it closes the current
! month's data file and opens the next month's data file.
!
! precip, snowfall and radiation data are aggregated data over 6 hours
! They have been stored such that the data at time x represent the
! aggregation over the preceeding 6 hours.  Thus data at 
!     0 hours are the sum of the values 18-24 hours of the previous day
!     6 hours are the sum of the values  0- 6 hours of the current day
!    12 hours are the sum of the values  6-12 hours of the current day
!    18 hours are the sum of the values 12-18 hours of the current day
!
! Consequently new data are read for all variables at the same time.
! For all point data variables, e.g. temp., the data are for
! now + 6 hours. For precip and radiation data the new data
! are the aggregation over the next six hours.
!
! Modifications:
!  Kevin Schaefer created subroutine from sibdrv_read_ncep1(8/12/04)
!  Kevin Schaefer added check on zero humidity (8/16/04)
!-------------------------------------------------------
!
use sib_const_module, only:   &
    nsib,               &
    subset,             &
    subcount
use sib_io_module, only:   &
    dr_format,           &
    driver_id
use kinds
use netcdf
use typeSizes
use sibtype
use timetype

!
! define local variables
type(sib_t), dimension(subcount), intent(inout) :: sib  ! main sib variable tree
type(time_struct), intent(in) :: time         ! sibdrive time variables
integer(kind=int_kind) ::  i                  ! generic index
character*80 filename                                   ! netcdf driver data file name
integer(kind=int_kind) :: varid              ! netcdf variable id number
integer(kind=int_kind), dimension(2) :: mstart ! starting index location for driver data
integer(kind=int_kind), dimension(2) :: mcount ! number of driver data points to read
real(kind=real_kind), dimension(nsib) :: var  ! generic driver data variable
real(kind=real_kind), dimension(nsib) :: uwd  ! u (zonal) wind component
real(kind=real_kind), dimension(nsib) :: vwd  ! v (meridional) wind component
character(len=13) :: subname                            ! subroutine name
!
! set subroutine name
data subname/'read_ncep1'/
!
! print message
!print*, 'read new NCEP1 Driver data', time%hour, time%driver_recnum
!
!switch previous and next driver data
    do i=1,subcount
        sib(i)%prog%ps1     = sib(i)%prog%ps2
        sib(i)%prog%tm1     = sib(i)%prog%tm2
        sib(i)%prog%sh1     = sib(i)%prog%sh2
        sib(i)%prog%spdm1   = sib(i)%prog%spdm2
        sib(i)%prog%lspr1   = sib(i)%prog%lspr2
        sib(i)%prog%cupr1   = sib(i)%prog%cupr2
        sib(i)%prog%dlwbot1 = sib(i)%prog%dlwbot2
        sib(i)%prog%sw_dwn1 = sib(i)%prog%sw_dwn2
    enddo
!
! switch files if needed
    if (time%switch_driver) then
!
! close old file
        CHECK( nf90_close(driver_id) )
!
! new file name
        write(filename, dr_format) time%driver_year, time%driver_month
        print*, '\tswitch drvr to ', trim(filename)
!
! open new file
        CHECK( nf90_open(trim(filename), nf90_nowrite, driver_id) )
    endif
!
! set starting point for reading driver data file
    mstart=(/1,time%driver_recnum/); mcount=(/nsib,1/)
!
!-------------------------------------------------------
! Temperature at 2 m
    ENSURE_VAR(driver_id, 'tmp', varid)
    CHECK( nf90_get_var(driver_id, varid, var, mstart, mcount) )
!
! subgrid the driver data
    do i=1,subcount
        sib(i)%prog%tm2 = var(subset(i))
    enddo
!
!-------------------------------------------------------
! Surface solar radiation downwards
    ENSURE_VAR(driver_id, 'dswrf', varid)
    CHECK( nf90_get_var(driver_id, varid, var, mstart, mcount) )
!
! subgrid the driver data
    do i=1,subcount
        sib(i)%prog%sw_dwn2 = var(subset(i))
    enddo
!
!-------------------------------------------------------
! Surface thermal (infrared) radiation down
    ENSURE_VAR(driver_id, 'dlwrf', varid)
    CHECK( nf90_get_var(driver_id, varid, var, mstart, mcount) )
!
! subgrid the driver data
    do i=1,subcount
        sib(i)%prog%dlwbot2 = var(subset(i))
    enddo
!
!-------------------------------------------------------
! total wind speed
!
! U-wind at 10 m
    ENSURE_VAR(driver_id, 'ugrd', varid)
    CHECK( nf90_get_var(driver_id, varid, uwd, mstart, mcount) )
!
! V-wind at 10 m
    ENSURE_VAR(driver_id, 'vgrd', varid)
    CHECK( nf90_get_var(driver_id, varid, vwd, mstart, mcount) )
!
! subgrid the driver data
! combine winds into total wind
    do i=1,subcount
        sib(i)%prog%spdm2 = uwd(subset(i))*uwd(subset(i))+ vwd(subset(i))*vwd(subset(i))
        sib(i)%prog%spdm2 = sqrt(sib(i)%prog%spdm2)
    enddo
!
!-------------------------------------------------------
! specific humidity at 2 m
    ENSURE_VAR(driver_id, 'spfh', varid)
    CHECK( nf90_get_var(driver_id, varid, var, mstart, mcount) )
!
! subgrid the driver data
    do i=1,subcount
        sib(i)%prog%sh2 = var(subset(i))
        if(sib(i)%prog%sh2==0.) sib(i)%prog%sh2=1.e-4
    enddo
!
!-------------------------------------------------------
! Surface Pressure
    ENSURE_VAR(driver_id, 'pres', varid)
    CHECK( nf90_get_var(driver_id, varid, var, mstart, mcount) )
!
! subgrid the driver data
! convert pressure from pascals to millibars
    do i=1,subcount
        sib(i)%prog%ps2 = var(subset(i))*.01
    enddo
!
!-------------------------------------------------------
! Large Scale Precipitation
    ENSURE_VAR(driver_id, 'prate', varid)
    CHECK( nf90_get_var(driver_id, varid, var, mstart, mcount) )
!
! subgrid the driver data
! convert precipitation from kg m-2 s-2 to milimeters
    do i=1,subcount
        sib(i)%prog%lspr2 = var(subset(i))*time%driver_step
    enddo
!
!-------------------------------------------------------
! Convective Precipitation
    ENSURE_VAR(driver_id, 'cprat', varid)
    CHECK( nf90_get_var(driver_id, varid, var, mstart, mcount) )
!
! subgrid the driver data
! convert precipitation from kg m-2 s-2 to milimeters
    do i=1,subcount
        sib(i)%prog%cupr2 = var(subset(i))*time%driver_step
    enddo
!
!-------------------------------------------------------
! some checks on the driver data
    do i=1,subcount
!
! check for positive radiation
        if ( sib(i)%prog%sw_dwn2 < 0 ) sib(i)%prog%sw_dwn2 = 0.0
        if ( sib(i)%prog%dlwbot2 < 0 ) sib(i)%prog%dlwbot2 = 0.0
!
! check for positive precipitation
        if ( sib(i)%prog%lspr2 < 0.0 ) sib(i)%prog%lspr2 = 0.0
        if ( sib(i)%prog%cupr2 < 0.0 ) sib(i)%prog%cupr2 = 0.0
    enddo
!
end subroutine sibdrv_read_ncep1

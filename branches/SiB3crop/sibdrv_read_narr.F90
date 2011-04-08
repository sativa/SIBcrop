subroutine sibdrv_read_narr( sib, time )
!--------------------------------------------------------------
!    This subroutines reads the forcing surface meteorological data 
!    for the next time step.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.

!    precip, snowfall and radiation data are aggregated data over 6 hours
!    They have been stored such that the data at time x represent the
!    aggregation over the preceeding 6 hours.
!    Thus data at 
!     0 hours are the sum of the values 21-24 hours of the previous day
!     3 hours are the sum of the values  0- 3 hours of the current day
!     6 hours are the sum of the values  3- 6 hours of the current day
!     9 hours are the sum of the values  6- 9 hours of the current day
!    12 hours are the sum of the values  9-12 hours of the current day
!    15 hours are the sum of the values 12-15 hours of the current day
!    18 hours are the sum of the values 15-18 hours of the current day
!    21 hours are the sum of the values 18-21 hours of the current day
!
!    Consequently new data are read for all variables at the same time.
!    For all point data variables, e.g. temp., the data are for
!    now + 3 hours. For precip and radiation data the new data
!    are the aggregation over the next six hours.
!
! Modifications:
!  Kevin Schaefer moved conversion pascals to millibars from sibdrv_interp to here (8/16/04)
!  Kevin Schaefer changed from tdew1/2 to sh1/2 because thats whats there (8/17/04)
!--------------------------------------------------------------

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
        sib(i)%prog%sw_dwn1   = sib(i)%prog%sw_dwn2
    enddo

    ! switch files if needed
    if ( time%switch_driver ) then
        status = nf90_close( driver_id )

        write( filename, dr_format ) time%driver_year, time%driver_month
        print*,'      opening narr data: ',trim(filename)
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

!    print*,subname,'Time level in file: ',ihour,iday,idoy,imon,iyear

    !* Get variable id's
    ENSURE_VAR( driver_id, 'TMP_221_HTGL', nct2mid ) ! Temperature at 2 m
    ENSURE_VAR( driver_id, 'T_CDC_221_EATM', nctccid ) ! Total Cloud Cover
    ENSURE_VAR( driver_id, 'DSWRF_221_SFC_ave3h', ncswdid ) ! Surface solar rad downwards
    ENSURE_VAR( driver_id, 'DLWRF_221_SFC_ave3h', ncldwid ) ! Surface thermal rad down
    ENSURE_VAR( driver_id, 'U_GRD_221_HTGL', ncuwdid ) ! U-wind at 10 m
    ENSURE_VAR( driver_id, 'V_GRD_221_HTGL', ncvwdid ) ! V-wind at 10 m
    ENSURE_VAR( driver_id, 'SPF_H_221_HTGL', ncshid ) ! humidity at 2 m
    ENSURE_VAR( driver_id, 'PRES_221_HTGL', ncsfpid ) ! Log Surface Pressure
    ENSURE_VAR( driver_id, 'A_PCP_221_SFC_acc3h', nclspid ) ! Large Scale Precipitation
    ENSURE_VAR( driver_id, 'ACPCP_221_SFC_acc3h', nccvpid ) ! Convective Precipitation
    ENSURE_VAR( driver_id, 'SNO_D_221_SFC', ncsflid ) ! Snow Fall

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
    
! ncp, convert snowfall from meters to millimeters
    xx = 1000. * xx
! ncp, convert precip from kg/m2 to kg/m2/s
    lsp = lsp / (3600. * 3.)
    cvp = cvp / (3600. * 3.)

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
        sib(i)%prog%sw_dwn2 = sib(i)%prog%sw_dwn2
        sib(i)%prog%dlwbot2 = sib(i)%prog%dlwbot2
        if ( sib(i)%prog%sw_dwn2 < 0 ) sib(i)%prog%sw_dwn2 = 0.0
        if ( sib(i)%prog%dlwbot2 < 0 ) sib(i)%prog%dlwbot2 = 0.0
        
        ! convert total cloud cover to fraction
        sib(i)%prog%tcc2 = sib(i)%prog%tcc2 * 0.01

        ! 10 m wind
        sib(i)%prog%spdm2=SQRT(uwd(subset(i))*uwd(subset(i))+vwd(subset(i))*vwd(subset(i)))

        ! add snowfall to large scale precip and let SiB decide about snow.
        sib(i)%prog%lspr2 = sib(i)%prog%lspr2!+xx(subset(i))
        ! convert to mm
        sib(i)%prog%lspr2 = (sib(i)%prog%lspr2-sib(i)%prog%cupr2)*time%driver_step
        sib(i)%prog%cupr2 = sib(i)%prog%cupr2*time%driver_step
        ! make sure precip > 0
        if ( sib(i)%prog%lspr2 < 0.0 ) sib(i)%prog%lspr2 = 0.0
        if ( sib(i)%prog%cupr2 < 0.0 ) sib(i)%prog%cupr2 = 0.0

        !kdcorbin, 07/09 - make sure that humidity at 2m is greater than zero
        if(sib(i)%prog%sh2<=0.) sib(i)%prog%sh2=1.e-4

        !kdcorbin, 07/09 - remove any -999's in temperature
        if(sib(i)%prog%tm2 < 0.0) then
           print *,'encountered 2 meter temp = ',sib(i)%prog%tm2,'K'
           print*,'...changing to 285K'
           sib(i)%prog%tm2 = 285.0
        endif

        ! Conversion Pa -> hPa (pascals to millibars)
        sib(i)%prog%ps2 = sib(i)%prog%ps2 * 0.01
    enddo

!    print*,subname,'New driver data read ',ihour,iday,imon,iyear
!    print*,'------------------------------------------------------'
!    print*,'Extrema of new input data'
!    print*, minval(sib%prog%tm2      ),maxval(sib%prog%tm2  ),' Temperature'
!    print*, minval(sib%prog%tcc2     ),maxval(sib%prog%tcc2),' Total cloudiness'
!    print*, minval(sib%prog%sh2    ),maxval(sib%prog%sh2 ),' dew point'
!    print*, minval(sib%prog%spdm2    ),maxval(sib%prog%spdm2),' Surface wind'
!    print*, minval(sib%prog%ps2      ),maxval(sib%prog%ps2 ),' Pressure'
!    print*, minval(sib%prog%dlwbot2  ),maxval(sib%prog%dlwbot2),  &
!        ' Long wave down'
!    print*, minval(sib%prog%lspr2    ),maxval(sib%prog%lspr2),' Large sc precip'
!    print*, minval(sib%prog%cupr2    ),maxval(sib%prog%cupr2),' Convective '
!    print*, minval(sib%prog%sw_dwn2),maxval(sib%prog%sw_dwn2),  &
!        ' Short wave down'
!    print*,'-----------------------------------------------------'


end subroutine sibdrv_read_narr

subroutine sibdrv_read_geos4( sib, time )

!****--------------------------------------------------------------------
!    This subroutines reads the forcing surface meteorological data 
!    for the next time step.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.

!    precip, snowfall and radiation data are aggregated data over 3 hours
!    They have been stored such that the data at time x represent the
!    aggregation over the preceeding 3 hours.
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
!  Kevin Schaefer removed conversion from millibars to pascals (8/16/04)
!  Kevin Schaefer changed from tdew1/2 to sh1/2 because thts whtas read in (8/17/04)
!  Kevin Schaefer added calls to error handling routine (11/12/04)
!  Kevin Schaefer deleted retreival of undefimed time variable (11/12/04)
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
#ifdef PGF
use netcdf
use typeSizes
#endif 
use sibtype
use timetype

type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time

real(kind=dbl_kind) :: pid180 
real(kind=dbl_kind) :: cosz(nsib)

integer(kind=int_kind) :: i, iyear, imon, iday, idoy, ihour, imin

integer(kind=int_kind) ::  n

character*100 filename
character*7 gchar
integer(kind=int_kind) :: nctimeid,ncyid,ncmid,nctdid,ncdoyid,nchid
integer(kind=int_kind), dimension(2) :: mstart,mcount

integer(kind=int_kind) :: nct2mid ! Temperature at 2 m
integer(kind=int_kind) :: ncswdid ! Surface solar radiation downwards           
integer(kind=int_kind) :: albedoid! Albedo
integer(kind=int_kind) :: ncldwid ! Surface thermal radiation downwards
integer(kind=int_kind) :: ncuwdid ! U-wind at 10 m
integer(kind=int_kind) :: ncvwdid ! V-wind at 10 m
integer(kind=int_kind) :: ncshid ! humidity at 2 m
integer(kind=int_kind) :: ncsfpid ! Log Surface Pressure
integer(kind=int_kind) :: nclspid ! Large Scale Precipitation
integer(kind=int_kind) :: nccvpid ! Convective Precipitation

real(kind=real_kind), dimension(nsib) :: t2m ! Total Cloud Cover
real(kind=real_kind), dimension(nsib) :: swd ! Surface solar radiation downwards           
real(kind=real_kind), dimension(nsib) :: alb ! Total Cloud Cover
real(kind=real_kind), dimension(nsib) :: ldw ! Surface thermal radiation downwards
real(kind=real_kind), dimension(nsib) :: sh ! humidity at 2 m
real(kind=real_kind), dimension(nsib) :: sfp ! Log Surface Pressure
real(kind=real_kind), dimension(nsib) :: lsp ! Large Scale Precipitation
real(kind=real_kind), dimension(nsib) :: cvp ! Convective Precipitation

integer(kind=int_kind) :: status
real(kind=real_kind) :: xtime,xyear,xmonth,xdoy,xday,xhour
real(kind=dbl_kind), dimension(nsib) :: xx,uwd,vwd

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
        status = nf90_open( trim(filename), nf90_nowrite, driver_id )
        if(status/=nf90_noerr) call handle_err(status,'read_geos4',1)
        print *, 'drvr file switched to ',trim(filename)
    endif

    ! read new driver data
    ! read time values from driver data file
    status = nf90_inq_varid( driver_id,   'year', ncyid )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',3)
    status = nf90_inq_varid( driver_id,   'month',ncmid )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',4)
    status = nf90_inq_varid( driver_id,   'doy',  ncdoyid )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',5)
    status = nf90_inq_varid( driver_id,   'day',  nctdid )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',6)
    status = nf90_inq_varid( driver_id,   'hour', nchid )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',7)

    ! read time
    mstart(1) = time%driver_recnum
    status = nf90_get_var( driver_id, ncyid,    xyear, mstart(1:1) )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',9)
    status = nf90_get_var( driver_id, ncmid,   xmonth, mstart(1:1) )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',10)
    status = nf90_get_var( driver_id, ncdoyid,   xdoy, mstart(1:1) )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',11)
    status = nf90_get_var( driver_id, nctdid,    xday, mstart(1:1) )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',12)
    status = nf90_get_var( driver_id, nchid,    xhour, mstart(1:1) )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',13)

    ! get veriable id's
    status=nf90_inq_varid( driver_id, 't2m', nct2mid ) ! Temperature at 2 m
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',14)
    status=nf90_inq_varid( driver_id, 'radswg', ncswdid ) ! Surface solar rad downwards
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',15)
    status=nf90_inq_varid( driver_id, 'albedo', albedoid ) ! Albedo
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',16)
    status=nf90_inq_varid( driver_id, 'lwgdown', ncldwid ) ! Surface thermal rad down
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',17)
    status=nf90_inq_varid( driver_id, 'u10m', ncuwdid ) ! U-wind at 10 m
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',18)
    status=nf90_inq_varid( driver_id, 'v10m', ncvwdid ) ! V-wind at 10 m
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',19)
    status=nf90_inq_varid( driver_id, 'q2m', ncshid ) ! humidity at 2 m
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',20)
    status=nf90_inq_varid( driver_id, 'ps', ncsfpid ) ! Log Surface Pressure
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',21)
    status=nf90_inq_varid( driver_id, 'preacc', nclspid ) ! Large Scale Precipitation
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',22)
    status=nf90_inq_varid( driver_id, 'precon', nccvpid ) ! Convective Precipitation
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',23)

    ! get data
    mstart=(/1,time%driver_recnum/); mcount=(/nsib,1/)
    status = nf90_get_var( driver_id, nct2mid, t2m,     & !Temperature at 2 m
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',24)
    status = nf90_get_var( driver_id, ncswdid, swd,     & !Surface solar rad downwards
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',25)
    status = nf90_get_var( driver_id, albedoid, alb,    & !Albedo
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',26)
    status = nf90_get_var( driver_id, ncldwid, ldw,     & !Surface thermal rad downwards
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',27)
    status = nf90_get_var( driver_id, ncuwdid, uwd,     & ! U-wind at 10 m
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',28)
    status = nf90_get_var( driver_id, ncvwdid, vwd,     & ! V-wind at 10 m
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',29)
    status = nf90_get_var( driver_id, ncshid, sh,     & ! specific humidity
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',30)
    status = nf90_get_var( driver_id, ncsfpid, sfp,     & ! Surface Pressure
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',31)
    status = nf90_get_var( driver_id, nclspid, lsp,     & ! Large Scale Precipitation
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',32)
    status = nf90_get_var( driver_id, nccvpid, cvp,     & ! Convective Precipitation
         mstart,  mcount )
    if(status/=nf90_noerr) call handle_err(status,'read_geos4',33)


    do i=1,subcount
        ! pull out landpoints in subdomain
        sib(i)%prog%tm2 = t2m(subset(i))
        sib(i)%prog%sw_dwn2 = swd(subset(i))/(1.-alb(subset(i)))
        sib(i)%prog%dlwbot2 = ldw(subset(i))
        sib(i)%prog%sh2 = sh(subset(i))
        sib(i)%prog%ps2 = sfp(subset(i))
        sib(i)%prog%lspr2 = lsp(subset(i))
        sib(i)%prog%cupr2 = cvp(subset(i))
        sib(i)%prog%tcc2 = 0.0_dbl_kind
    
        ! 10 m wind
        sib(i)%prog%spdm2=SQRT(uwd(subset(i))*uwd(subset(i))+vwd(subset(i))*vwd(subset(i)))

        ! convert to mm
        sib(i)%prog%lspr2 = sib(i)%prog%lspr2/8.0_dbl_kind
        sib(i)%prog%cupr2 = sib(i)%prog%cupr2/8.0_dbl_kind
        sib(i)%prog%lspr2 = sib(i)%prog%lspr2 - sib(i)%prog%cupr2

        !convert specific humidity from g kg-1 to kg kg-1
        sib(i)%prog%sh2 = sib(i)%prog%sh2*0.001_dbl_kind

    enddo

    print*,'read driver data ',xyear,xdoy,xhour
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


end subroutine sibdrv_read_geos4

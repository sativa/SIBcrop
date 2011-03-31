subroutine create_qp2( out_path, numvars, subcount, ihr, jhr, year,  &
                       month, longitude, latitude, lonindex, latindex,  &
                       doqpsib, nameqpsib, listqpsib, qp2id, qp2varid,   &
                       qp2timeid, qp2charid, qp2startid, qp2endid,   &
                       qp2periodid, drvr_type, biome_source, soil_source,  &
                       soref_source, ndvi_source, c4_source, d13cresp_source,  &
                       rank )
use netcdf
use typeSizes

#include "nc_util.h"

! parameters
character(len=256), intent(in) :: out_path      ! directory to write file to
integer, intent(in) :: numvars                  ! max # of variables written to file
integer, intent(in) :: subcount                 ! number of landpoints
integer, intent(in) :: ihr                      ! number of longitude indices
integer, intent(in) :: jhr                      ! number of latitude indices
integer, intent(in) :: year                     ! current year (used in filename)
integer, intent(in) :: month                    ! current month ( "  " )
real, dimension(ihr), intent(in) :: longitude   ! array of longitude coordinates
real, dimension(jhr), intent(in) :: latitude    ! array of latitude coordinates
integer, dimension(subcount), intent(in) :: lonindex ! array of longitude indices
integer, dimension(subcount), intent(in) :: latindex ! array of latitude indices
logical, dimension(numvars), intent(in) :: doqpsib    ! defines which vars. to write to file
character(len=16), dimension(numvars), intent(in) :: nameqpsib    ! names of vars.
character(len=80), dimension(numvars), intent(in) :: listqpsib    ! variable descriptions
integer, intent(out) :: qp2id                   ! qp file id#
integer, dimension(numvars), intent(out):: qp2varid     ! variable id#s 
integer, intent(out) :: qp2timeid               ! time variable id#
integer, intent(out) :: qp2charid               ! char_time variable id#
integer, intent(out) :: qp2startid              ! start_period variable id#
integer, intent(out) :: qp2endid                ! end_period variable id#
integer, intent(out) :: qp2periodid             ! period_length variable id#
character(len=8), intent(in) :: drvr_type       ! type of driver data used
character(len=100), intent(in) :: biome_source
character(len=100), intent(in) :: soil_source
character(len=100), intent(in) :: soref_source
character(len=100), intent(in) :: ndvi_source
character(len=100), intent(in) :: c4_source
character(len=100), intent(in) :: d13cresp_source
integer, intent(in) :: rank

! netcdf variables
integer :: status           ! return status of netcdf functions
integer :: latid            ! latitude dimension id #
integer :: lonid            ! longitude dimension id #
integer :: timeid           ! time dimension id #
integer :: charid           ! char_len dimension id #
integer :: latitudeid       ! latitude variable id #
integer :: longitudeid      ! longitude variable id #
integer :: subcountid       ! landpoints variable id #
integer :: lonindexid      ! longitude indexing array id #
integer :: latindexid      ! latitude indexing array id #

! local variables
integer :: n                        ! index variable
character(len=256) :: filename      ! filename string
character(len=40) :: units          ! variable units
character(len=80) :: longname       ! variable description
integer :: unit_len, long_len       ! not used, returned by get_units()

    ! make sure qp2id is not tied to any open file
    status = nf90_close( qp2id )

    ! create file
    write( filename, '(a,i4.4,i2.2,a,i3.3,a)' ) trim(out_path)//'hsib_',   &
        year, month, 'p', rank, '.qp2.nc'
    CHECK( nf90_create( trim(filename), nf90_clobber, qp2id) )
    
    ! define global attributes
    call global_atts( qp2id, 'sib3', 'lat/lon', '1.0', drvr_type,  &
        biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
        d13cresp_source, rank )

    ! define dimensions
    CHECK( nf90_def_dim( qp2id, 'time', nf90_unlimited, timeid ) )
    CHECK( nf90_def_dim( qp2id, 'char_len', 10, charid ) )
    CHECK( nf90_def_dim( qp2id, 'latitude', jhr, latid ) )
    CHECK( nf90_def_dim( qp2id, 'longitude', ihr, lonid ) )
    CHECK( nf90_def_dim( qp2id, 'landpoints', subcount, subcountid ) )
    
    ! define control variables
    CHECK( nf90_def_var( qp2id, 'time', nf90_double, (/timeid/), qp2timeid ) )
    CHECK( nf90_put_att( qp2id, qp2timeid, 'quantity', 'time' ) )
    CHECK( nf90_put_att( qp2id, qp2timeid, 'units', 'days since 1-1-1' ) )
    CHECK( nf90_put_att( qp2id, qp2timeid, 'calender', 'noleap' ) )
    
    CHECK( nf90_def_var( qp2id, 'char_time', nf90_char, (/charid,timeid/), qp2charid ) )
    CHECK( nf90_put_att( qp2id, qp2charid, 'format', 'mm/dd/yyyy' ) )
    
    CHECK( nf90_def_var( qp2id, 'start_period', nf90_double, (/timeid/), qp2startid ) )
    CHECK( nf90_put_att( qp2id, qp2startid, 'long_name', 'start of averaged period' ) )
    CHECK( nf90_put_att( qp2id, qp2startid, 'units', 'days since 1-1-1' ) )
    
    CHECK( nf90_def_var( qp2id, 'end_period', nf90_double, (/timeid/), qp2endid ) )
    CHECK( nf90_put_att( qp2id, qp2endid, 'long_name', 'end of averaged period' ) )
    CHECK( nf90_put_att( qp2id, qp2endid, 'units', 'days since 1-1-1' ) )
    
    CHECK( nf90_def_var( qp2id, 'period_length', nf90_double, (/timeid/), qp2periodid ) )
    CHECK( nf90_put_att( qp2id, qp2periodid, 'long_name', 'length of averaged period' ) )
    CHECK( nf90_put_att( qp2id, qp2periodid, 'units', 'days' ) )
    
    CHECK( nf90_def_var( qp2id, 'latitude', nf90_float, (/latid/), latitudeid ) )
    CHECK( nf90_put_att( qp2id, latitudeid, 'units', 'degrees_north' ) )
    CHECK( nf90_put_att( qp2id, latitudeid, 'quantity', 'latitude' ) )
    
    CHECK( nf90_def_var( qp2id, 'longitude', nf90_float, (/lonid/), longitudeid ) )
    CHECK( nf90_put_att( qp2id, longitudeid, 'units', 'degrees_east' ) )
    CHECK( nf90_put_att( qp2id, longitudeid, 'quantity', 'longitude' ) )
    
    CHECK( nf90_def_var( qp2id, 'lonindex', nf90_int, (/subcountid/), lonindexid ) )
    CHECK( nf90_put_att( qp2id, lonindexid, 'long_name', 'Longitude array index' ) )
    CHECK( nf90_put_att( qp2id, lonindexid, 'units', 'index-integer' ) )

    CHECK( nf90_def_var( qp2id, 'latindex', nf90_int, (/subcountid/), latindexid ) )
    CHECK( nf90_put_att( qp2id, latindexid, 'long_name', 'Latitude array index' ) )
    CHECK( nf90_put_att( qp2id, latindexid, 'units', 'index-integer' ) )

    ! define data variables
    do n = 1, numvars
        if ( doqpsib(n) ) then
           status = nf90_def_var( qp2id, trim(nameqpsib(n)), nf90_float, &
                (/subcountid,timeid/), qp2varid(n) )
            CHECK( status )
            call get_units( listqpsib(n), longname, long_len, units, unit_len )
            CHECK( nf90_put_att( qp2id, qp2varid(n), 'long_name', trim(longname) ) )
            CHECK( nf90_put_att( qp2id, qp2varid(n), 'title', trim(longname) ) )
            CHECK( nf90_put_att( qp2id, qp2varid(n), 'units', trim(units) ) )
            CHECK( nf90_put_att( qp2id, qp2varid(n), 'missing_value', 1.e36 ) )
        endif
    enddo
    
    ! switch from definition mode to data mode
    CHECK( nf90_enddef( qp2id ) )
    
    ! assign values to variables not variant with time
    CHECK( nf90_put_var( qp2id, latitudeid, latitude ) )
    CHECK( nf90_put_var( qp2id, longitudeid, longitude ) )
    CHECK( nf90_put_var( qp2id, lonindexid, lonindex ) )
    CHECK( nf90_put_var( qp2id, latindexid, latindex ) )

end subroutine create_qp2


!-----------------------------------------------------------------------

subroutine write_qp2( qp2id, qp2timeid, qp2startid, qp2endid, qp2periodid,  &
                      qp2charid, numvars, subcount, qp2varid, qpsib,  &
                      doqpsib, indxqpsib, year, month, day,  &
                      seconds, end_period, period_length )
use netcdf
use typeSizes
use kinds
use sib_const_module, only: dtsib

! parameters
integer, intent(in) :: qp2id            ! file id#
integer, intent(in) :: qp2timeid        ! time variable id #
integer, intent(in) :: qp2startid       ! start_period id #
integer, intent(in) :: qp2endid         ! end_period id #
integer, intent(in) :: qp2periodid      ! period_length id #
integer, intent(in) :: qp2charid        ! char_time id #
integer, intent(in) :: numvars          ! max # of variables written to file
integer, intent(in) :: subcount         ! number of landpoints in subdomain
integer, dimension(numvars), intent(in) :: qp2varid   ! variable id #s
real(kind=dbl_kind), dimension(subcount,numvars), intent(in) :: qpsib     ! variable values
logical, dimension(numvars), intent(in) :: doqpsib    ! defines which variables to write out
integer, dimension(numvars), intent(in) :: indxqpsib  ! ??????????
integer, intent(in) :: year     ! current year (used in char_time)
integer, intent(in) :: month    ! current month (used in char_time)
integer, intent(in) :: day      ! current day (used in char_time)
integer, intent(in) :: seconds  ! current second of the year
double precision, intent(in) :: end_period      ! end of averaged period
double precision, intent(in) :: period_length   ! length of averaged period

! netcdf variables
integer :: status       ! return status of netcdf functions
integer :: dimid        ! dimension id #

! local variables
integer :: n, i                      ! index variables
integer :: step                      ! next time step in qp file
double precision :: dyear
character(len=10) :: char_time       ! mm/dd/yyyy
character(len=10) :: name
double precision :: secyear = 86400.

    ! find next time step
    CHECK( nf90_inq_dimid( qp2id, 'time', dimid ) )
    CHECK( nf90_inquire_dimension( qp2id, dimid, name, step ) )
    step = step + 1
    
    ! write out time variables
    dyear = seconds/secyear
    CHECK( nf90_put_var( qp2id, qp2timeid, dyear, (/step/) ) )
    CHECK( nf90_put_var( qp2id, qp2startid, end_period - period_length, (/step/) ) )
    CHECK( nf90_put_var( qp2id, qp2endid, end_period, (/step/) ) )
    CHECK( nf90_put_var( qp2id, qp2periodid, period_length, (/step/) ) )
    
    write(char_time, '(i2.2,a1,i2.2,a1,i4.4)') month, '/', day, '/', year
    CHECK( nf90_put_var( qp2id, qp2charid, char_time, (/1,step/), (/10,1/) ) )
    
    ! write out data variables
    do n = 1, numvars
        if ( doqpsib(n) ) then
           status = nf90_put_var( qp2id, qp2varid(n), qpsib(:,indxqpsib(n)), &
                (/1,step/), (/subcount,1/) )
            CHECK( status )
        endif
    enddo

end subroutine write_qp2

!-----------------------------------------------------------------------

subroutine create_qp3( out_path, numvars, subcount, ihr, jhr,  &
                       year, month, nsoil, longitude, latitude,  &
                       lonindex, latindex, doqp3sib, nameqp3sib,   &
                       listqp3sib, drvr_type, biome_source, soil_source,  &
                       soref_source, ndvi_source, c4_source, d13cresp_source,  &
                       qp3id, qp3varid, qp3timeid, qp3charid, qp3startid,  &
                       qp3endid, qp3periodid, rank )
use netcdf
use typeSizes

! parameters
character(len=256), intent(in) :: out_path
integer, intent(in) :: numvars          ! # of variables to write out
integer, intent(in) :: subcount         ! # of landpoints
integer, intent(in) :: ihr              ! # of longitude indices
integer, intent(in) :: jhr              ! # of latitude indices
integer, intent(in) :: year             ! current year of simulation
integer, intent(in) :: month            ! current month
integer, intent(in) :: nsoil            ! # of soil layers
real, dimension(ihr), intent(in) :: longitude   ! array of longitude coordinates
real, dimension(jhr), intent(in) :: latitude    ! array of latitude coordinates
integer, dimension(subcount), intent(in) :: lonindex
integer, dimension(subcount), intent(in) :: latindex
logical, dimension(numvars), intent(in) :: doqp3sib  ! defines which variables are written out
character(len=16), dimension(numvars), intent(in) :: nameqp3sib  ! variable names
character(len=80), dimension(numvars), intent(in) :: listqp3sib  ! variable descriptions
character(len=8), intent(in) :: drvr_type   ! driver data type used for simulation
character(len=100), intent(in) :: biome_source
character(len=100), intent(in) :: soil_source
character(len=100), intent(in) :: soref_source
character(len=100), intent(in) :: ndvi_source
character(len=100), intent(in) :: c4_source
character(len=100), intent(in) :: d13cresp_source
integer, intent(out) :: qp3id           ! file id #
integer, dimension(numvars), intent(out) :: qp3varid     ! variable id #s
integer, intent(out) :: qp3timeid       ! time variable id #
integer, intent(out) :: qp3charid       ! char_time id #
integer, intent(out) :: qp3startid      ! start_period id #
integer, intent(out) :: qp3endid        ! end_period id #
integer, intent(out) :: qp3periodid     ! period_length id #
integer, intent(in) :: rank

! netcdf variables
integer :: status           ! return status of netcdf functions
integer :: latid            ! latitude dimension id #
integer :: lonid            ! longitude dimension id #
integer :: timeid           ! time dimension id #
integer :: charid           ! char_len dimension id #
integer :: levelid          ! level dimension id #
integer :: latitudeid       ! latitude variable id #
integer :: longitudeid      ! longitude variable id #
integer :: levid            ! level variable id #
integer :: landpointsid     ! landpoints dimension id #
integer :: lonindexid      ! lonindex variable id #
integer :: latindexid      ! latindex variable id #


! local variables
integer :: n                        ! index variable
character(len=40) :: units          ! variable units
character(len=80) :: longname       ! variable description
integer :: unit_len, long_len       ! not used, used by get_units()
real, dimension(nsoil) :: levels    ! numbering of levels
character(len=256) :: filename      ! file name

    ! make sure qp3id is not tied to any open file
    status = nf90_close( qp3id )
    
    ! create file
    write( filename, '(a,i4.4,i2.2,a,i3.3,a)' ) trim(out_path)//'hsib_',   &
        year, month, 'p', rank, '.qp3.nc'
    CHECK( nf90_create( trim(filename), nf90_clobber, qp3id) )
    
    ! define global attributes
    call global_atts( qp3id, 'sib3', 'lat/lon', '1.0', drvr_type,  &
        biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
        d13cresp_source, rank )
    
    ! define dimensions
    CHECK( nf90_def_dim( qp3id, 'time', nf90_unlimited, timeid ) )
    CHECK( nf90_def_dim( qp3id, 'char_len', 10, charid ) )
    CHECK( nf90_def_dim( qp3id, 'latitude', jhr, latid ) )
    CHECK( nf90_def_dim( qp3id, 'longitude', ihr, lonid ) )
    CHECK( nf90_def_dim( qp3id, 'level', nsoil, levelid ) )
    CHECK( nf90_def_dim( qp3id, 'landpoints', subcount, landpointsid ) )
    
    ! define control variables
    CHECK( nf90_def_var( qp3id, 'time', nf90_double, (/timeid/), qp3timeid ) )
    CHECK( nf90_put_att( qp3id, qp3timeid, 'quantity', 'time' ) )
    CHECK( nf90_put_att( qp3id, qp3timeid, 'units', 'days since 1-1-1' ) )
    CHECK( nf90_put_att( qp3id, qp3timeid, 'calender', 'noleap' ) )
    
    CHECK( nf90_def_var( qp3id, 'char_time', nf90_char, (/charid,timeid/), qp3charid ) )
    CHECK( nf90_put_att( qp3id, qp3charid, 'format', 'mm/dd/yyyy' ) )
    
    CHECK( nf90_def_var( qp3id, 'start_period', nf90_double, (/timeid/), qp3startid ) )
    CHECK( nf90_put_att( qp3id, qp3startid, 'long_name', 'start of averaged period' ) )
    CHECK( nf90_put_att( qp3id, qp3startid, 'units', 'days since 1-1-1' ) )
    
    CHECK( nf90_def_var( qp3id, 'end_period', nf90_double, (/timeid/), qp3endid ) )
    CHECK( nf90_put_att( qp3id, qp3endid, 'long_name', 'end of averaged period' ) )
    CHECK( nf90_put_att( qp3id, qp3endid, 'units', 'days since 1-1-1' ) )
    
    CHECK( nf90_def_var( qp3id, 'period_length', nf90_double, (/timeid/), qp3periodid ) )
    CHECK( nf90_put_att( qp3id, qp3periodid, 'long_name', 'length of averaged period' ) )
    CHECK( nf90_put_att( qp3id, qp3periodid, 'units', 'days' ) )
    
    CHECK( nf90_def_var( qp3id, 'latitude', nf90_float, (/latid/), latitudeid ) )
    CHECK( nf90_put_att( qp3id, latitudeid, 'units', 'degrees_north' ) )
    CHECK( nf90_put_att( qp3id, latitudeid, 'quantity', 'latitude' ) )
    
    CHECK( nf90_def_var( qp3id, 'longitude', nf90_float, (/lonid/), longitudeid ) )
    CHECK( nf90_put_att( qp3id, longitudeid, 'units', 'degrees_east' ) )
    CHECK( nf90_put_att( qp3id, longitudeid, 'quantity', 'longitude' ) )
    
    CHECK( nf90_def_var( qp3id, 'level', nf90_float, (/levelid/), levid ) )

    CHECK( nf90_def_var( qp3id, 'lonindex', nf90_int, (/landpointsid/), lonindexid ) )
    CHECK( nf90_put_att( qp3id, lonindexid, 'long_name', 'longitude index array' ) )
    CHECK( nf90_put_att( qp3id, lonindexid, 'units', 'index-integer' ) )

    CHECK( nf90_def_var( qp3id, 'latindex', nf90_int, (/landpointsid/), latindexid ) )
    CHECK( nf90_put_att( qp3id, latindexid, 'long_name', 'latitude index array' ) )
    CHECK( nf90_put_att( qp3id, latindexid, 'units', 'index-integer' ) )
    
    ! define data variables
    do n = 1, numvars
        if ( doqp3sib(n) ) then
           status = nf90_def_var( qp3id, trim(nameqp3sib(n)), nf90_float, &
                (/landpointsid,levelid,timeid/), qp3varid(n) )
            CHECK( status )
            call get_units( listqp3sib(n), longname, long_len, units, unit_len )
            CHECK( nf90_put_att( qp3id, qp3varid(n), 'long_name', trim(longname) ) )
            CHECK( nf90_put_att( qp3id, qp3varid(n), 'title', trim(longname) ) )
            CHECK( nf90_put_att( qp3id, qp3varid(n), 'units', trim(units) ) )
            CHECK( nf90_put_att( qp3id, qp3varid(n), 'missing_value', 1.e36 ) )
        endif
    enddo
    
    ! switch from definition mode to data mode
    CHECK( nf90_enddef( qp3id ) )
    
    ! assign values to variables not variant with time
    CHECK( nf90_put_var( qp3id, latitudeid, latitude ) )
    CHECK( nf90_put_var( qp3id, longitudeid, longitude ) )
    CHECK( nf90_put_var( qp3id, latindexid, latindex ) )
    CHECK( nf90_put_var( qp3id, lonindexid, lonindex ) )
    
    do n = 1, nsoil
        levels(n) = n
    enddo
    CHECK( nf90_put_var( qp3id, levid, levels ) )

end subroutine create_qp3

!-----------------------------------------------------------------------

subroutine write_qp3( qp3id, qp3timeid, qp3startid, qp3endid, qp3periodid,  &
                      qp3charid, numvars, subcount, nsoil, qp3varid,  &
                      qp3sib, doqp3sib, indxqp3sib, year, month,  &
                      day, seconds, end_period, period_length )
!
! Modifications:
!  Kevin Schaefer removed _date print statement (8/18/04)
use netcdf
use typeSizes
use kinds

! parameters
integer, intent(in) :: qp3id            ! file id #
integer, intent(in) :: qp3timeid        ! time variable id #
integer, intent(in) :: qp3startid       ! start_period id #
integer, intent(in) :: qp3endid         ! end_period id #
integer, intent(in) :: qp3periodid      ! period_length id #
integer, intent(in) :: qp3charid        ! char_time id #
integer, intent(in) :: numvars          ! # of possible vars. written to file
integer, intent(in) :: subcount         ! # landpoints in subdomain
integer, intent(in) :: nsoil            ! # soil layers
integer, dimension(numvars), intent(in) :: qp3varid    ! variable id #s
real(kind=dbl_kind), dimension(subcount,nsoil,numvars), intent(in) :: qp3sib ! variable values
logical, dimension(numvars), intent(in) :: doqp3sib      ! defines which variables to write out
integer, dimension(numvars), intent(in) :: indxqp3sib    ! ???????????
integer, intent(in) :: year     ! current year, used for char_time
integer, intent(in) :: month    ! current month, used for char_time
integer, intent(in) :: day      ! current day, used for char_time
integer, intent(in) :: seconds  ! current second of the year
double precision, intent(in) :: end_period      ! end of averaged period
double precision, intent(in) :: period_length   ! length of averaged period

! netcdf variables
integer :: status   ! return status of netcdf functions
integer :: dimid    ! dimension id #

! local variables
integer :: n, i, l      ! index variables
integer :: step         ! next time step in qp file
double precision :: dyear
character(len=10) :: char_time      ! mm/dd/yyyy
double precision :: secyear = 86400.
character(len=10) :: name
    ! find next time step
    CHECK( nf90_inq_dimid( qp3id, 'time', dimid ) )
    CHECK( nf90_inquire_dimension( qp3id, dimid, name, step ) )
    step = step + 1
    
    ! write out time variables
    dyear = seconds/secyear
    CHECK( nf90_put_var( qp3id, qp3timeid, dyear, (/step/) ) )
    CHECK( nf90_put_var( qp3id, qp3startid, end_period - period_length, (/step/) ) )
    CHECK( nf90_put_var( qp3id, qp3endid, end_period, (/step/) ) )
    CHECK( nf90_put_var( qp3id, qp3periodid, period_length, (/step/) ) )

    write(char_time, '(i2.2,a1,i2.2,a1,i4.4)') month, '/', day, '/', year
    CHECK( nf90_put_var( qp3id, qp3charid, char_time, (/1,step/), (/10,1/) ) )
    
    ! write out data variables
    do n = 1, numvars
        if ( doqp3sib(n) ) then
           status = nf90_put_var( qp3id, qp3varid(n),  &
                qp3sib(:,:,indxqp3sib(n)), (/1,1,step/),  &
                (/subcount,nsoil,1/) )
            CHECK( status )
        endif
    enddo

end subroutine write_qp3

!-----------------------------------------------------------------------

subroutine global_atts (fileID, runname, grid, version, driver,  &
    biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
    d13cresp_source, rank )

!-----------------------------------------------------------------------
! Purpose:
!   sets global attributes 
!
! Scope: 
!   Module variables used:
!
! Bugs:
!   netcdf file (fileID) must be in define mode.
!
!-----------------------------------------------------------------------
use netcdf
use typeSizes

! input parameters
integer, intent(in) :: fileID
character(len=*), intent(in) :: runname
character(len=*), intent(in) :: grid
character(len=*), intent(in) :: version
character(len=*), intent(in) :: driver
character(len=100), intent(in) :: biome_source
character(len=100), intent(in) :: soil_source
character(len=100), intent(in) :: soref_source
character(len=100), intent(in) :: ndvi_source
character(len=100), intent(in) :: c4_source
character(len=100), intent(in) :: d13cresp_source
integer, intent(in) :: rank

! local variables
integer :: status
character(len=30) :: current_time
character(len=8) :: t_date
character(len=10) :: t_time
character(len=5) :: zone
character(len=4) :: c_rank
integer, dimension(8) :: values

    call date_and_time(t_date, t_time, zone, values)

    current_time = t_date(5:6) // "/" // t_date(7:8) // "/" // t_date(1:4)   &
        //" at "// t_time(1:2) // ":" //t_time(3:4) // " "      &
        // zone // " GMT "

    write( c_rank, '(i4.4)' ) rank

    !   add standard global attributes
    CHECK( nf90_put_att ( fileID, nf90_global, 'calendar', 'noleap' ) )
    CHECK( nf90_put_att ( fileID, nf90_global, 'institution', 'Colorado State University' ) )
    CHECK( nf90_put_att ( fileID, nf90_global, 'history', 'Created: '//current_time ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'run', runname ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'rank', c_rank )  )
    CHECK( nf90_put_att( fileID, nf90_global, 'grid', grid ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'version', version ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'Driver_Data', driver ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'biome_source', trim(biome_source) ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'soil_source', trim(soil_source) ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'soref_source', trim(soref_source) ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'ndvi_source', trim(ndvi_source) ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'c4_source', trim(c4_source) ) )
    CHECK( nf90_put_att( fileID, nf90_global, 'd13cresp_source', trim(d13cresp_source) ) )

end subroutine global_atts


!-----------------------------------------------------------------------

subroutine get_units(description, longname, long_len, units, unit_len)

!-----------------------------------------------------------------------
!  Purpose:
!   extracts a string enclosed within parentheses - used for 
!   units which are contained in a general description string.
!   returns the units string (units) and its length (unit_len),
!   the description string with the units removed (longname),
!   and its length (long_len).
!   note: embedded parentheses are ok
!
!  Variables:
!   Module parameters used:
!   MAXUNITS, MAXLONGN
!   
!  Bugs:
!   1) if the rightmost parenthesis is unmatched, units will be 
!      set to " " (one space) - this is to be interpreted as "none"
!   2) if a "(" is unmatched, it will be part of the returned
!      longname string
!       3) strings of only units (i.e. entire string enclosed in 
!      parentheses) do not work.
!-----------------------------------------------------------------------

character(len=*), intent(in) :: description
character(len=*), intent(out) :: units
character(len=*), intent(out) :: longname
integer, intent(out) :: unit_len, long_len

integer :: n, start_paren, end_paren, paren_count

    paren_count = 0
    start_paren = len_trim(description)
    end_paren = len_trim(description)

    do n = len(description), 1, -1
        if (description(n:n)==")") then
            if (paren_count == 0) then
                end_paren = n
            endif
            paren_count = paren_count + 1
        else if (description(n:n) == "(") then
            paren_count = paren_count - 1
            if (paren_count == 0) then
                start_paren = n
                exit
            endif
        end if
    end do

    !   in case of confusion, clear units and return unaltered description
    !   note: start_paren > end_paren should not be possible, but just in case...
    !         start_paren = end_paren occurs when there are no units
    !         start_paren = end_paren-1 occurs when units are "()"
    !   FIXME: n==1 is too limiting - what if I wanted only units?
    if (n == 1 .or. start_paren >= end_paren) then   ! no units
        units = " "
        unit_len = 1
        longname = trim(description)
        long_len = len_trim(longname)
    else if (start_paren == (end_paren-1)) then      ! "()" case
        units = " "
        unit_len = 1
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
        long_len = len_trim(longname)
    else                                             ! normal units
        units = description(start_paren+1:end_paren-1)
        unit_len = len_trim(units)
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
        long_len = len_trim(longname)
    end if

end subroutine get_units

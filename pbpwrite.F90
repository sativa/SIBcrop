
subroutine create_pbp( npoints, year, month, numvars, totnumvars,  &
                       latitude, longitude, dopbpsib, namepbpsib, listpbpsib,  &
                       indxpbpsib, drvr_type, biome_source, soil_source,  &
                       soref_source, ndvi_source, c4_source, d13cresp_source,  &
                       outpath,pbptimeid, pbpcharid, pbpvarid,  &
                        rank1 )
#ifdef PGF
use netcdf 
use typeSizes
#endif


! parameters
integer, intent(in) :: npoints
integer, intent(in) :: year
integer, intent(in) :: month
integer, intent(in) :: numvars
integer, intent(in) :: totnumvars
real, dimension(npoints), intent(in) :: latitude
real, dimension(npoints), intent(in) :: longitude
logical, dimension(totnumvars), intent(in) :: dopbpsib
character(len=16), dimension(totnumvars), intent(in) :: namepbpsib
character(len=80), dimension(totnumvars), intent(in) :: listpbpsib
integer, dimension(totnumvars), intent(in) :: indxpbpsib
character(len=8), intent(in) :: drvr_type
character(len=100), intent(in) :: biome_source
character(len=100), intent(in) :: soil_source
character(len=100), intent(in) :: soref_source
character(len=100), intent(in) :: ndvi_source
character(len=100), intent(in) :: c4_source
character(len=100), intent(in) :: d13cresp_source
character(len=256), intent(in) :: outpath

integer, intent(out) :: pbptimeid
integer, intent(out) :: pbpcharid
integer, dimension(numvars), intent(out) :: pbpvarid
integer, intent(in) :: rank1

! netcdf variables
integer :: status
integer :: tid
integer :: clid
integer :: npid
integer :: npointsid
integer :: latid
integer :: lonid
integer :: pbpid
! local variables
integer :: n
character(len=40) :: units          ! variable units
character(len=80) :: longname       ! variable description
integer :: unit_len, long_len       ! not used, returned by get_units()
integer, dimension(npoints) :: npoints_array
character(len =256) ::filename

    ! create file name
    write( filename, '(a,i4.4,i2.2,a,i3.3,a)' ) trim(outpath)//'psib_',  &
        year, month, 'p', rank1, '.pbp1.nc'

    ! create file and define dimensions
    status = nf90_create( trim(filename), nf90_clobber, pbpid )
    status = nf90_def_dim( pbpid, 'time', nf90_unlimited, tid )
    status = nf90_def_dim( pbpid, 'char_len', 10, clid )
    status = nf90_def_dim( pbpid, 'npoints', npoints, npointsid )
    
    ! define global atts
    call global_atts( pbpid, 'sib3', 'lat/lon', '1.0', drvr_type,  &
        biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
        d13cresp_source, rank1 )
    
    ! define variables
    status = nf90_def_var( pbpid, 'time', nf90_double, (/tid/), pbptimeid )
    status = nf90_put_att( pbpid, pbptimeid, 'quantity', 'time' )
    status = nf90_put_att( pbpid, pbptimeid, 'units', 'days since 1-1-1' )
    status = nf90_put_att( pbpid, pbptimeid, 'calendar', 'noleap' )
    
    status = nf90_def_var( pbpid, 'char_time', nf90_char, (/clid,tid/), pbpcharid )
    status = nf90_put_att( pbpid, pbpcharid, 'format', 'mm/dd/yyyy' )
    
    status = nf90_def_var( pbpid, 'npoints', nf90_int, (/npointsid/), npid )

    status = nf90_def_var( pbpid, 'latitude', nf90_float, (/npointsid/), latid )
    status = nf90_put_att( pbpid, latid, 'units', 'degrees_north' )
    status = nf90_put_att( pbpid, latid, 'quantity', 'latitude' )
    
    status = nf90_def_var( pbpid, 'longitude', nf90_float, (/npointsid/), lonid )
    status = nf90_put_att( pbpid, lonid, 'units', 'degrees_east' )
    status = nf90_put_att( pbpid, lonid, 'quantity', 'longitude' )

    do n = 1, totnumvars
        if ( dopbpsib(n) ) then
            status = nf90_def_var( pbpid, trim(namepbpsib(n)), nf90_float,  &
                (/npointsid,tid/), pbpvarid(indxpbpsib(n)) )
            call get_units( listpbpsib(n), longname, long_len, units, unit_len )
            status = nf90_put_att( pbpid, pbpvarid(indxpbpsib(n)),  &
                'long_name', trim(longname) )
            status = nf90_put_att( pbpid, pbpvarid(indxpbpsib(n)),  &
                'title', trim(longname) )
            status = nf90_put_att( pbpid, pbpvarid(indxpbpsib(n)),  &
                'units', trim(units) )
            status = nf90_put_att( pbpid, pbpvarid(indxpbpsib(n)),  &
                'missing_value', 1.e36 )
        endif
    enddo

    ! switch from definition mode to data mode
    status = nf90_enddef( pbpid )

    ! assign values to variables not variant with time
    status = nf90_put_var( pbpid, latid, latitude )
    status = nf90_put_var( pbpid, lonid, longitude )
    
    do n = 1, npoints
        npoints_array(n) = n
    enddo
    status = nf90_put_var( pbpid, npid, npoints_array(:) )
    
    status = nf90_close( pbpid )
    
end subroutine create_pbp


!-----------------------------------------------------------------------

subroutine write_pbp( npoints, year, month, day, seconds,  &
                      numvars, pbp, pbptid, pbpcid,  &
                      pbpvid,outpath, rank1 )
#ifdef PGF
use netcdf
use typeSizes
#endif
use kinds

! parameters
integer, intent(in) :: npoints
integer, intent(in) :: year
integer, intent(in) :: month
integer, intent(in) :: day
integer, intent(in) :: seconds
integer, intent(in) :: numvars
real(kind=dbl_kind), dimension(numvars+1,npoints), intent(in) :: pbp
integer, intent(in) :: pbptid
integer, intent(in) :: pbpcid
integer, intent(in), dimension(numvars) :: pbpvid
character(len=256), intent(in) :: outpath
integer, intent(in) :: rank1

! netcdf variables
integer :: status
integer :: dimid

! local variables
integer :: i, n
integer :: step
double precision :: dyear
character(len=10) :: char_time
character(len=256) :: filename
character(len=10) :: name
double precision :: secyear = 86400.
integer :: pbpid

        !open file
        write( filename, '(a,i4.4,i2.2,a,i3.3,a)' ) trim(outpath)//'psib_',  &
           year, month, 'p', rank1, '.pbp1.nc'
        status = nf90_open( trim(filename), nf90_write, pbpid )
 

    ! find next time step
    status = nf90_inq_dimid( pbpid, 'time', dimid )
    status = nf90_inquire_dimension( pbpid, dimid, name,step )
    step = step + 1
    
    ! write out time variables
    dyear = seconds/secyear
    status = nf90_put_var( pbpid, pbptid, dyear, (/step/) )

    write( char_time, '(i2.2,a1,i2.2,a1,i4.4)' ) month, '/', day, '/', year
    status = nf90_put_var( pbpid, pbpcid, char_time, (/1,step/), (/10,1/) )
  
    ! write out data variables
    do i = 1, numvars
        status = nf90_put_var( pbpid, pbpvid(i), pbp(i,:),  &
            (/1,step/), (/npoints,1/) )
    enddo
    

   status = nf90_close( pbpid )

end subroutine write_pbp


!-----------------------------------------------------------------------

subroutine create_pbp2( npoints, levels, year, month, numvars,  &
                        totnumvars, latitude, longitude, dopbp2sib,  &
                        namepbp2sib, listpbp2sib, indxpbp2sib, drvr_type,  &
                        biome_source, soil_source, soref_source, ndvi_source,  &
                        c4_source, d13cresp_source, out_path, &
                        pbp2timeid, pbp2charid, pbp2varid,rank )
#ifdef PGF
use netcdf
use typeSizes
#endif

! parameters
integer, intent(in) :: npoints
integer, intent(in) :: levels
integer, intent(in) :: year
integer, intent(in) :: month
integer, intent(in) :: numvars
integer, intent(in) :: totnumvars
real, dimension(npoints), intent(in) :: latitude
real, dimension(npoints), intent(in) :: longitude
logical, dimension(totnumvars), intent(in) :: dopbp2sib
character(len=16), dimension(totnumvars), intent(in) :: namepbp2sib
character(len=80), dimension(totnumvars), intent(in) :: listpbp2sib
integer, dimension(totnumvars), intent(in) :: indxpbp2sib
character(len=8), intent(in) :: drvr_type
character(len=100), intent(in) :: biome_source
character(len=100), intent(in) :: soil_source
character(len=100), intent(in) :: soref_source
character(len=100), intent(in) :: ndvi_source
character(len=100), intent(in) :: c4_source
character(len=100), intent(in) :: d13cresp_source
character(len=256), intent(in) :: out_path
integer, intent(out) :: pbp2timeid
integer, intent(out) :: pbp2charid
integer, dimension(numvars), intent(out) :: pbp2varid
integer, intent(in) :: rank

! netcdf variables
integer :: status
integer :: tid
integer :: charid
integer :: npointsid
integer :: levelid
integer :: npid
integer :: levid
integer :: latid
integer :: lonid
character(len=256) :: filename
integer :: pbp2id

! local variables
integer :: x
integer, dimension(npoints) :: npoints_array
integer, dimension(levels) :: levels_array
character(len=40) :: units          ! variable units
character(len=80) :: longname       ! variable description
integer :: unit_len, long_len       ! not used, returned by get_units()

    ! create file name
    write( filename, '(a,i4.4,i2.2,a,i3.3,a)' ) trim(out_path)//'psib_',  &
        year, month, 'p', rank, '.pbp2.nc'
    
    ! make sure pbp2id is not tied to any open file
    status = nf90_close( pbp2id )
    
    ! create file and define dimensions
    status = nf90_create( trim(filename), nf90_clobber, pbp2id )
    status = nf90_def_dim( pbp2id, 'time', nf90_unlimited, tid )
    status = nf90_def_dim( pbp2id, 'char_len', 10, charid )
    status = nf90_def_dim( pbp2id, 'npoints', npoints, npointsid )
    status = nf90_def_dim( pbp2id, 'level', levels, levelid )
    
    ! define global atts
    call global_atts( pbp2id, 'sib3', 'lat/lon', '1.0', drvr_type,  &
        biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
        d13cresp_source, rank )

    ! define variables
    status = nf90_def_var( pbp2id, 'time', nf90_double, (/tid/), pbp2timeid )
    status = nf90_put_att( pbp2id, pbp2timeid, 'quantity', 'time' )
    status = nf90_put_att( pbp2id, pbp2timeid, 'units', 'days since 1-1-1' )
    status = nf90_put_att( pbp2id, pbp2timeid, 'calendar', 'noleap' )
    
    status = nf90_def_var( pbp2id, 'char_time', nf90_char, (/charid,tid/), pbp2charid )
    status = nf90_put_att( pbp2id, pbp2charid, 'format', 'mm/dd/yyyy' )
    
    status = nf90_def_var( pbp2id, 'npoints', nf90_int, (/npointsid/), npid )
    
    status = nf90_def_var( pbp2id, 'latitude', nf90_float, (/npointsid/), latid )
    status = nf90_put_att( pbp2id, latid, 'units', 'degrees_north' )
    status = nf90_put_att( pbp2id, latid, 'quantity', 'latitude' )
    
    status = nf90_def_var( pbp2id, 'longitude', nf90_float, (/npointsid/), lonid )
    status = nf90_put_att( pbp2id, lonid, 'units', 'degrees_east' )
    status = nf90_put_att( pbp2id, lonid, 'quantity', 'longitude' )
    
    status = nf90_def_var( pbp2id, 'level', nf90_int, (/levelid/), levid )

    do x = 1, numvars
        if ( dopbp2sib(x) ) then
            status = nf90_def_var( pbp2id, trim(namepbp2sib(x)),  nf90_float,  &
                (/npointsid,levelid,tid/), pbp2varid(indxpbp2sib(x)) )
            call get_units( listpbp2sib(x), longname, long_len, units, unit_len )
            status = nf90_put_att( pbp2id, pbp2varid(indxpbp2sib(x)),  &
                'long_name', trim(longname) )
            status = nf90_put_att( pbp2id, pbp2varid(indxpbp2sib(x)),  &
                'title', trim(longname) )
            status = nf90_put_att( pbp2id, pbp2varid(indxpbp2sib(x)),  &
                'units', trim(units) )
            status = nf90_put_att( pbp2id, pbp2varid(indxpbp2sib(x)),  &
                'missing_value', 1.e36 )
        endif
    enddo

    ! switch from define mode to data mode
    status = nf90_enddef( pbp2id )
    
    ! assign values to variables not variant with time
    status = nf90_put_var( pbp2id, latid, latitude )
    status = nf90_put_var( pbp2id, lonid, longitude )
    
    do x = 1, npoints
        npoints_array(x) = x
    enddo
    status = nf90_put_var( pbp2id, npid, npoints_array(:) )
    
    do x = 1, levels
        levels_array(x) = x
    enddo
    status = nf90_put_var( pbp2id, levid, levels_array(:) )

	status = nf90_close( pbp2id )
    
    
end subroutine create_pbp2


!-----------------------------------------------------------------------

subroutine write_pbp2( npoints, levels, year, month, day,  &
                       seconds, numvars, pbp2, pbp2timeid,  &
                       pbp2charid, pbp2varid,out_path, rank )
#ifdef PGF
use netcdf
use typeSizes
#endif
use kinds


! parameters
integer, intent(in) :: npoints
integer, intent(in) :: levels
integer, intent(in) :: year
integer, intent(in) :: month
integer, intent(in) :: day
integer, intent(in) :: seconds
integer, intent(in) :: numvars
real(kind=dbl_kind), dimension(levels,numvars+1,npoints), intent(in) :: pbp2
integer, intent(in) :: pbp2timeid
integer, intent(in) :: pbp2charid
integer, dimension(numvars), intent(in) :: pbp2varid
character(len=256), intent(in) :: out_path
integer, intent(in) :: rank

! netcdf variables
integer :: status
integer :: dimid


! local variables
integer :: i
integer :: step
double precision :: dyear
character(len=10) :: char_time
character(len=256) :: filename
double precision :: secyear = 86400.
character(len=10) :: name
integer :: pbp2id

        ! create file name
        write( filename, '(a,i4.4,i2.2,a,i3.3,a)' ) trim(out_path)//'psib_',  &
            year, month, 'p', rank, '.pbp2.nc'
        status = nf90_open( trim(filename), nf90_write, pbp2id )


    ! find next time step
    status = nf90_inq_dimid( pbp2id, 'time', dimid )
    status = nf90_inquire_dimension( pbp2id, dimid, name,step )
    step = step + 1
    
    ! write out time variables
    dyear = seconds/secyear
    status = nf90_put_var( pbp2id, pbp2timeid, dyear, (/step/) )

    write( char_time, '(i2.2,a1,i2.2,a1,i4.4)' ) month, '/', day, '/', year
    status = nf90_put_var( pbp2id, pbp2charid, char_time, (/1,step/), (/10,1/) )
    
    ! write out data variables
    do i = 1, numvars
        status = nf90_put_var( pbp2id, pbp2varid(i), pbp2(:,i,:),  &
            (/1,1,step/), (/npoints,levels,1/) )
    enddo
 status = nf90_close( pbp2id )


end subroutine write_pbp2

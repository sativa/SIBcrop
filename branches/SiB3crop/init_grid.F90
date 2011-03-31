subroutine init_grid( rank, nchunks )
!-------------------------------------------------------------
! reads in sibdrv control variables and inputs
! sets up grid
!
! modifications:
! kdcorbin, 03/11 - removed grid_path file

use netcdf
use typeSizes

use kinds
use sib_const_module
use sib_io_module

#include "nc_util.h"

! parameters
integer(kind=int_kind), intent(in) :: rank
integer(kind=int_kind), intent(in) :: nchunks

! netcdf variables
integer(kind=int_kind) :: status    ! return value of netcdf functions
integer(kind=int_kind) :: ncid      ! netcdf file id#
integer(kind=int_kind) :: varid     ! netcdf variable id#
integer(kind=int_kind) :: dimid     ! netcdf dimension id#
integer(kind=int_kind) :: dimlen    ! netcdf dimension length
character(len=12)      :: dim_name  ! netcdf dimension name
! local grid variables
! kdcorbin, 02/11 - commenting newmap
integer(kind=int_kind), allocatable, dimension(:,:) :: newmap 
         ! map containing subdomain landpoints indexed to nsib vector
integer(kind=int_kind) :: i,j,k                 ! index variables
integer(kind=int_kind) :: ntest1, ntest2, ntest3
integer(kind=int_kind) :: lowerlon, upperlon    ! longitude subdomain limits
integer(kind=int_kind) :: lowerlat, upperlat    ! latitude subdomain limits
integer(kind=int_kind) :: lat_index             ! index value for pbp's
integer(kind=int_kind) :: lon_index             ! index value for pbp's
integer(kind=int_kind), allocatable, dimension(:,:) :: temp_pbp ! temporary array of
                                                            ! pbp coordinates

! subgridding and parallelization variables
real(kind=dbl_kind) ::     dlat   ! latitude gridcell spacing
real(kind=dbl_kind) ::     dlon   ! longitude gridcell spacing
real(kind=dbl_kind) ::     lllat  ! lower left corner of domain
real(kind=dbl_kind) ::     lllon  ! lower left corner of domain
real(kind=real_kind) :: minlon    ! subdomain limits
real(kind=real_kind) :: maxlon
real(kind=real_kind) :: minlat
real(kind=real_kind) :: maxlat
integer(kind=int_kind) :: init_subcount     ! initial subcount before parallel
integer(kind=int_kind), dimension(:), allocatable :: init_subset  ! initial
                                                     ! subset before parallel
integer(kind=int_kind) :: olength
integer(kind=int_kind) :: start_index   ! start index for parallelization
integer(kind=int_kind) :: end_index     ! end index for parallelization
real(kind=dbl_kind), dimension(:,:), allocatable :: lonlatpbp ! temporary array of
                                                    ! pbp coordinates

!     NAMELISTS
!     kdcorbin, 03/11 - removed grid_path
namelist /inlist_sibdrv/ & ! USER DEFINED PARAMETERS
    nsib,ztemp,zwind
namelist /IOLIST_SIBDRV/ & !jk USER SELETED I/O OPTIONS
    param_path, ic_path, dr_format, out_path, qp_path,  &
    pbp_path, co2_path, drvr_type, param_type
namelist /SUBGRID_SIBDRV/ &
    minlon, maxlon, minlat, maxlat
namelist /PBPLIST_SIBDRV/ & ! USER DEFINED PBP DIAGNOSTIC LOCATIONS
    IJTLENsib
namelist /SIBDRV_CONTROL_LIST/ &
    starttime, startyear, endtime, endyear, dtsib, dtsibmetin,  &
    dtsibout, dtsibres, ndtsibpbp, dtsibbcin, roll_respf 

    print *, 'INIT_GRID:'

    !------------------------------------------------------------
    ! read in namel_sibdrv
    !------------------------------------------------------------
    open(unit=2,file='namel_sibdrv',form='formatted')  !jk
    print *,'   reading sib inlist'
    read (2,INLIST_SIBDRV)
    print *,'   reading sib i/olst'      !jk
    read (2,IOLIST_SIBDRV)                         !jk
    print *,'   reading subgrid values'
    read (2,SUBGRID_SIBDRV)
    print *,'   reading sib pbplst'
    read (2,PBPLIST_SIBDRV)
    allocate (lonlatpbp(2,ijtlensib))
    lonlatpbp = 0.0
    read(2,*,err=919)lonlatpbp
    919  continue
    print *,'   reading sib_control_lst'
    read (2,SIBDRV_CONTROL_LIST)
    close(2)
    print *, '   SiB time step (s) = ',dtsib
    if(dtsibout > 0) then
        print *, '   SiB out written (s) = ',dtsibout
    else
        print *, '   SiB out written (months) = ',-dtsibout
    endif
    if(dtsibres > 0) then
        print *, '   SiB restart written (s) = ',dtsibres
    else
        print *, '   SiB restart written (months) = ',-dtsibres
    endif

    histpp = ndtsibpbp /= 0
    !----------------------------------------------------------------
    ! read in grid information from TI file
    !----------------------------------------------------------------
    allocate( latsib(nsib) )
    allocate( lonsib(nsib) )
    allocate( latindex(nsib) )
    allocate( lonindex(nsib) )

    print*,'   nsib= ',nsib
    print*, '   drvr_type= ',drvr_type

    print*,'   reading parameter file: ',trim(param_path)//'_TI.nc'
    status = nf90_open( trim(param_path)//'_TI.nc', nf90_nowrite, ncid )
    if ( status /= nf90_noerr ) call handle_err( status )
    status = nf90_inq_dimid( ncid, 'nsib', dimid )
    if ( status /= nf90_noerr ) call handle_err( status )
    status = nf90_inquire_dimension( ncid, dimid, dim_name,dimlen )
    if ( status /= nf90_noerr ) call handle_err( status )
    if ( dimlen /= nsib ) print *, dimlen, 'and', nsib, "don\'t match"

    ENSURE_VAR( ncid, 'latsib', varid )
    status = nf90_get_var( ncid, varid, latsib )
    ENSURE_VAR( ncid, 'lonsib', varid )
    status = nf90_get_var( ncid, varid, lonsib )

     !kdcorbin, 03/11 - modified reading TI file for single point runs
     if (drvr_type == 'single') then 
         latindex=1
         lonindex=1
     else
         ENSURE_VAR(ncid, 'latindex', varid)
         status = nf90_get_var( ncid, varid, latindex )
         ENSURE_VAR(ncid, 'lonindex', varid)
         status = nf90_get_var( ncid, varid, lonindex )
    endif

    ENSURE_VAR( ncid, 'numlat', varid )
    status = nf90_get_var( ncid, varid, jhr )
    ENSURE_VAR( ncid, 'numlon', varid )
    status = nf90_get_var( ncid, varid, ihr )
    ENSURE_VAR( ncid, 'dlat', varid )
    status = nf90_get_var( ncid, varid, dlat )
    ENSURE_VAR( ncid, 'dlon', varid )
    status = nf90_get_var( ncid, varid, dlon )
    ENSURE_VAR( ncid, 'lllat', varid )
    status = nf90_get_var( ncid, varid, lllat )
    ENSURE_VAR( ncid, 'lllon', varid )
    status = nf90_get_var( ncid, varid, lllon )
    
    status = nf90_close( ncid )

    nhr = ihr * jhr

    !kdcorbin, 02/11 - redefined latitude and longitude arrays
    !allocate( latitude(jhr) )
    !allocate( longitude(ihr) )
    !longitude(1) = lllon
    !do i = 2, ihr
    !    longitude(i) = longitude(i-1) + dlon
    !enddo
    !latitude(1) = lllat
    !do i = 2, jhr
    !    latitude(i) = latitude(i-1) + dlat
    !enddo

    allocate( latitude(nsib) )
    allocate( longitude(nsib) )
    do i=1,nsib
         longitude(i) = lonsib(i)
         latitude(i) = latsib(i)
    enddo

    !---------------------------------------------------------------
    ! calculate subset
    ! kdcorbin, 02/11 - commenting subset for now
    !---------------------------------------------------------------
    !allocate( newmap(ihr,jhr) )
    !newmap(:,:) = 0
    !do i = 1, nsib
    !    newmap( lonindex(i), latindex(i) ) = i
    !enddo
    
    ! convert domain limits to indices
    !lowerlon = int( (minlon-lllon)/dlon + 1 )
    !upperlon = int( (maxlon-lllon)/dlon + 1 )
    !lowerlat = int( (minlat-lllat)/dlat + 1 )
    !upperlat = int( (maxlat-lllat)/dlat + 1 )
    
    ! make sure we stay within the domain
    !if ( lowerlon < 1 ) lowerlon = 1
    !if ( upperlon > ihr ) upperlon = ihr
    !if ( lowerlat < 1 ) lowerlat = 1
    !if ( upperlat > jhr ) upperlat = jhr
    
    ! count number of landpoints in subdomain
    !init_subcount = 0
    !do j = lowerlon, upperlon
    !   do i = lowerlat, upperlat
    !        if ( newmap(j,i) > 0 ) init_subcount = init_subcount + 1
    !    enddo
    !enddo
    

    ! create vector indexing landpoints in subdomain
    !allocate( init_subset(init_subcount) )
    !init_subcount = 0
    !do i = lowerlat, upperlat
    !    do j = lowerlon, upperlon
    !        if ( newmap(j,i) > 0 ) then
    !            init_subcount = init_subcount + 1
    !            init_subset(init_subcount) = newmap(j,i)
    !        endif
    !    enddo
    !enddo
    
    ! calculate subcount for parallelization
    !olength = init_subcount / nchunks
    !subcount = olength
    !if ( nchunks == 1 ) then
    !    subcount = init_subcount
    !elseif ( rank == nchunks ) then
    !    if ( (rank-1)*olength + olength <= init_subcount ) then
    !        subcount = init_subcount - (rank-1)*olength
    !    else
    !        subcount = mod( init_subcount, olength*(nchunks-1) )
    !    endif
    !endif
    !print*, '   nsib=',subcount, 'nsibmax=',nsib
    
    ! calculate starting and ending vertices
    !start_index = (rank-1) * olength + 1
    !end_index = start_index + subcount - 1
   
    ! allocate subset and assign values
    !allocate( subset(subcount) )
    !subset(:) = init_subset( start_index : end_index )
    !deallocate( init_subset )

    ! fill 2d map with new landmask
    !allocate( sublat(subcount) )
    !allocate( sublon(subcount) )
    !newmap(:,:) = 0
    !do i = 1, subcount
    !    newmap(lonindex(subset(i)),latindex(subset(i))) = i
    !    sublat(i) = latindex(subset(i))
    !    sublon(i) = lonindex(subset(i))
    !enddo
    
    !kdcorbin, 02/11 - added following code to maintain original grid
    subcount = nsib
    allocate ( subset(subcount) )
    allocate ( sublat(subcount) )
    allocate ( sublon(subcount) )
    subset = latindex
    !fill 2D map with new landmask
    do i = 1, subcount
         sublat(i) = latindex(subset(i))
         sublon(i) = lonindex(subset(i))
     enddo

    !---------------------------------------------------------------
    ! Find pbp indices and remove duplicates
    !---------------------------------------------------------------
    allocate(temp_pbp(2,ijtlensib))
    temp_pbp(:,:) = 0
    do i = 1, ijtlensib
        ! find latitude index
        lat_index = int( (lonlatpbp(2,i)-lllat)/dlat + 1 )

        ! find longitude index
        lon_index = int( (lonlatpbp(1,i)-lllon)/dlon + 1 )

        if ( lat_index < 1 .or. lat_index > jhr .or.  &
             lon_index < 1 .or. lon_index > ihr ) then
             
            print *, 'Point ', lonlatpbp(1,i), lonlatpbp(2,i),  &
                'is not inside the grid, please fix this.'
            stop
        else  !kdcorbin, 02/11 - removed newmap check
             temp_pbp(1,i) = lon_index
             temp_pbp(2,i) = lat_index
        endif
    enddo

    ! remove duplicates
    do i = 1, ijtlensib-1
        do j = i+1, ijtlensib
            if ( temp_pbp(1,i) == temp_pbp(1,j) .and.  &
                 temp_pbp(2,i) == temp_pbp(2,j) .and.  &
                 temp_pbp(1,i) /= 0 .and. temp_pbp(2,i) /= 0 ) then

                ! duplicate, set second instance to zero
                temp_pbp(1,j) = 0
                temp_pbp(2,j) = 0
                print *, 'point', lonlatpbp(1,j), lonlatpbp(2,j),  &
                    ' is a duplicate of point', lonlatpbp(1,i), lonlatpbp(2,i)
                print *, 'This point has been removed'
            endif
        enddo
    enddo

    ! count remaining number of pbp's, allocate new array and copy subset index
    j = 0
    do i = 1, ijtlensib
        if ( temp_pbp(1,i) /= 0 .and. temp_pbp(2,i) /= 0 )  j = j + 1
    enddo

    if ( j > 0 ) then
        ! there is at least one pbp in the subdomain
        allocate(imultpbpsib(j))
        allocate( latpbp(j) )
        allocate( lonpbp(j) )

        j = 0
        do i = 1, ijtlensib
            if ( temp_pbp(1,i) /= 0 .and. temp_pbp(2,i) /= 0 )  then
                j = j + 1
                !kdcorbin, 02/11 - changed from newmap(temp_pbp(1,i),temp_pbp(2,i))
                imultpbpsib(j) = j
                latpbp(j) = lonlatpbp(2,i)
                lonpbp(j) = lonlatpbp(1,i)
            endif
        enddo
        ijtlensib = j
    else
        histpp = .false.
    endif

    deallocate( temp_pbp )
    deallocate( lonlatpbp )
    !deallocate( newmap )

end subroutine init_grid

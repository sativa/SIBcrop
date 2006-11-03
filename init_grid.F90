subroutine init_grid( rank, nchunks )
!-------------------------------------------------------------
! reads in sibdrv control variables and inputs
! sets up grid
!
! modofications:
#ifdef PGF
use netcdf
use typeSizes
#endif
use kinds
use sib_const_module
use sib_io_module


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
integer(kind=int_kind), allocatable, dimension(:,:) :: newmap ! map containing only
                                                          ! subdomain landpoints
                                                          ! indexed to nsib vector
integer(kind=int_kind) :: i,j,k                 ! index variables
integer(kind=int_kind) :: ntest1, ntest2, ntest3
integer(kind=int_kind) :: lowerlon, upperlon    ! longitude subdomain limits
integer(kind=int_kind) :: lowerlat, upperlat    ! latitude subdomain limits
integer(kind=int_kind) :: lat_index             ! index value for pbp's
integer(kind=int_kind) :: lon_index             ! index value for pbp's
integer(kind=int_kind), allocatable, dimension(:,:) :: temp_pbp ! temporary array of
                                                            ! pbp coordinates

real(kind=real_kind), dimension(:), allocatable ::  &
    areasib,      & !  SiB gridpoint area
    weight_sib,   & !  weight factor
    wgt2_sib,     & !  weight factor
    lat_hr,       & !  used for output
    lon_hr          !  used for output
logical (kind=log_kind), dimension(:), allocatable :: hr_sib        
integer (kind=int_kind), dimension(:), allocatable   ::    &
    link_hr_sib !  link between sib and global points

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
namelist /inlist_sibdrv/ & ! USER DEFINED PARAMETERS
    nsib,ztemp,zwind
namelist /IOLIST_SIBDRV/ & !jk USER SELETED I/O OPTIONS
    param_path, ic_path, dr_format, out_path, qp_path,  &
    pbp_path, co2_path, grid_path, drvr_type
namelist /SUBGRID_SIBDRV/ &
    minlon, maxlon, minlat, maxlat
namelist /PBPLIST_SIBDRV/ & ! USER DEFINED PBP DIAGNOSTIC LOCATIONS
    IJTLENsib
namelist /SIBDRV_CONTROL_LIST/ &
    starttime, startyear, endtime, endyear, dtsib, dtsibmetin,  &
    dtsibout, dtsibres, ndtsibpbp, dtsibbcin, roll_respf 

    print *, 'INIT_GRID:'

    !-----------------------------------------------------------------------
    ! read in namel_sibdrv
    !-----------------------------------------------------------------------
    open(unit=2,file='namel_sibdrv',form='formatted')  !jk
    print *,'\t reading sib inlist'
    read (2,INLIST_SIBDRV)
    print *,'\t reading sib i/olst'      !jk
    read (2,IOLIST_SIBDRV)                         !jk
    print *,'\t reading subgrid values'
    read (2,SUBGRID_SIBDRV)
    print *,'\t reading sib pbplst'
    read (2,PBPLIST_SIBDRV)
    allocate (lonlatpbp(2,ijtlensib))
    lonlatpbp = 0.0
    read(2,*,err=919)lonlatpbp
    919  continue
    print *,'\t reading sib_control_lst'
    read (2,SIBDRV_CONTROL_LIST)
    close(2)
    print *, '\t SiB time step (s) = ',dtsib
    if(dtsibout > 0) then
        print *, '\t SiB out written (s) = ',dtsibout
    else
        print *, '\t SiB out written (months) = ',-dtsibout
    endif
    if(dtsibres > 0) then
        print *, '\t SiB restart written (s) = ',dtsibres
    else
        print *, '\t SiB restart written (months) = ',-dtsibres
    endif

    histpp = ndtsibpbp /= 0
    !-----------------------------------------------------------------------
    ! read in grid information
    !-----------------------------------------------------------------------
    allocate( latsib(nsib) )
    allocate( lonsib(nsib) )
 print*, 'drvr_type=',drvr_type
    if(drvr_type=='single')then
        allocate( areasib(1) )
        allocate(link_hr_sib(1) )
        allocate(weight_sib(1) )
        allocate(hr_sib(1) )
        open(unit=3,file=grid_path,form='formatted') !jk
        read(3,*)ntest1
        if(ntest1 /= 1) stop ' sib_gridmap file no match with model nsib'
        read(3,*)ntest2
        if(ntest2 /= 1) stop ' sib_gridmap file no match with model ihr'
        read(3,*)ntest3
        if(ntest3 /= 1) stop ' sib_gridmap file no match with model jhr'
        read(3,*)lonsib(1)
        read(3,*)latsib(1)
        read(3,*)areasib(1)
        read(3,*)link_hr_sib(1)
        read(3,*)weight_sib(1)         ! jk
        read(3,*)hr_sib(1)        
        ! assign some grid information and exit subroutine
        subcount = 1
        allocate(subset(1))
        allocate(imultpbpsib(1))
        allocate(newmap(1,1))
        allocate(lat_hr(1))
        allocate(lon_hr(1))
        allocate(latpbp(1))
        allocate(lonpbp(1))
        allocate(latindex(1))
        allocate(lonindex(1))
        allocate(sublon(1))
        allocate(sublat(1))
        allocate(latitude(1))
        allocate(longitude(1))
        subset(1) = 1
        imultpbpsib(1) = 1
        ijtlensib = 1
        newmap(1,1) = 1
        lat_hr(1) = lonlatpbp(2,1)
        lon_hr(1) = lonlatpbp(1,1)
        ihr = 1
        jhr = 1
        nhr = 1
        latpbp(1) = lonlatpbp(2,1)
        lonpbp(1) = lonlatpbp(1,1)
        latindex(1) = 1
        lonindex(1) = 1
        sublon(1) = 1
        sublat(1) = 1
        latitude(1) = latsib(1)
        longitude(1) = lonsib(1)
        return
    endif

    allocate( latindex(nsib) )
    allocate( lonindex(nsib) )
    print*, trim(param_path)//'_TI.nc'
    status = nf90_open( trim(param_path)//'_TI.nc', nf90_nowrite, ncid )
    if ( status /= nf90_noerr ) call handle_err( status )
    status = nf90_inq_dimid( ncid, 'nsib', dimid )
    if ( status /= nf90_noerr ) call handle_err( status )
    status = nf90_inquire_dimension( ncid, dimid, dim_name,dimlen )
    if ( status /= nf90_noerr ) call handle_err( status )
    if ( dimlen /= nsib ) print *, dimlen, 'and', nsib, "don\'t match"
    
    status = nf90_inq_varid( ncid, 'latsib', varid )
    status = nf90_get_var( ncid, varid, latsib )
    status = nf90_inq_varid( ncid, 'lonsib', varid )
    status = nf90_get_var( ncid, varid, lonsib )
    status = nf90_inq_varid( ncid, 'latindex', varid )
    status = nf90_get_var( ncid, varid, latindex )
    status = nf90_inq_varid( ncid, 'lonindex', varid )
    status = nf90_get_var( ncid, varid, lonindex )
    status = nf90_inq_varid( ncid, 'numlat', varid )
    status = nf90_get_var( ncid, varid, jhr )
    status = nf90_inq_varid( ncid, 'numlon', varid )
    status = nf90_get_var( ncid, varid, ihr )
    status = nf90_inq_varid( ncid, 'dlat', varid )
    status = nf90_get_var( ncid, varid, dlat )
    status = nf90_inq_varid( ncid, 'dlon', varid )
    status = nf90_get_var( ncid, varid, dlon )
    status = nf90_inq_varid( ncid, 'lllat', varid )
    status = nf90_get_var( ncid, varid, lllat )
    status = nf90_inq_varid( ncid, 'lllon', varid )
    status = nf90_get_var( ncid, varid, lllon )
    
    status = nf90_close( ncid )

    nhr = ihr * jhr

    allocate( latitude(jhr) )
    allocate( longitude(ihr) )
    longitude(1) = lllon
    do i = 2, ihr
        longitude(i) = longitude(i-1) + dlon
    enddo
    latitude(1) = lllat
    do i = 2, jhr
        latitude(i) = latitude(i-1) + dlat
    enddo

    !-----------------------------------------------------------------------
    ! calculate subset
    !-----------------------------------------------------------------------
    allocate( newmap(ihr,jhr) )
    newmap(:,:) = 0
    do i = 1, nsib
        newmap( lonindex(i), latindex(i) ) = i
    enddo
    
    ! convert domain limits to indices
    lowerlon = int( (minlon-lllon)/dlon + 1 )
    upperlon = int( (maxlon-lllon)/dlon + 1 )
    lowerlat = int( (minlat-lllat)/dlat + 1 )
    upperlat = int( (maxlat-lllat)/dlat + 1 )
    
    ! make sure we stay within the domain
    if ( lowerlon < 1 ) lowerlon = 1
    if ( upperlon > ihr ) upperlon = ihr
    if ( lowerlat < 1 ) lowerlat = 1
    if ( upperlat > jhr ) upperlat = jhr
    
    ! count number of landpoints in subdomain
    init_subcount = 0
    do j = lowerlon, upperlon
        do i = lowerlat, upperlat
            if ( newmap(j,i) > 0 ) init_subcount = init_subcount + 1
        enddo
    enddo
    

    ! create vector indexing landpoints in subdomain
    allocate( init_subset(init_subcount) )
    init_subcount = 0
    do i = lowerlat, upperlat
        do j = lowerlon, upperlon
            if ( newmap(j,i) > 0 ) then
                init_subcount = init_subcount + 1
                init_subset(init_subcount) = newmap(j,i)
            endif
        enddo
    enddo
    
    ! calculate subcount for parallelization
    olength = init_subcount / nchunks
    subcount = olength
    if ( nchunks == 1 ) then
        subcount = init_subcount
    elseif ( rank == nchunks ) then
        if ( (rank-1)*olength + olength <= init_subcount ) then
            subcount = init_subcount - (rank-1)*olength
        else
            subcount = mod( init_subcount, olength*(nchunks-1) )
        endif
    endif
    print*, '\t nsib=',subcount, 'nsibmax=',nsib
    
    ! calculate starting and ending vertices
    start_index = (rank-1) * olength + 1
    end_index = start_index + subcount - 1
   
    ! allocate subset and assign values
    allocate( subset(subcount) )
    subset(:) = init_subset( start_index : end_index )
    deallocate( init_subset )

    ! fill 2d map with new landmask
    allocate( sublat(subcount) )
    allocate( sublon(subcount) )
    newmap(:,:) = 0
    do i = 1, subcount
        newmap(lonindex(subset(i)),latindex(subset(i))) = i
        sublat(i) = latindex(subset(i))
        sublon(i) = lonindex(subset(i))
    enddo
    
   
    !-----------------------------------------------------------------------
    ! Find pbp indices and remove duplicates
    !-----------------------------------------------------------------------
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
        endif
        if ( newmap(lon_index,lat_index) /= 0 ) then
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
                imultpbpsib(j) = newmap(temp_pbp(1,i),temp_pbp(2,i))
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
    deallocate( newmap )

end subroutine init_grid

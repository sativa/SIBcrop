program sibmerge
!-------------------------------------------------------------------------------
! Author  :  Owen Leonard
! Date    :  May 19, 2004
! Purpose :
!   This program merges the output files from SiB3 after a simulation has run 
! using multiple processes.  There are eight command line arguments that must 
! be given in the following order.
!  1. # of processes  (second command line argument to SiBD3)
!  2. first year of simulation (YYYY - four digits)
!  3. first month of simulation (MM - two digits)
!  4. last year of simulation (YYYY - four digits)
!  5. last month of simulation (MM - two digits)
!  6. path of input files ending in '/'
!  7. path of output files ending in '/' (./ works for current directory)
!  8. merger restart files only?  value must be 'restart', otherwise it merges
!     all files
! Notes   :
!   This program requires the netcdf library.
!   
!   I put in the restart_only flag because restart files are written at the
!  beginning of the month, while qp and pbp files are finished at the end of
!  the month.  If I stop in October, I will use the October restart files.
!  However, the qp and pbp files don't have any values saved to them, so I 
!  can't merge them.  If I tried to, I would get a segmentation fault.  So,
!  in this case, I would only merge the restart files.
!-------------------------------------------------------------------------------

implicit none

! input variables
integer :: nprocs               ! # of processes
integer :: start_year           ! beginning year of simulation
integer :: start_month          ! beginning month of simulation
integer :: end_year             ! ending year of simulation
integer :: end_month            ! ending month of simulation
character(len=256) :: in_path   ! directory where input files reside
character(len=256) :: out_path  ! directory where output files will be written
character(len=7)   :: restart   ! merge restart files only? must equal 'restart'

! file variables
character(len=256) :: outfile   ! name of output file
character(len=256) :: var       ! temp variable holding value of command line args
character(len=256), dimension(:,:), allocatable :: filenames
integer :: p, m, y          ! index variables
logical :: restart_only     ! flag determining if restart files are the only
                            ! files to be merged or not
integer*1 :: iargc          ! intrinsic function that counts # command line args

    ! check that all of the command line arguments exist
    if ( IArgC() < 7 ) then
        print *, 'There must be at least 7 command line arguments in the following order'
        print *, "  1. # of processes  (second command line argument to SiBD3)"
        print *, "  2. first year of simulation (YYYY - four digits)"
        print *, "  3. first month of simulation (MM - two digits)"
        print *, "  4. last year of simulation (YYYY - four digits)"
        print *, "  5. last month of simulation (MM - two digits)"
        print *, "  6. path of input files ending in '/'"
        print *, "  7. path of output files ending in '/'"
        print *, "\t(./ works for current directory)"
        print *, "  8. (optional) merge restart files only?  value must be "
        print *, "\t 'restart', otherwise it merges qp and pbp files as well"
        stop
    endif

    ! retrieve required information from command line
    call getarg( 1, var )
    read( var, * ) nprocs
    call getarg( 2, var )
    read( var, * ) start_year
    call getarg( 3, var )
    read( var, * ) start_month
    call getarg( 4, var )
    read( var, * ) end_year
    call getarg( 5, var )
    read( var, * ) end_month
    call getarg( 6, in_path )
    call getarg( 7, out_path )
    if ( IArgC() >= 8 ) then
        call getarg( 8, restart )
    else
        restart = 'no'
    endif
    
    ! check input 
    if ( nprocs < 1 ) stop 'nprocs must be >= 1'
    if ( start_year > end_year ) stop "can\'t run backward in time"
    if ( start_year == end_year .and. start_month > end_month )  &
        stop "can\'t run backward in time"
    if ( len(trim(in_path)) == 0 )  &
        print *, 'in_path was empty, assuming that files are located'//  &
                 ' in current directory'
    if ( len(trim(out_path)) == 0 )  &
        print *, 'out_path was empty, resulting files will be '//  &
                 'written to current directory'
    if ( restart == 'restart' ) then
        restart_only = .true.
    else
        restart_only = .false.
    endif

    allocate( filenames(nprocs,2) )

    do y = start_year, end_year
        do m = 1, 12

            ! don't start run until start_month
            if ( y == start_year .and. m < start_month ) cycle
            
            ! end simulation after end_month
            if ( y == end_year .and. m > end_month ) exit
            
            print *, 'Processing ', m, y
            
            ! create restart filenames and combine them
            print *, '  combining restart files'
            do p = 1, nprocs
                write( filenames(p,1), '(a,i4.4,i2.2,a,i3.3,a)' )  &
                    trim(in_path)//'sib_r', y, m, 'p', p, '.nc'
            enddo
            write( outfile, '(a,i4.4,i2.2,a)' ) trim(out_path)//'sib_r',  &
                y, m, '.nc'
            call restartmerge( nprocs, filenames(:,1), outfile )

            if ( .not. restart_only ) then
                ! create qp filenames
                do p = 1, nprocs
                    write( filenames(p,1), '(a,i4.4,i2.2,a,i3.3,a)' )  &
                        trim(in_path)//'hsib_', y, m, 'p', p, '.qp2.nc'
                    write( filenames(p,2), '(a,i4.4,i2.2,a,i3.3,a)' )  &
                        trim(in_path)//'hsib_', y, m, 'p', p, '.qp3.nc'
                enddo
                write( outfile, '(a,i4.4,i2.2)' ) trim(out_path)//'hsib_', y, m

                ! combine qp2 files
                print *, '  combining qp2 files'
                call qp2merge( nprocs, filenames(:,1), outfile )

                ! combine qp3 files
                print *, '  combining qp3 files'
                call qp3merge( nprocs, filenames(:,2), outfile )

                ! create pbp filenames and combine them
                do p = 1, nprocs
                    write( filenames(p,1), '(a,i4.4,i2.2,a,i3.3,a)' )  &
                        trim(in_path)//'psib_', y, m, 'p', p, '.pbp2.nc'
                    write( filenames(p,2), '(a,i4.4,i2.2,a,i3.3,a)' )  &
                        trim(in_path)//'psib_', y, m, 'p', p, '.pbp3.nc'
                enddo

                ! combine pbp2 files
                print *, '  combining pbp2 files'
                write( outfile, '(a,i4.4,i2.2,a)' ) trim(out_path)//'psib_',  &
                    y, m, '.pbp2.nc'
                call pbpmerge( nprocs, filenames(:,1), outfile )

                ! combine pbp3 files
                print *, '  combining pbp3 files'
                write( outfile, '(a,i4.4,i2.2,a)' ) trim(out_path)//'psib_',  &
                    y, m, '.pbp3.nc'
                call pbpmerge( nprocs, filenames(:,2), outfile )
                
                if ( m == 12 ) then
                    ! create respfactor filenames and combine them
                    do p = 1, nprocs
                        write( filenames(p,1), '(a,i4.4,a,i3.3,a)' )  &
                            trim(in_path)//'CO2_respf_', y, 'p', p, '_new'
                    enddo

                    ! combine respfactor files
                    print *, '  combining respfactor files'
                    write( outfile, '(a,i4.4,a)') trim(out_path)//'CO2_respf_',  &
                        y, '_new'
                    call respfmerge( nprocs, filenames(:,1), outfile )
                endif
            endif
            
        enddo
    enddo

    deallocate( filenames )


end program sibmerge

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine qp2merge( num, filenames, out_path )

use netcdf
use typeSizes

#include "nc_util.h"

implicit none

! parameters
integer, intent(in) :: num                                   ! # processes
character(len=256), dimension(num), intent(in) :: filenames  ! input file names
character(len=256), intent(in) :: out_path                   ! output file name

! netcdf variables
integer :: status                              ! return value of netcdf function
integer, dimension(num) :: ncid                ! input file id#s
integer :: outid                               ! output file id#
integer :: unlimitid                           ! unlimited dimension id#
integer :: numvars                             ! # variables in input files
integer :: numatts                             ! # attributes in input files
character(len=20) :: name                      ! variable name
integer, dimension(:), allocatable :: dimid    ! dimension id#s
integer, dimension(:), allocatable :: vardims  ! variable dimensions
integer :: dimlen                              ! dimension length
integer :: timelen                             ! dimension length
integer :: landlen                             ! dimension length
integer :: latlen                              ! dimension length
integer :: lonlen                              ! dimension length
integer :: xtype                               ! variable data type
integer :: numdims                             ! # dimensions in input files
integer :: varid                               ! variable id#

! local variables
character(len=256) :: outfile                  ! output file name
character(len=256) :: command                  ! output file name
double precision, dimension(:), allocatable :: dvalues         ! input data values
character(len=10), dimension(:), allocatable :: cvalues        ! output data values
real, dimension(:,:), allocatable :: values   ! input data values
real, dimension(:,:), allocatable :: result   ! output data values
integer, dimension(:), allocatable :: lonindex
integer, dimension(:), allocatable :: latindex
integer :: f, n, a, d                          ! index variables
integer :: landpoints

    ! create output file name and create file
    write( outfile, '(a,i1.1,a)' ) trim(out_path)//'.qp', 2, '.nc'

    ! copy file if only one process, then return out of subroutine
    if ( num == 1 ) then
        command = "cp "//trim(filenames(1))//" "//trim(outfile)
        call system( trim(command) )
        return
    endif

    ! open all input files
    do f = 1, num
        CHECK( nf90_open( trim(filenames(f)), nf90_nowrite, ncid(f) ) )
    enddo
    status = nf90_inquire( ncid(1), nDimensions=numdims, nVariables=numvars,  &
        nAttributes=numatts, unlimitedDimId=unlimitid )
    CHECK( status )
    
    CHECK( nf90_create( trim(outfile), nf90_clobber, outid ) )
    ! copy global attributes and dimensions
    !   over to output file
    do a = 1, numatts
        CHECK( nf90_inq_attname( ncid(1), nf90_global, a, name ) )
        status = nf90_copy_att( ncid(1), nf90_global, trim(name),  &
            outid, nf90_global )
        CHECK( status )
    enddo
    allocate( dimid(numdims) )
    do d = 1, numdims - 1
        if ( d /= unlimitid ) then
            CHECK( nf90_inquire_dimension( ncid(1), d, name=name, len=dimlen ) )
            CHECK( nf90_def_dim( outid, trim(name), dimlen, dimid(d) ) )
        else
            CHECK( nf90_def_dim( outid, 'time', nf90_unlimited, dimid(d) ) )
        endif
    enddo
    landpoints = 0
    do f = 1, num
        CHECK( nf90_inquire_dimension( ncid(f), numdims, name=name, len=dimlen ) )
        landpoints = landpoints + dimlen
    enddo
    CHECK( nf90_def_dim( outid, trim(name), landpoints, dimid(numdims) ) )
    
    ! define variables in the files (excluding lonindex and latindex)
    allocate( vardims(numdims) )
    do n = 1, 7
        status = nf90_inquire_variable( ncid(1), n, name=name, xtype=xtype,  &
            ndims=numdims, dimids=vardims, natts=numatts )
        CHECK( status )
        status = nf90_def_var( outid, trim(name), xtype,  &
            vardims(1:numdims), varid )
        CHECK( status )

        ! copy attributes over
        do a = 1, numatts
            CHECK( nf90_inq_attname( ncid(1), n, a, name=name ) )
            CHECK( nf90_copy_att( ncid(1), n, trim(name), outid, n ) )
        enddo
    enddo
    status = nf90_inquire_variable( ncid(1), 8, name=name, xtype=xtype,  &
        natts=numatts )
    CHECK( status )
    CHECK( nf90_def_var( outid, trim(name), xtype, (/5/), varid ) )
    do a = 1, numatts
        CHECK( nf90_inq_attname( ncid(1), 8, a, name=name ) )
        CHECK( nf90_copy_att( ncid(1), 8, trim(name), outid, varid ) )
    enddo
    status = nf90_inquire_variable( ncid(1), 9, name=name, xtype=xtype,  &
        natts=numatts )
    CHECK( status )
    CHECK( nf90_def_var( outid, trim(name), xtype, (/5/), varid ) )
    do a = 1, numatts
        CHECK( nf90_inq_attname( ncid(1), 9, a, name=name ) )
        CHECK( nf90_copy_att( ncid(1), 9, trim(name), outid, varid ) )
    enddo
    do n = 10, numvars
        status = nf90_inquire_variable( ncid(1), n, name=name, xtype=xtype,  &
            natts=numatts )
        CHECK( status )
        CHECK( nf90_def_var( outid, trim(name), xtype, (/5,1/), varid ) )
        ! copy attributes over
        do a = 1, numatts
            CHECK( nf90_inq_attname( ncid(1), n, a, name=name ) )
            CHECK( nf90_copy_att( ncid(1), n, trim(name), outid, varid ) )
        enddo
    enddo
    
    CHECK( nf90_enddef( outid ) )
    
    ! copy 1-D and 2-D variables over
    do n = 1, 7
        status = nf90_inquire_variable( ncid(1), n, xtype=xtype,  &
            ndims=numdims, dimids=vardims )
        CHECK( status )
        if ( xtype == nf90_double ) then
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=dimlen ) )
            allocate( dvalues(dimlen) )
            CHECK( nf90_get_var( ncid(1), n, dvalues ) )
            CHECK( nf90_put_var( outid, n, dvalues ) )
            deallocate( dvalues )
        else if ( xtype == nf90_char ) then
            CHECK( nf90_inquire_dimension( ncid(1), vardims(2), len=dimlen ) )
            allocate( cvalues(dimlen) )
            CHECK( nf90_get_var( ncid(1), n, cvalues ) )
            CHECK( nf90_put_var( outid, n, cvalues ) )
            deallocate( cvalues )
        else 
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=dimlen ) )
            allocate( values(dimlen,1) )
            CHECK( nf90_get_var( ncid(1), n, values ) )
            CHECK( nf90_put_var( outid, n, values ) )
            deallocate( values )
        endif
    enddo
    
    ! copy lonindex and latindex values over
    allocate( lonindex(landpoints) )
    allocate( latindex(landpoints) )
    a = 1
    do f = 1, num
        CHECK( nf90_inquire_dimension( ncid(f), 5, len=landlen ) )
        CHECK( nf90_get_var( ncid(f), 8, lonindex(a:a+landlen-1) ) )
        CHECK( nf90_get_var( ncid(f), 9, latindex(a:a+landlen-1) ) )
        a = a + landlen
    enddo
    
    CHECK( nf90_put_var( outid, 8, lonindex ) )
    CHECK( nf90_put_var( outid, 9, latindex ) )
    deallocate( lonindex )
    deallocate( latindex )
    
    ! copy variables over to output file
    CHECK( nf90_inquire_dimension( ncid(1), 1, len=timelen ) )
    allocate( result(landpoints,timelen) )
    do n = 10, numvars
        a = 1
        do f = 1, num
            CHECK( nf90_inquire_dimension( ncid(f), 5, len=landlen ) )
            CHECK( nf90_get_var( ncid(f), n, result(a:a+landlen-1,:) ) )
            a = a + landlen
        enddo
        CHECK( nf90_put_var( outid, n, result ) )
    enddo
    deallocate( result )
    deallocate( dimid )
    deallocate( vardims )

    CHECK( nf90_close( outid ) )

end subroutine qp2merge

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine qp3merge( num, filenames, out_path )

use netcdf
use typeSizes
implicit none

! parameters
integer, intent(in) :: num                                   ! # processes
character(len=256), dimension(num), intent(in) :: filenames  ! input file names
character(len=256), intent(in) :: out_path                   ! output file name

! netcdf variables
integer :: status                              ! return value of netcdf function
integer, dimension(num) :: ncid                ! input file id#s
integer :: outid                               ! output file id#
integer :: unlimitid                           ! unlimited dimension id#
integer :: numvars                             ! # variables in input files
integer :: numatts                             ! # attributes in input files
character(len=20) :: name                      ! variable name
integer, dimension(:), allocatable :: dimid    ! dimension id#s
integer, dimension(:), allocatable :: vardims  ! variable dimensions
integer :: dimlen                              ! dimension length
integer :: timelen                             ! dimension length
integer :: landlen                             ! dimension length
integer :: latlen                              ! dimension length
integer :: lonlen                              ! dimension length
integer :: levlen                              ! dimension length
integer :: xtype                               ! variable data type
integer :: numdims                             ! # dimensions in input files
integer :: varid                               ! variable id#

! local variables
character(len=256) :: outfile                  ! output file name
double precision, dimension(:), allocatable :: dvalues         ! input data values
character(len=10), dimension(:), allocatable :: cvalues        ! output data values
real, dimension(:,:,:), allocatable :: values   ! input data values
real, dimension(:,:,:), allocatable :: result   ! output data values
integer, dimension(:), allocatable :: lonindex
integer, dimension(:), allocatable :: latindex
integer :: f, n, a, d                          ! index variables
integer :: landpoints
character(len=256) :: command

    ! create output file name and create file
    write( outfile, '(a,i1.1,a)' ) trim(out_path)//'.qp', 3, '.nc'

    ! copy file if only one process, then return out of subroutine
    if ( num == 1 ) then
        command = "cp "//trim(filenames(1))//" "//trim(outfile)
        call system( trim(command) )
        return
    endif

    ! open all input files
    do f = 1, num
        CHECK( nf90_open( trim(filenames(f)), nf90_nowrite, ncid(f) ) )
    enddo
    status = nf90_inquire( ncid(1), nDimensions=numdims, nVariables=numvars,  &
        nAttributes=numatts, unlimitedDimId=unlimitid )
    CHECK( status )

    
    CHECK( nf90_create( trim(outfile), nf90_clobber, outid ) )
    ! copy global attributes and dimensions 
    !   over to output file
    do a = 1, numatts
        CHECK( nf90_inq_attname( ncid(1), nf90_global, a, name ) )
        status = nf90_copy_att( ncid(1), nf90_global, trim(name),  &
            outid, nf90_global )
        CHECK( status )
    enddo
    allocate( dimid(numdims) )
    do d = 1, numdims - 1
        if ( d /= unlimitid ) then
            status = nf90_inquire_dimension( ncid(1), d, name=name, len=dimlen )
            CHECK( status )
            CHECK( nf90_def_dim( outid, trim(name), dimlen, dimid(d) ) )
        else
            CHECK( nf90_def_dim( outid, 'time', nf90_unlimited, dimid(d) ) )
        endif
    enddo
    landpoints = 0
    do f = 1, num
        status = nf90_inquire_dimension( ncid(f), numdims, name=name, len=dimlen )
        CHECK( status )
        landpoints = landpoints + dimlen
    enddo
    CHECK( nf90_def_dim( outid, trim(name), landpoints, dimid(numdims) ) )
    
    ! define variables in the files 
    allocate( vardims(numdims) )
    do n = 1, 8
        status = nf90_inquire_variable( ncid(1), n, name=name, xtype=xtype,  &
            ndims=numdims, dimids=vardims, natts=numatts )
        CHECK( status )
        status = nf90_def_var( outid, trim(name), xtype,  &
            vardims(1:numdims), varid )
        CHECK( status )

        ! copy attributes over
        do a = 1, numatts
            CHECK( nf90_inq_attname( ncid(1), n, a, name=name ) )
            CHECK( nf90_copy_att( ncid(1), n, trim(name), outid, n ) )
        enddo
    enddo
    status = nf90_inquire_variable( ncid(1), 9, name=name, xtype=xtype,  &
        natts=numatts )
    CHECK( status )
    CHECK( nf90_def_var( outid, trim(name), xtype, (/6/), varid ) )
    do a = 1, numatts
        CHECK( nf90_inq_attname( ncid(1), 9, a, name=name ) )
        CHECK( nf90_copy_att( ncid(1), 9, trim(name), outid, varid ) )
    enddo
    status = nf90_inquire_variable( ncid(1), 10, name=name, xtype=xtype,  &
        natts=numatts )
    CHECK( status )
    CHECK( nf90_def_var( outid, trim(name), xtype, (/6/), varid ) )
    do a = 1, numatts
        CHECK( nf90_inq_attname( ncid(1), 10, a, name=name ) )
        CHECK( nf90_copy_att( ncid(1), 10, trim(name), outid, varid ) )
    enddo
    do n = 11, numvars
        status = nf90_inquire_variable( ncid(1), n, name=name, xtype=xtype,  &
            natts=numatts )
        CHECK( status )
        CHECK( nf90_def_var( outid, trim(name), xtype, (/6,5,1/), varid ) )
        ! copy attributes over
        do a = 1, numatts
            CHECK( nf90_inq_attname( ncid(1), n, a, name=name ) )
            CHECK( nf90_copy_att( ncid(1), n, trim(name), outid, varid ) )
        enddo
    enddo
    
    CHECK( nf90_enddef( outid ) )
    
    ! copy 1-D and 2-D variables over
    do n = 1, 8
        status = nf90_inquire_variable( ncid(1), n, xtype=xtype,  &
            ndims=numdims, dimids=vardims )
        CHECK( status )
        if ( xtype == nf90_double ) then
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=dimlen ) )
            allocate( dvalues(dimlen) )
            CHECK( nf90_get_var( ncid(1), n, dvalues ) )
            CHECK( nf90_put_var( outid, n, dvalues ) )
            deallocate( dvalues )
        else if ( xtype == nf90_char ) then
            CHECK( nf90_inquire_dimension( ncid(1), vardims(2), len=dimlen ) )
            allocate( cvalues(dimlen) )
            CHECK( nf90_get_var( ncid(1), n, cvalues ) )
            CHECK( nf90_put_var( outid, n, cvalues ) )
            deallocate( cvalues )
        else 
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=dimlen ) )
            allocate( values(dimlen,1,1) )
            CHECK( nf90_get_var( ncid(1), n, values ) )
            CHECK( nf90_put_var( outid, n, values ) )
            deallocate( values )
        endif
    enddo
    
    ! copy lonindex and latindex values over
    allocate( lonindex(landpoints) )
    allocate( latindex(landpoints) )
    a = 1
    do f = 1, num
        CHECK( nf90_inquire_dimension( ncid(f), 6, len=landlen ) )
        CHECK( nf90_get_var( ncid(f), 9, lonindex(a:a+landlen-1) ) )
        CHECK( nf90_get_var( ncid(f), 10, latindex(a:a+landlen-1) ) )
        a = a + landlen
    enddo
    
    CHECK( nf90_put_var( outid, 9, lonindex ) )
    CHECK( nf90_put_var( outid, 10, latindex ) )
    deallocate( lonindex )
    deallocate( latindex )
    
    ! copy variables over to output file
    CHECK( nf90_inquire_dimension( ncid(1), 1, len=timelen ) )
    CHECK( nf90_inquire_dimension( ncid(1), 5, len=levlen ) )
    allocate( result(landpoints,levlen,timelen) )
    do n = 11, numvars
        a = 1
        do f = 1, num
            CHECK( nf90_inquire_dimension( ncid(f), 6, len=landlen ) )
            CHECK( nf90_get_var( ncid(f), n, result(a:a+landlen-1,:,:) ) )
            a = a + landlen
        enddo
        CHECK( nf90_put_var( outid, n, result ) )
    enddo
    deallocate( result )
    deallocate( dimid )
    deallocate( vardims )

    CHECK( nf90_close( outid ) )

end subroutine qp3merge

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine restartmerge( num, filenames, out_path )

use netcdf
use typeSizes
implicit none

! parameters
integer, intent(in) :: num                                   ! # processes
character(len=256), dimension(num), intent(in) :: filenames  ! input file names
character(len=256), intent(in) :: out_path                   ! output directory

! netcdf variables
integer :: status                              ! return value of netcdf function
integer, dimension(num) :: ncid                ! input file id#s
integer :: outid                               ! output file id#
integer :: varid                               ! variable id#
integer, dimension(:), allocatable :: dimid    ! dimension id#s
integer :: dimlen                              ! dimension length
integer :: numdims                             ! # dimensions in input files
integer :: numvars                             ! # variables in input files
integer, dimension(:), allocatable :: vardims  ! variable dimensions
character(len=20) :: name                      ! variable name
integer :: xtype                               ! variable data type

! local variables
integer :: d,v,f,n                                  ! index variables
character(len=256) :: outfile                       ! output file name
integer, dimension(:), allocatable :: lengths       ! dimension lengths
integer, dimension(:), allocatable :: intvalues     ! input data values
integer, dimension(:), allocatable :: intresults    ! output data values
real, dimension(:,:,:), allocatable :: realvalues     ! input data values
real, dimension(:,:,:), allocatable :: realresults    ! output data values
character(len=256) :: command


    ! copy file if only one process, then return out of subroutine
    if ( num == 1 ) then
        command = "cp "//trim(filenames(1))//" "//trim(out_path)
        call system( trim(command) )
        return
    endif

    ! open first file and find out how many dimensions it has
    status = nf90_open( trim(filenames(1)), nf90_nowrite, ncid(1) )
    if ( status == 2 ) then
        print *, "restart files don\'t exist for specified year and month"
        return
    endif
    CHECK( nf90_inquire( ncid(1), nDimensions=numdims, nVariables=numvars ) )

    ! create output file
    CHECK( nf90_create( trim(out_path), nf90_clobber, outid ) )

    ! copy dimensions over to output file
    allocate( dimid(numdims) )
    do d = 1, numdims
        CHECK( nf90_inquire_dimension( ncid(1), d, name=name, len=dimlen ) )
        CHECK( nf90_def_dim( outid, trim(name), dimlen, dimid(d) ) )
    enddo

    ! define variables in the files 
    allocate( vardims(numdims) )
    allocate( lengths(numdims) )
    do n = 1, numvars
        status = nf90_inquire_variable( ncid(1), n, name=name, xtype=xtype,  &
            ndims=numdims, dimids=vardims )
        CHECK( status )
        status = nf90_def_var( outid, trim(name), xtype,  &
            vardims(1:numdims), varid )
        CHECK( status )
    enddo
    
    CHECK( nf90_enddef( outid ) )

    ! copy non-variant values to output file
    allocate( intvalues(1) )
    allocate( realvalues(1,1,1) )
    CHECK( nf90_get_var( ncid(1), 1, intvalues ) )
    CHECK( nf90_put_var( outid, 1, intvalues ) )
    CHECK( nf90_get_var( ncid(1), 2, intvalues ) )
    CHECK( nf90_put_var( outid, 2, intvalues ) )
    CHECK( nf90_get_var( ncid(1), 3, intvalues ) )
    CHECK( nf90_put_var( outid, 3, intvalues ) )
    CHECK( nf90_get_var( ncid(1), 4, realvalues ) )
    CHECK( nf90_put_var( outid, 4, realvalues ) )
    CHECK( nf90_get_var( ncid(1), 5, intvalues ) )
    CHECK( nf90_put_var( outid, 5, intvalues ) )
    
    
    ! open remaining restart files to merge
    do f = 2, num
        CHECK( nf90_open( trim(filenames(f)), nf90_nowrite, ncid(f) ) )
    enddo
    
    ! add up subcount and write to file
    allocate( intresults(1) )
    intresults = 0
    do f = 1, num
        CHECK( nf90_get_var( ncid(f), 6, intvalues ) )
        intresults = intresults + intvalues
    enddo
    CHECK( nf90_put_var( outid, 6, intresults ) )
  
    ! add up arrays and write out to file
    deallocate( intvalues )
    deallocate( realvalues )
    deallocate( intresults )
    do n = 7, numvars
        status = nf90_inquire_variable( ncid(1), n, ndims=numdims,  &
             dimids=vardims, xtype=xtype, name=name )
        CHECK( status )
        if ( xtype == nf90_int ) then
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=lengths(1) ) )
            allocate( intvalues(lengths(1)) )
            allocate( intresults(lengths(1)) )
            intresults(:) = 0
            do f = 1, num
                CHECK( nf90_get_var( ncid(f), n, intvalues(:) ) )
                intresults = intresults + intvalues
            enddo
            CHECK( nf90_put_var( outid, n, intresults(:) ) )
            deallocate( intvalues )
            deallocate( intresults )
        else if ( trim(name) == 'tot_an' ) then
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=lengths(1) ) )
            CHECK( nf90_inquire_dimension( ncid(1), vardims(2), len=lengths(2) ) )
            allocate( realvalues(lengths(1), lengths(2), 1) )
            allocate( realresults(lengths(1), lengths(2), 1) )
            realresults(:,:,:) = 1.e36
            do f = 1, num
                CHECK( nf90_get_var( ncid(f), n, realvalues ) )
                do v = 1, lengths(2)
                    if ( realvalues(1,v,1) /= 1.e36 )  &
                        realresults(:,v,1) = realvalues(:,v,1)
                enddo
            enddo
            CHECK( nf90_put_var( outid, n, realresults ) )
            deallocate( realvalues )
            deallocate( realresults )
        else if ( trim(name) == 'tot_ss' ) then
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=lengths(1) ) )
            CHECK( nf90_inquire_dimension( ncid(1), vardims(2), len=lengths(2) ) )
            CHECK( nf90_inquire_dimension( ncid(1), vardims(3), len=lengths(3) ) )
            allocate( realvalues(lengths(1), lengths(2), lengths(3)) )
            allocate( realresults(lengths(1), lengths(2), lengths(3)) )
            realresults(:,:,:) = 1.e36
            do f = 1, num
                CHECK( nf90_get_var( ncid(f), n, realvalues ) )
                do v = 1, lengths(2)
                    if ( realvalues(1,v,1) /= 1.e36 )  &
                        realresults(:,v,:) = realvalues(:,v,:)
                enddo
            enddo
            CHECK( nf90_put_var( outid, n, realresults ) )
            deallocate( realvalues )
            deallocate( realresults )
        else
            CHECK( nf90_inquire_dimension( ncid(1), vardims(1), len=lengths(1) ) )
            if ( numdims == 1 ) then
                lengths(2) = 1
            else
               status = nf90_inquire_dimension( ncid(1), vardims(2),  &
                    len=lengths(2) )
               CHECK( status )
            endif
            allocate( realvalues(lengths(1),lengths(2),1) )
            allocate( realresults(lengths(1),lengths(2),1) )
            realresults(:,:,:) = 1.e36
            do f = 1, num
                CHECK( nf90_get_var( ncid(f), n, realvalues ) )
                do v = 1, lengths(1)
                    if ( realvalues(v,1,1) /= 1.e36 )  &
                        realresults(v,:,1) = realvalues(v,:,1)
                enddo
            enddo
            CHECK( nf90_put_var( outid, n, realresults ) )
            deallocate( realvalues )
            deallocate( realresults )
        endif
    enddo
    
    ! close all files
    do f = 1, num
        status = nf90_close( ncid(f) )
    enddo
    CHECK( nf90_close( outid ) )

    deallocate( dimid )
    deallocate( vardims )
    deallocate( lengths )

end subroutine restartmerge


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine pbpmerge( num, filenames, outfile )

use netcdf
use typeSizes
implicit none

! parameters
integer, intent(in) :: num                                   ! # processes
character(len=256), dimension(num), intent(in) :: filenames  ! input file names
character(len=256), intent(in) :: outfile                    ! output directory

! netcdf variables
integer :: status
integer, dimension(num) :: ncid
integer :: outid
integer :: varid
integer :: numatts
integer :: numdims
integer :: vdims
integer :: npoints
integer :: numvars
integer :: unlimitid
integer :: dimlen
character(len=20) :: name
integer, dimension(:), allocatable :: dimid
integer, dimension(:), allocatable :: vardims
integer :: dimcount
integer :: xtype

! local variables
integer :: a, d, f, e, n, v
logical, dimension(num) :: exists
integer, dimension(:), allocatable :: ivalues
character(len=10), dimension(:), allocatable :: cvalues        ! output data values
double precision, dimension(:), allocatable :: dvalues
real, dimension(:,:,:), allocatable :: fvalues
real, dimension(:,:,:), allocatable :: fresults
character(len=256) :: command

    ! copy file if only one process, then return out of subroutine
    if ( num == 1 ) then
        command = "cp "//trim(filenames(1))//" "//trim(outfile)
        call system( trim(command) )
        return
    endif

    ! open existing pbp files
    e = 0
    do f = 1, num
        status = nf90_open( trim(filenames(f)), nf90_nowrite, ncid(f) )
        if ( status == 2 ) then
            exists(f) = .false.
        else
            exists(f) = .true.
            e = f
        endif
    enddo

    ! make sure there is at least one existing pbp file
    if ( e == 0 ) then
        print *, 'pbp files do not exists for this month'
        return
    endif
    
    ! find out information from input files
    status = nf90_inquire( ncid(e), nDimensions=numdims, nVariables=numvars,  &
        nAttributes=numatts, unlimitedDimId=unlimitid )
    CHECK( status )

    ! create output file
    CHECK( nf90_create( trim(outfile), nf90_clobber, outid ) )

    ! copy global attributes over to output file
    do a = 1, numatts
        CHECK( nf90_inq_attname( ncid(e), nf90_global, a, name ) )
        status = nf90_copy_att( ncid(e), nf90_global, trim(name),  &
            outid, nf90_global )
        CHECK( status )
    enddo
    
    ! copy dimensions over to output file
    allocate( dimid(numdims) )
    CHECK( nf90_def_dim( outid, 'time', nf90_unlimited, dimid(1) ) )
    CHECK( nf90_def_dim( outid, 'char_len', 10, dimid(2) ) )
    dimcount = 0
    do f = 1, num
        if ( exists(f) ) then
            status = nf90_inquire_dimension( ncid(f), 3, name=name,  &
                len=dimlen )
            CHECK( status )
            dimcount = dimcount + dimlen
        endif
    enddo
    CHECK( nf90_def_dim( outid, trim(name), dimcount, dimid(3) ) )
    if ( numdims == 4 )  &
        CHECK( nf90_def_dim( outid, 'level', 10, dimid(4) ) )
    
    ! define variables in the files 
    allocate( vardims(numdims) )
    do n = 1, numvars

        status = nf90_inquire_variable( ncid(e), n, name=name, xtype=xtype,  &
            ndims=vdims, dimids=vardims, natts=numatts )
        CHECK( status )
        status = nf90_def_var( outid, trim(name), xtype,  &
            vardims(1:vdims), varid )
        CHECK( status )

        ! copy attributes over
        do a = 1, numatts
            CHECK( nf90_inq_attname( ncid(e), n, a, name=name ) )
            CHECK( nf90_copy_att( ncid(e), n, trim(name), outid, n ) )
        enddo
    enddo
    
    CHECK( nf90_enddef( outid ) )
    
    ! copy variables that don't depend on npoints
    CHECK( nf90_inquire_dimension( ncid(e), 1, len=dimlen ) )
    allocate( dvalues(dimlen) )
    CHECK( nf90_get_var( ncid(e), 1, dvalues ) )
    CHECK( nf90_put_var( outid, 1, dvalues ) )
    deallocate( dvalues )
    allocate( cvalues(dimlen) )
    CHECK( nf90_get_var( ncid(e), 2, cvalues ) )
    CHECK( nf90_put_var( outid, 2, cvalues ) )
    deallocate( cvalues )
    if ( numdims == 4 ) then
        allocate( ivalues(10) )
        CHECK( nf90_get_var( ncid(e), 6, ivalues ) )
        CHECK( nf90_put_var( outid, 6, ivalues ) )
        deallocate( ivalues )
    endif
    
    ! write npoints out to file
    allocate( ivalues(dimcount) )
    do n = 1, dimcount
        ivalues(n) = n
    enddo
    CHECK( nf90_put_var( outid, 3, ivalues ) )
    deallocate( ivalues )
    
    ! copy latitude and longitude variables
    allocate( fvalues(dimcount,1,1) )
    do v = 4, 5
        d = 1
        do f = 1, num
            if ( exists(f) ) then
                CHECK( nf90_inquire_dimension( ncid(f), 3, len=e ) )
                CHECK( nf90_get_var( ncid(f), v, fvalues(d:d+e-1,:,:) ) )
                d = d + e
            endif
        enddo
        CHECK( nf90_put_var( outid, v, fvalues ) )
    enddo
    deallocate( fvalues )
    
    ! copy the rest of the variables
    if ( numdims == 3 ) then
        e = 6
        allocate( fvalues(dimcount,dimlen,1) )
    else
        e = 7
        allocate( fvalues(dimcount,10,dimlen) )
    endif
    do v = e, numvars
        d = 1
        do f = 1, num
            if ( exists(f) ) then
                CHECK( nf90_inquire_dimension( ncid(f), 3, len=e ) )
                CHECK( nf90_get_var( ncid(f), v, fvalues(d:d+e-1,:,:) ) )
                d = d + e
            endif
        enddo
        CHECK( nf90_put_var( outid, v, fvalues ) )
    enddo

    CHECK( nf90_close( outid ) )

    deallocate( fvalues )
    deallocate( vardims )
    deallocate( dimid )
    
    
    

end subroutine pbpmerge

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine respfmerge( nprocs, filenames, outfile )

implicit none

! parameters
integer, intent(in) :: nprocs
character(len=256), dimension(nprocs), intent(in) :: filenames
character(len=256), intent(in) :: outfile

! local variables
integer :: nsib
integer :: nsoil
double precision, dimension(:,:), allocatable :: respfactor
double precision, dimension(:,:), allocatable :: resptotal
integer :: p, i, j
integer :: status
character(len=256) :: command

    ! copy file if only one process, then return out of subroutine
    if ( nprocs == 1 ) then
        command = "cp "//trim(filenames(1))//" "//trim(outfile)
        call system( trim(command) )
        return
    endif

    open( unit=3, file=trim(filenames(1)), form='unformatted',  &
        iostat=status ) !jk
    read( 3, iostat=status ) nsib
    read( 3, iostat=status ) nsoil
    close( 3 )

    allocate( respfactor(nsib, nsoil) )
    allocate( resptotal(nsib,nsoil) )
    
    resptotal(:,:) = 0.0
    do p = 1, nprocs
        open( unit=3, file=trim(filenames(p)), form='unformatted' ) !jk
        read( 3, iostat=status ) nsib
        read( 3, iostat=status ) nsoil
        read( 3, iostat=status ) respfactor(:,:)
        close( 3 )
        if ( status > 0 ) then
            stop 'Error reading respfactor'
        endif
        do j = 1, nsoil
            do i = 1, nsib
                resptotal(i,j) = resptotal(i,j) + respfactor(i,j)
            enddo
        enddo
    enddo
    
    open( unit=3, file=trim(outfile), form='unformatted' )
    write( 3 ) nsib
    write( 3 ) nsoil
    write( 3 ) resptotal(:,:)
    close( 3 )

    deallocate( respfactor )
    deallocate( resptotal )
    
end subroutine respfmerge













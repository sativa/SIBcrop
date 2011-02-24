program sibdrive

use kinds
use timetype
use sib_const_module
use sib_io_module
use sibtype
use sib_bc_module

implicit none

! local variables
type(sib_t), dimension(:), allocatable :: sib
type(time_struct) :: time
integer(kind=long_kind) :: t
integer(kind=int_kind) :: i,k
integer(kind=int_kind) :: rank
integer(kind=int_kind) :: nchunks
integer, external :: iargc
real(kind=dbl_kind) del_store
real(kind=dbl_kind) sum_flux
real(kind=dbl_kind) residual

!variables for timing
 real etime          ! Declare the type of etime()
 real elapsed(2)     ! For receiving user and system time
 real total          ! For receiving total time
      
! command line variables
character(len=4) :: one, two
character(len=4) :: dfdfd

!Crop variables:
type time_dep_var
   real (kind=real_kind) :: fPAR    ! Canopy absorbed fraction of PAR
   real (kind=real_kind) :: LAI     ! Leaf-area index
   real (kind=real_kind) :: Green   ! Canopy greeness fraction of LAI
   real (kind=real_kind) :: zo      ! Canopy roughness coeff 
   real (kind=real_kind) :: zp_disp ! Zero plane displacement
   real (kind=real_kind) :: RbC     ! RB Coefficient (c1)
   real (kind=real_kind) :: RdC     ! RC Coefficient (c2)
   real (kind=real_kind) :: gmudmu  ! Time-mean leaf projection
end type time_dep_var

type(time_dep_var) TimeVar

    ! read in parallelization values from command line
    call getarg( 1, one )
    if ( one == '' .or. one == '>' ) then
        rank = 1
        nchunks = 1
    else
        call getarg( 2, two )
        if ( two == '' .or. two == '>' ) then
            stop 'Command line arguments incorrect:  SiBD3 rank nchunks'
        else
            read( one, * ) rank
            read( two, * ) nchunks
            if ( rank > nchunks ) stop 'nchunks greater than rank'
            if ( rank < 1 .or. nchunks < 1 ) stop 'rank or nchunks < 1'
        endif
    endif

    ! read in namel_sibdrv
    call init_grid( rank, nchunks )
    
    ! allocate sib structure
    allocate( sib(subcount) )
    
    call init_var(sib)

    ! initialize all values and prepare for timestep loop
    call init_sibdrv( sib,time )

    ! set time varaibles for initial time step
    call time_check( time)
         
    ! call output_control
    call output_control( sib, time, rank )

    ! timestep loop
    do t = time%start_second, time%end_second - time%dtsib, time%dtsib

	! print out date information once a day
        if ( time%sec_day == time%dtsib ) then
          print*, time%month_names(time%month),time%day,time%year
        endif
	
        ! calculate solar declination
        if ( time%new_day ) call solar_dec( time )

        ! read in driver data needed
        if ( time%read_driver ) then
            if ( drvr_type == 'ecmwf' ) then
                call sibdrv_read_ecmwf( sib, time )
            elseif ( drvr_type == 'ncep1' ) then
                call sibdrv_read_ncep1( sib, time )
            elseif ( drvr_type == 'ncep2' ) then
                call sibdrv_read_ncep2( sib, time )
            elseif ( drvr_type == 'geos4' ) then
                call sibdrv_read_geos4( sib, time )
            elseif ( drvr_type == 'single' ) then
                call sibdrv_read_single( sib, time )
            !kdcorbin, 02/11
            elseif ( drvr_type == 'narr' ) then
                call sibdrv_read_narr( sib,time )
            else
                stop 'Invalid drvr_type specified'
            endif
            
            ! calculate mean cosine zenith angle
            call mean_zenith_angle( sib, time )
        endif

        ! read in bc data needed
        if ( time%read_bc ) then
             call update_bc( sib, time)
        endif

        ! update crop variables and interpolate bc variables once a day
        ! kdcorbin, 01/11 
        if ( time%sec_day == time%dtsib ) then
            call crop_accum(sib,time,timevar)
            call bc_interp( sib, time )
         endif

        ! interpolate
        call sibdrv_interp( sib, time )

        ! call sib_main()
        dtt = time%dtsib
        dti = 1./dtt
        tau = time%sec_year

        do i = 1, subcount
            call sib_main( sib(i),time )
        enddo

        ! call respfactor_control
        call respfactor_control( sib, time, rank )
	
	! call time_manager()
        call time_manager( time, drvr_type, roll_respf )
        sib(:)%stat%julday = time%doy 

        ! call output_control
        call output_control( sib, time, rank )
	
        ! write restart
        if ( time%write_restart ) then
            call rtape_sib( sib, time, rank )
        endif

    ! cas CO2 conservation
        do i = 1, subcount
          del_store=sib(i)%prog%cas-sib(i)%prog%cas_old
          sum_flux=sib(i)%diag%respg*dtt
          sum_flux=sum_flux-sib(i)%diag%assimn(6)*dtt
          sum_flux=sum_flux-sib(i)%diag%cflux*dtt
          residual=del_store-sum_flux-sib(i)%prog%expand

         !
         ! code output to test carbon conservation
         ! write(unit=85,11) t, residual,del_store, sum_flux,      &
         !      sib(i)%prog%expand,sib(i)%prog%cas, sib(i)%prog%cas_old,      &
         !      sib(i)%diag%respg*dtt, -sib(i)%diag%assimn(6)*dtt,            &
         !      -sib(i)%diag%cflux*dtt, sib(i)%diag%ra, sib(i)%prog%pco2ap
         !11     format(i14,11(2x,e15.8))
      enddo
	
    enddo


    ! make sure all files have been closed
    call file_closer
    !
    ! print message
    print*, 'end simulation'
      total = etime(elapsed)
    print *, 'End: total=', total, ' user=', elapsed(1),' system=', elapsed(2)

!itb_crop...close diagnostic files
   close(20)


end program sibdrive

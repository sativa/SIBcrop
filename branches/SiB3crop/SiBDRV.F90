program sibdrive

use kinds
use timetype
use sib_const_module
use sib_io_module
use sibtype
implicit none

! local variables
type(sib_t), dimension(:), allocatable :: sib
type(time_struct) :: time
integer(kind=long_kind) :: t
integer(kind=int_kind) :: i
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
    call init_sibdrv( sib, time )
    
    ! set time varaibles for initial time step
    call time_check( time)
    
      
    ! call output_control
    call output_control( sib, time, rank )
 
!
! test file to test CO2 conservation
!      open(unit=85, file='test', form='formatted')
    
    ! timestep loop
    do t = time%start_second, time%end_second - time%dtsib, time%dtsib
     
	! print out date information once a day
        if ( time%sec_day == time%dtsib ) then
          print*, time%month_names(time%month),time%day,time%year
        endif
	

        ! calculate solar declination
        if ( time%new_day ) call solar_dec( time )

!itb_crop...calculate Ta_bar (and any other crop stuff we 
!itb_crop...need )
        if ( time%new_day ) call crop_accum(sib)

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
            else
                stop 'Invalid drvr_type specified'
            endif
            
            ! calculate mean cosine zenith angle
            call mean_zenith_angle( sib, time )
        endif

        ! read in bc data needed
        if ( time%read_bc ) then
            call new_bc( sib, time )
        endif

        ! interpolate bc variables once a day      
        if ( time%new_day ) call bc_interp( sib, time )

        ! interpolate
        call sibdrv_interp( sib, time )
	
        ! call sib_main()
        dtt = time%dtsib
        dti = 1./dtt
        tau = time%sec_year
        !$OMP DO

        do i = 1, subcount
            call sib_main( sib(i) )
        enddo
        !$OMP END DO
        
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
        del_store=sib(1)%prog%cas-sib(1)%prog%cas_old
        sum_flux=sib(1)%diag%respg*dtt
        sum_flux=sum_flux-sib(1)%diag%assimn(6)*dtt
        sum_flux=sum_flux-sib(1)%diag%cflux*dtt
        residual=del_store-sum_flux-sib(1)%prog%expand
!
! code output to test carbon conservation
!        write(unit=85,11) t, residual,del_store, sum_flux, sib(1)%prog%expand,sib(1)%prog%cas, sib(1)%prog%cas_old, sib(1)%diag%respg*dtt, -sib(1)%diag%assimn(6)*dtt, -sib(1)%diag%cflux*dtt, sib(1)%diag%ra, sib(1)%prog%pco2ap
!11     format(i14,11(2x,e15.8))
	
    enddo


    ! make sure all files have been closed
    call file_closer
    !
    ! print message
    print*, 'end simulation'
      total = etime(elapsed)
    print *, 'End: total=', total, ' user=', elapsed(1),' system=', elapsed(2)

end program sibdrive

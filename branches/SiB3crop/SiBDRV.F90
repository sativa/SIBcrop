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


!itb_crop

! begin time dependant, output variables
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

type(aero_var),dimension(50,50) :: tempaerovar
real(kind=real_kind),dimension(2,2) :: temptran,tempref
integer(kind=int_kind) :: temp_biome

!itb_crop_end



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

        if ( time%new_day .AND. time%doy > time%init_doy)  then
  
          call crop_accum(sib,time,timevar)


        endif
!itb_crop_end



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

print*,'bc:',sib(1)%diag%phen_switch,sib(1)%param%zlt

!itb_crop...bypassing the call to new_bc; don't want to 
!itb_crop...use the NDVI-derived parameters
!EL.........reading in monthly-varying physfracs


            call read_physfrac( sib, time )
			

!itb_crop...what we want to do is use the min_ndvi_crop value
!itb_crop...set in sibtype, unless the phenology model is
!itb_crop...being utilized

    if(sib(1)%param%biome >= 20.0) temp_biome = 12

        do i = 1, subcount
             sib(i)%param%aparc1       = sib(i)%param%aparc2
             sib(i)%param%zlt1         = sib(i)%param%zlt2
             sib(i)%param%green1       = sib(i)%param%green2
             sib(i)%param%z0d1         = sib(i)%param%z0d2
             sib(i)%param%zp_disp1     = sib(i)%param%zp_disp2
             sib(i)%param%rbc1         = sib(i)%param%rbc2
             sib(i)%param%rdc1         = sib(i)%param%rdc2
             sib(i)%param%gmudmu1      = sib(i)%param%gmudmu2
             sib(i)%param%d13cresp1    = sib(i)%param%d13cresp2
             do k = 1, physmax
                 sib(i)%param%physfrac1(k) = sib(i)%param%physfrac2(k)
             enddo
         enddo

    temptran(1,1) = sib(1)%param%tran(1,1)
    temptran(1,2) = sib(1)%param%tran(1,2)
    temptran(2,1) = sib(1)%param%tran(2,1)
    temptran(2,2) = sib(1)%param%tran(2,2)

    tempref(1,1) = sib(1)%param%ref(1,1)
    tempref(1,2) = sib(1)%param%ref(1,2)
    tempref(2,1) = sib(1)%param%ref(2,1)
    tempref(2,2) = sib(1)%param%ref(2,2)

    tempaerovar = aerovar(:,:,temp_biome)
	
	If (sib(1)%diag%phen_switch==0) then


         call mapper(                              &
            latsib(1),                             &
            time%mid_month(time%pmonth),           &
            sib(1)%diag%min_ndvi_crop,             &
            sib(1)%diag%min_ndvi_crop,             &
            sib(1)%diag%min_fvcov_crop,            &
            sib(1)%param%chil,                     &
            temptran,                              &
            tempref,                               & 
            morphtab(temp_biome),                  &
            tempaerovar,                           &
            laigrid,                               &
            fvcovergrid,                           &
            timevar)

            sib(1)%param%aparc2 = timevar%fpar
            sib(1)%param%zlt2 = timevar%lai
            sib(1)%param%green2 = timevar%green
            sib(1)%param%z0d2 = timevar%zo
            sib(1)%param%zp_disp2 = timevar%zp_disp
            sib(1)%param%rbc2 = timevar%rbc
            sib(1)%param%rdc2 = timevar%rdc
            sib(1)%param%gmudmu2 = timevar%gmudmu
			

        endif

       endif

!itb_crop_end






        ! interpolate bc variables once a day      
        if ( time%new_day ) call bc_interp( sib, time )

        ! interpolate
        call sibdrv_interp( sib, time )



!itb_crop...
       ! sib(:)%diag%year = time%year
		!sib(:)%diag%doy=time%doy	!to calculate planting dates- EL
!itb_crop_end...
	
        ! call sib_main()
        dtt = time%dtsib
        dti = 1./dtt
        tau = time%sec_year
        !$OMP DO

        do i = 1, subcount

!print*,'LAI:',sib(i)%diag%phen_lai,sib(i)%diag%leafwt_c

            call sib_main( sib(i),time )
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

!itb_crop...close diagnostic files
   close(20)


end program sibdrive

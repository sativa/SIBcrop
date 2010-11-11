subroutine time_init( time )

use timetype
use sib_const_module
use sib_io_module
implicit none

! parameters
type(time_struct), intent(out) :: time

! local variables
integer :: x

    ! set some non-varying values
    time%days_per_year = 365
    time%sec_per_day = 86400
    time%mid_month(:) = (/15.5,45.,74.5,105.,135.5,166.,196.5,227.5,258.,288.5,319.,349.5/)
    time%days_per_month(:) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    time%doy1_month(:) = (/1,32,60,91,121,152,182,213,244,274,305,335/)
    time%month_names(:) = (/'January   ',  &
                            'February  ',  &
                            'March     ',  &
                            'April     ',  &
                            'May       ',  &
                            'June      ',  &
                            'July      ',  &
                            'August    ',  &
                            'September ',  &
                            'October   ',  &
                            'November  ',  &
                            'December  '/)

    ! if starttime > 0, then ntinital already in units of seconds
    ! if starttime < 0, convert days to seconds 
    if ( starttime < 0 ) starttime = (-starttime-1)*time%sec_per_day

    ! if endtime > 0, then endtime already in units of seconds
    ! if endtime < 0, convert days to seconds 
    if ( endtime < 0 ) endtime = (-endtime)*time%sec_per_day

    ! make sure endtime doesn't occur before starttime
    if ( endyear == startyear ) then
      if ((starttime >= endtime) )  &
          stop 'simulation ends before it starts, check starttime,'  &
              //' endtime, startyear, and endyear'

    endif
    
    ! make sure endyear doesn't occur before startyear
    if ( endyear < startyear )  &
        stop 'simulation ends before it starts, check startyear and endyear'
    
    ! make sure number of seconds in simulation is evenly divisible by dtsib
    if ( mod( endtime-starttime, dtsib) /= 0 )  &
        stop 'dtsib does not divide evenly into the total simulation'

    ! make sure dtsib divides evenly into a day
    if (  mod( time%sec_per_day, dtsib ) /= 0 )  &
        stop 'dtsib does not divide evenly into a day'
        
    ! make sure driver data timestep divides evenly into a day
    if ( mod( time%sec_per_day, dtsibmetin) /= 0 )  &
        stop 'dtsibmetin does not divide evenly into a day'
    
    ! make sure dtsibout is evenly divisible by dtsib
    if ( dtsibout > 0 .and. mod( dtsibout, dtsib ) /= 0 )  &
        stop 'dtsib does not divide evenly into dtsibout'
    
    ! make sure dtsibres is evenly divisible by dtsib
    if ( dtsibres > 0 .and. mod( dtsibres, dtsib ) /= 0 )  &
        stop 'dtsib does not divide evenly into dtsibres'
    
    ! make sure ndtsibpbp is positive
    if ( ndtsibpbp < 0 ) stop 'ndtsibpbp must be >= 0'
    
    
    ! set initial values
    time%init_year = startyear
    time%init_doy = starttime / time%sec_per_day + 1
    time%init_second = starttime
    do x = 1, 12
        if ( time%init_doy >= time%doy1_month(x) ) then
            time%init_month = x
            time%init_day = time%init_doy - time%doy1_month(x) + 1
        endif
    enddo

    
    !   else set them to init_* values
        time%start_year = time%init_year
        time%start_month = time%init_month
        time%start_day = time%init_day
        time%start_doy = time%init_doy
        time%start_second = time%init_second

    
    ! set ending values
    time%end_year = endyear
    ! end_second must be calculated incrementaly to avoid overflow that
    !  occurs when there are >= 58 years of simulation
    time%end_second = time%days_per_year * time%sec_per_day
    time%end_second = time%end_second * (time%end_year - time%start_year)
    time%end_second = endtime + time%end_second
    time%end_doy = endtime / time%sec_per_day + 1
    do x = 1, 12
        if ( time%end_doy >= time%doy1_month(x) ) then
            time%end_month = x
            time%end_day = time%end_doy - time%doy1_month(x) + 1
        endif
    enddo
    
    ! set values from namel to variables in structure
    time%dtsib = dtsib
    time%driver_step = dtsibmetin
    time%bc_step = dtsibbcin
    time%restart_step = dtsibres
    time%qp_step = dtsibout
    time%pbp_step = time%dtsib * ndtsibpbp
    
    ! set current times to get ready for simulation
    time%year = time%start_year
    time%month = time%start_month
    time%doy = time%start_doy
    time%day = time%start_day
    time%hour = time%sec_day/3600.
    time%sec_day = 0 
    time%sec_year = time%start_second
    time%sec_tot = time%start_second

   
    ! set previous times for boundary condition data
    time%pmonth = time%start_month - 1
    if ( time%pmonth == 0 ) time%pmonth = 12
    time%pyear = time%year
    if ( time%pmonth > time%start_month ) time%pyear = time%year - 1
    time%ppmonth = time%pmonth - 1
    if ( time%ppmonth == 0 ) time%ppmonth = 12
    time%ppyear = time%pyear
    if ( time%ppmonth > time%pmonth ) time%ppyear = time%pyear - 1
    
    ! set next times for boundary condition data
    time%nmonth = time%month
    time%nyear = time%year

    ! set pbp averaging variable
    time%pbp_incnt = 1
    time%qp_incnt = 1

    ! set some driver data timing variables
    time%driver_times = time%sec_per_day/time%driver_step
    time%driver_month = time%month
    time%driver_year = time%year
    time%driver_recnum = (time%day-1) * time%driver_times + 1
    time%driver_hour = int( time%hour, long_kind )
    time%driver_day = time%day
    time%start_period = (time%start_year-1) * time%days_per_year +  &
        real(time%start_doy-1) + real(time%start_second) / real(time%sec_per_day)
    
    ! initialize flags
    time%write_qp = .false.
    time%switch_qp = .true.
    time%pbp_offset = 0
    time%write_pbp = .false.
    time%pbp_count = 0    
    time%qp_count = 1
    time%switch_pbp = .true.
    time%read_driver = .true.
    time%read_bc = .true.
    time%switch_bc = .false.
    time%write_restart = .false.
    time%interp_bc = .true.
    time%interp_driver = .true.
    time%switch_driver = .true.
    time%new_day = .true.
    time%calc_respf = .false.
    time%write_respf = .false.
    



end subroutine time_init


!-------------------------------------------------------------------------------
subroutine time_check( time )
!-------------------------------------------------------------------------------

use timetype
use sib_io_module
implicit none

! parameters
type(time_struct), intent(inout) :: time 

! local variables
    

    time%switch_driver = .false.
    if (drvr_type == 'single' ) then
            time%driver_recnum = 1
    else
        time%driver_recnum =  time%driver_recnum + 1
    endif
    time%driver_hour = int( time%hour + time%driver_step/3600., long_kind )
    if ( time%driver_hour > 24 )  &
        time%driver_day = mod(time%driver_day,time%days_per_month(time%month))+1
    time%driver_hour = mod( time%driver_hour, 24 )

   if ( time%driver_recnum >  &
            time%days_per_month(time%month) * time%driver_times ) then
         
        time%driver_recnum = 1
        time%driver_month = time%nmonth
        time%driver_year = time%nyear   
        time%switch_driver = .true.
    endif
    

end subroutine time_check

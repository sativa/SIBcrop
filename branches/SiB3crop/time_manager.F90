subroutine time_manager( time, drvr_type, roll_respf )
!----------------------------------------------------------------------
! calculates all time related variables and sets all operations control flags
!
! Modifications:
!  Kevin Schaefer shifted read driver back 1 time step (8/17/04)
!  Kevin Schaefer shifted switch driver data back 1 time step (8/19/04)
!  Kevin Schaefer corrected driver month/year calculation (8/19/04)

use kinds
use timetype
implicit none

! parameters
type(time_struct), intent(inout) :: time
character (len=8), intent(in) :: drvr_type 
logical(kind=log_kind), intent(in) :: roll_respf

! local variables
integer :: x
integer :: low_bound
integer :: high_bound


    !---------------------------------------------------------------------------
    ! TIME VARIABLES
    !---------------------------------------------------------------------------

    ! update time variables
    time%sec_tot = time%sec_tot + time%dtsib
    time%sec_year = time%sec_year + time%dtsib
    time%sec_day = time%sec_day+time%dtsib
    if(mod(time%sec_year,time%sec_per_day)==600 .and. &
       time%sec_year>time%init_second+600)then
    	   time%doy=time%doy+1
	   time%sec_day =600
    endif
    do x = 1, 12
        if ( time%doy >= time%doy1_month(x) .and. time%sec_day > 0) then
            time%month = x
            time%day = time%doy - time%doy1_month(x) + 1
        endif
    enddo
    time%hour = mod( time%sec_day, time%sec_per_day ) / 3600.
    if( time%hour == 0)time%hour =24.0
    if (time%sec_year == time%dtsib+(time%sec_per_day*time%days_per_year))  then
       time%year = time%year + 1
       time%month = 1
       time%doy = 1
       time%day = 1
       time%sec_day = time%dtsib
       time%sec_year = time%dtsib
       
    endif
    !---------------------------------------------------------------------------
    ! DRIVER DATA
    !---------------------------------------------------------------------------
    
    ! check to see if it is time to read Driver Data
    if ( mod( time%sec_tot, time%driver_step ) == 0  &
        .and. time%sec_tot > time%init_second + time%dtsib ) then
        time%read_driver = .true.
        if ( drvr_type == 'single' ) then
            time%driver_recnum = 1
        else
            time%driver_recnum =  time%driver_recnum + 1
        endif
        time%driver_hour = int( time%hour + time%driver_step/3600. )
        if ( time%driver_hour >= 24 ) then
            time%driver_day =  &
                mod( time%driver_day, time%days_per_month(time%month) ) + 1
        endif
        time%driver_hour = mod( time%driver_hour, 24 )
    else
        time%read_driver = .false.
    endif
    
    ! check to see if it is time to switch Driver Data files
    if ( time%day == time%days_per_month(time%month) .and.  &
        time%sec_day == time%sec_per_day - time%driver_step ) then
        
        time%switch_driver = .true.
        time%driver_recnum = 1
        time%driver_month = time%driver_month+1
	!
	! KS check if at end of year
	if(time%driver_month>12) then
	   time%driver_month=1
	   time%driver_year=time%driver_year+1
	endif
    else
        time%switch_driver = .false.
    endif
    

    !---------------------------------------------------------------------------
    ! BOUNDARY CONDITIONS
    !---------------------------------------------------------------------------

    ! check to see if it is time to read boundary condition data
    low_bound = (time%mid_month(time%month)-1) * time%sec_per_day - time%dtsib
    high_bound = (time%mid_month(time%month)-1) * time%sec_per_day + time%dtsib
    if ( time%sec_year > low_bound .and. time%sec_year < high_bound ) then
        time%read_bc = .true.
    else
        time%read_bc = .false.
    endif
    
    ! check to see if it is time to switch boundary condition data files
    if ( time%month == size(time%mid_month) .and. time%read_bc ) then
        time%switch_bc = .true.
    else
        time%switch_bc = .false.
    endif
    
    ! if we are going to read in bc data, prepare for next month
    if ( time%read_bc ) then
        time%nmonth = time%month + 1
        if ( time%nmonth > 12 ) then
            time%nmonth = 1
            time%nyear = time%year + 1
        endif
    endif

    !---------------------------------------------------------------------------
    ! RESTART FILES
    !---------------------------------------------------------------------------

    ! check to see if it is time to write restart file
    if ( time%restart_step > 0 ) then
        ! time%restart_step units = seconds
        if ( mod( time%sec_year, time%restart_step ) == 0 ) then
            time%write_restart = .true.
        else
            time%write_restart = .false.
        endif
    else
        ! time%restart_step units = months
        if (time%sec_year == (time%doy1_month(time%nmonth)-1)*time%sec_per_day .or. &
        (time%sec_year == (365)*time%sec_per_day)) then
            time%write_restart = .true.
        else
            time%write_restart = .false.
        endif
    endif
    

    !---------------------------------------------------------------------------
    ! QP FILES
    !---------------------------------------------------------------------------

    ! check to see if it is time to write data to qp files
    time%qp_count = time%qp_count + 1
     if ( time%qp_step > 0 ) then
         ! time%qp_step units = seconds
         if ( mod( time%sec_year, time%qp_step ) == 0 ) then
             time%write_qp = .true.
             time%qp_incnt = 1. / real(time%qp_count)
             time%qp_count = 0
             time%end_period = real((time%year-1) * time%days_per_year) +  &
                 real(time%sec_year)/real(time%sec_per_day)
             time%period_length = time%end_period - time%start_period
             time%start_period = time%end_period
         else
             time%write_qp = .false.
         endif
     else
        ! time%qp_step units = months
        if ( time%sec_year==(time%doy1_month(time%nmonth)-1)*86400 .or. &
	    time%sec_year == time%sec_per_day*time%days_per_year) then
            
            time%write_qp = .true.
            time%qp_incnt = 1. /  real(time%qp_count)
            time%qp_count = 0
            time%end_period = real((time%year-1) * time%days_per_year) +  &
                real(time%sec_year+time%dtsib)/real(time%sec_per_day)
            time%period_length = time%end_period - time%start_period
            time%start_period = time%end_period
        else
            time%write_qp = .false.
        endif
     endif
    
    ! check to see if it is time to switch qp files
    if ( time%day == 1 .and. time%sec_day == time%dtsib) then
        time%switch_qp = .true.
    else
        time%switch_qp = .false.
    endif


    !---------------------------------------------------------------------------
    ! PBP FILES
    !---------------------------------------------------------------------------

    ! check to see if it is time to write data to pbp files
    if ( time%pbp_step > 0 ) then

        ! set pbp averaging variable
        time%pbp_count = time%pbp_count + 1

        if ( mod( time%sec_tot, time%pbp_step ) == time%pbp_offset ) then
            time%write_pbp = .true.
            time%pbp_incnt = 1. / real(time%pbp_count)
            time%pbp_count = 0
        else
            time%write_pbp = .false.
        endif
        
        if ( time%day == 1 .and. time%sec_day == time%dtsib) then
            time%switch_pbp = .true.
        else
            time%switch_pbp = .false.
        endif
    else
        ! do not write pbp data at all
        time%write_pbp = .false.
        time%switch_pbp = .false.
    endif
    
    
    !---------------------------------------------------------------------------
    ! INTERPOLATIONS
    !---------------------------------------------------------------------------
    
    if ( time%sec_day == time%dtsib ) then
        time%new_day = .true.
    else
        time%new_day = .false.
    endif
    
    
    !---------------------------------------------------------------------------
    ! RESPFACTOR
    !---------------------------------------------------------------------------
    if ( roll_respf ) then
        ! rolling respfactor
        !   calculate last timestep of the month
        if ( time%day == time%days_per_month(time%month) .and.   &
             time%sec_day + time%dtsib == time%sec_per_day ) then
        
            time%calc_respf = .true.
            time%write_respf = .true.
        else
            time%calc_respf = .false.
            time%write_respf = .false.
        endif
    else
        ! traditional respfactor
        !   calculate last timestep of a complete year
        if ( time%month == 12 .and.   &
             time%doy == 365  .and.   &
             time%sec_day + time%dtsib == time%sec_per_day .and.   &
             time%sec_year + time%dtsib ==   &
                time%sec_per_day * time%days_per_year ) then
             
            time%calc_respf = .true.
            time%write_respf = .true.
        else
            time%calc_respf = .false.
            time%write_respf = .false.
        endif
    endif
    

end subroutine time_manager

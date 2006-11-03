module timetype

!-------------------------------------------------------------------------------
! Author:  Owen Leonard
! Date:    April 17, 2004
! Purpose:
!   This modules contains a user defined type containing all of the time 
! variables assigned values in time_check, and updated by time_manager.
!-------------------------------------------------------------------------------

use kinds
implicit none

public time_struct

type time_struct
    
    
    ! time constants
    integer(kind=long_kind) :: init_year     ! initial year of simulation
    integer(kind=long_kind) :: init_month    ! initial month of simulation
    integer(kind=long_kind) :: init_day      ! initial day of month of simulation
    integer(kind=long_kind) :: init_doy      ! initial day of year of simulation
    integer(kind=long_kind) :: init_second  ! initial second of simulation

    integer(kind=long_kind) :: start_year    ! year of restart file
    integer(kind=long_kind) :: start_month   ! month of restart file
    integer(kind=long_kind) :: start_day     ! day of month of restart file
    integer(kind=long_kind) :: start_doy     ! day of year of restart file
    integer(kind=long_kind) :: start_second ! second of restart file
    
    integer(kind=long_kind) :: end_year      ! last year of simulation
    integer(kind=long_kind) :: end_month     ! last month of simulation
    integer(kind=long_kind) :: end_day       ! last day of month of simulation
    integer(kind=long_kind) :: end_doy       ! last day of year of simulation
    integer(kind=long_kind) :: end_second   ! last second of simulation
    
    integer(kind=long_kind) :: total_years   ! total number of years in simulation
    integer(kind=long_kind) :: total_months  ! total number of months in simulation
    integer(kind=long_kind) :: total_days    ! total number of days in simulation

    integer(kind=long_kind) :: dtsib         ! # seconds in simulation time step
    integer(kind=long_kind) :: bc_int_step   ! # seconds between bc interpolation
    integer(kind=long_kind) :: drvr_int_step ! # seconds between driver data interpolation
    
    real(kind=real_kind),   dimension(12) :: mid_month      ! mid-month day of year
    integer(kind=long_kind), dimension(12) :: days_per_month ! # days per month
    integer(kind=long_kind), dimension(12) :: doy1_month     ! day of year of first 
                                                            !   day in month
    character(len=10), dimension(12) :: month_names  ! names of the months
    
    integer(kind=long_kind) :: sec_per_day   ! # seconds in a day 
    integer(kind=long_kind) :: days_per_year ! number of days in the current year
    
    
    ! time variables
    integer(kind=long_kind) :: year          ! current year in the simulation
    integer(kind=int_kind) :: month         ! current month in the simulation
    real(kind=real_kind)   :: hour          ! current hour of day in the simulation
    integer(kind=long_kind) :: day           ! current day of the current month 
    integer(kind=long_kind) :: doy           ! current day of current year
    integer(kind=long_kind) :: sec_day       ! current second in the current day
    integer(kind=long_kind) :: sec_year      ! current second in the current year
    integer(kind=long_kind) :: sec_tot      ! current second in the whole simulation
    
    integer(kind=long_kind) :: pyear         ! year of previous month of simulation
    integer(kind=long_kind) :: pmonth        ! previous month of simulation
    integer(kind=long_kind) :: ppyear        ! year of two months previous of simulation
    integer(kind=long_kind) :: ppmonth       ! two months previous of simultion
    
    integer(kind=long_kind) :: nyear         ! year of next month of simulation
    integer(kind=long_kind) :: nmonth        ! next month of simulation
    
    integer(kind=long_kind) :: driver_times   ! # driver data timesteps per day
    integer(kind=int_kind) :: driver_recnum ! record # of data in driver data file
    integer(kind=long_kind) :: driver_month  ! month of driver data to read
    integer(kind=long_kind) :: driver_year   ! year of driver data to read
    integer(kind=int_kind) :: driver_hour   ! hour of driver data to be read
    integer(kind=long_kind) :: driver_day    ! day of driver data to be read
    
    real(kind=dbl_kind) :: start_period     ! start of averaged period
    real(kind=dbl_kind) :: end_period       ! end of averaged period
    real(kind=dbl_kind) :: period_length    ! length of averaged period
    
    ! I/O time variables
    integer(kind=long_kind) :: driver_step   ! # seconds in driver data timestep
    integer(kind=long_kind) :: bc_step       ! # seconds between bc file switching
    integer(kind=long_kind) :: restart_step  ! # seconds in restart file timestep
    integer(kind=long_kind) :: qp_step       ! # seconds in qp file timestep
    real(kind=dbl_kind)    :: qp_incnt      ! # to multiply qp's by to get average
    integer(kind=long_kind) :: qp_count      ! # to divide qp's by to get average
    integer(kind=long_kind) :: pbp_step      ! # seconds in pbp file timestep
    real(kind=dbl_kind)   :: pbp_incnt     ! # to multiply pbp's by to get average
    integer(kind=long_kind) :: pbp_count     ! # to divide pbp's by to get average
    integer(kind=long_kind) :: pbp_offset    ! # seconds to offset pbp writing
    
    
    ! Flags
    logical(kind=log_kind) :: write_qp      ! write data to qp files ?
    logical(kind=log_kind) :: switch_qp     ! switch qp files ?
    logical(kind=log_kind) :: write_pbp     ! write data to pbp files ?
    logical(kind=log_kind) :: switch_pbp    ! switch pbp files ?
    logical(kind=log_kind) :: read_driver   ! read driver data ?
    logical(kind=log_kind) :: switch_driver ! switch driver data file ?
    logical(kind=log_kind) :: read_bc       ! read boundary condition data ?
    logical(kind=log_kind) :: switch_bc     ! switch bc file ?
    logical(kind=log_kind) :: write_restart ! write restart file ?
    logical(kind=log_kind) :: interp_bc     ! interpolate boundary condition data ?
    logical(kind=log_kind) :: interp_driver ! interpolate driver data ?
    logical(kind=log_kind) :: new_day       ! new day ?  (solar dec and bc interp)
    logical(kind=log_kind) :: calc_respf    ! calculate new respfactor
    logical(kind=log_kind) :: write_respf   ! write respfactor file

end type time_struct

end module timetype

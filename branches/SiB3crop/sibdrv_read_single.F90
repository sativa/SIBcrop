subroutine sibdrv_read_single( sib, time )
!
! Modifications:
!  Kevin Schaefer moved conversion from pascals to millibars from sibdrv_interp to here (8/16/04)

use kinds
use sibtype
use timetype
use sib_const_module
use sib_io_module
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib ! NOTE: subcount should = 1
type(time_struct), intent(in) :: time

! local variables
integer(kind=int_kind) :: i
integer(kind=int_kind) :: status
real(kind=dbl_kind) :: yr, doy, hr,dy
real(kind=dbl_kind) :: temp_dpt  ! dew point
character(len=256) :: filename
character(len=13) :: subname
character(len=1025) :: record
data subname/'sibdrv_read '/

   !*** Storing previous time steps data
    do i=1,subcount
        sib(i)%prog%ps1       = sib(i)%prog%ps2
        sib(i)%prog%tm1       = sib(i)%prog%tm2
        sib(i)%prog%tcc1      = sib(i)%prog%tcc2
        sib(i)%prog%sh1       = sib(i)%prog%sh2
        sib(i)%prog%spdm1     = sib(i)%prog%spdm2
        sib(i)%prog%lspr1     = sib(i)%prog%lspr2
        sib(i)%prog%cupr1     = sib(i)%prog%cupr2
        sib(i)%prog%dlwbot1   = sib(i)%prog%dlwbot2
        sib(i)%prog%sw_dwn1   = sib(i)%prog%sw_dwn2
    enddo

    ! switch files if needed
    if ( time%switch_driver ) then
        close( 87, iostat = status )
               
        write(unit=filename,fmt=dr_format)time%driver_year, time%driver_month
        open( unit=87, file=trim(filename), form='formatted', iostat=status)
        if ( status > 0 ) then
            print *, 'SiBDRV_read_single'
            print *, 'Error opening file'
            stop
        endif
    endif
    
    ! read one line of driver data from file
!    print *,'SiBDRV_init_std'
!    print *,'opening drive files for ',time%sec_day
!    print *, 'dr_form=',trim(dr_format)
 !  print *, 'filename=',trim(filename)
    do i = 1, time%driver_recnum
        do  ! Read until not a comment.
            read( 87,'(a)', iostat=status ) record
            if ( status > 0 ) then
                print *, 'SiBDRV_read_single'
                print *, 'Error reading file'
                stop
            endif
!print*,record
!pause
            if ( record(1:1) .ne. '#' ) exit
        enddo

        	read(unit=record,fmt=*)yr,doy,hr,sib(1)%prog%tm2,temp_dpt, &
            sib(1)%prog%spdm2,sib(1)%prog%ps2,sib(1)%prog%dlwbot2,   &
            sib(1)%prog%sw_dwn2,sib(1)%prog%lspr2,sib(1)%prog%cupr2
!
! calculate large scale precipitation
        sib(1)%prog%lspr2 = sib(1)%prog%lspr2 !- sib(1)%prog%cupr2 
!print *,sib(1)%prog%lspr2
!EL the latter part of the above line was commented out since the driver data&
!for bondville had large scale ppt separately.
!
! KS comvert from pascals to millibars
        !sib(1)%prog%ps2=sib(1)%prog%ps2*0.01
!
! KS convert dew point to specific humidity
    	!  call qsat_eau(1,sib%prog%ps2*100.0,temp_dpt,sib%prog%sh2) 
			sib%prog%sh2=temp_dpt
! check for zero humidity
    if(sib(1)%prog%sh2==0.) sib(1)%prog%sh2=1.e-4
  
  enddo

end subroutine sibdrv_read_single

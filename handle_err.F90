subroutine handle_err( status, routine, number )
!==========================================================================
! handle_err is a generic subroutine that checks for errors in the netcdf 
! functions and prints out a description of the error.  it then terminates 
! the program
!
#ifdef PGF
use netcdf
use typeSizes
#endif

! parameters
integer, intent(in) :: status                       ! error status
character(len=*), intent(in), optional :: routine   ! routine where err occurred
integer, intent(in), optional :: number             ! error "line" number
                                                    ! used to pinpoint which line
                                                    ! caused the error

    print *, 'error', status, (nf90_strerror(status))
    if ( present(routine) ) print *, routine
    if ( present(number) ) print *, 'at',number
    stop 'stopped-netcdf error'

end subroutine handle_err

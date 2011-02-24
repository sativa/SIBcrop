! Subroutine to check netcdf routines
! kdcorbin, 02/11
 
subroutine check(status)

use netcdf

implicit none

integer, intent ( in ) :: status
 
     if (status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Error with netcdf.  Stopped."
      end if
end subroutine check   

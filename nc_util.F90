! check the return status of a netcdf function and stop and report the error
! if present.
! You should probably use the CHECK macro to automatically fill in file and
! line values.
subroutine nc_check(status, file, line)
  use netcdf

  implicit none

  integer, intent (in) :: status, line
  character(len=*), intent(in) :: file

  if (status /= nf90_noerr) then
     print "(A,' line ',I4,' ',A)", file, line, (nf90_strerror(status))
     stop
  end if
end subroutine nc_check

! Looks for varname in the netcdf file ncid.  If the variable
! does not exist, an error is printed and the program is terminated.
! You should probably use the ENSURE_VAR macro to automatically fill in
! file and line values.
subroutine nc_ensure_var(ncid, varname, varid, file, line)
  use netcdf

  implicit none

  integer, intent(in) :: ncid, line
  integer, intent(out) :: varid
  character(len=*), intent(in) :: varname, file
  integer :: status

  status = nf90_inq_varid(ncid, varname, varid)

  if (status /= nf90_noerr) then
     print "(A, ' line ',I4,' in nf90_inq_varid(',A,') ',A)", file, line, &
          varname, (nf90_strerror(status))
     stop
  end if
end subroutine nc_ensure_var

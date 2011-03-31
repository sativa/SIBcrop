! check_var_exists looks for varname in the netcdf file ncid.  If the variable
! does not exist, an error is printed and the program is terminated.
subroutine ensure_var(ncid, varname, varid, file, line)
  use netcdf

  integer, intent(in) :: ncid, line
  integer, intent(out) :: varid
  character(len=*), intent(in) :: varname, file
  integer :: status

  status = nf90_inq_varid(ncid, varname, varid)

  if (status /= nf90_noerr) then
     print *, file, ':', line, ': netcdf error', status, &
          (nf90_strerror(status)), 'reading variable', varname
     stop
  end if
end subroutine ensure_var

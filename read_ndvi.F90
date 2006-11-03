subroutine read_ndvi(filename, month, ndvidata, physfrac, d13cresp)
!--------------------------------------------------------------------
! opens current sib_bc file and checks if nsib matches, then reads data

use kinds
use sib_const_module, only: nsib,     & ! total number of SiB points
                            subset,   & ! array mapping subdomain to nsib vector
                            subcount, & ! actual number of SiB points being simulated
                            physmax    ! maximum # of physiology types
use sib_io_module, only: drvr_type,  & !jk flag to indicate driver data type
                         c4_source,  &
                         ndvi_source,  &
                         d13cresp_source
#ifdef PGF
use netcdf
use typeSizes
#endif

! define input variables
character *100, intent(in) :: filename ! Filename for sib_bc
integer(kind=int_kind), intent(in) :: month ! month to read
real(kind=real_kind), dimension(subcount), intent(inout) :: ndvidata ! a month of ndvi data
real(kind=dbl_kind), dimension(subcount,physmax), intent(out) :: physfrac
real(kind=dbl_kind), dimension(subcount), intent(inout) :: d13cresp

! define local variables
integer(kind=int_kind) :: i,j
integer(kind=int_kind) :: ntest1     ! value of nsib read from sib_bc file
integer(kind=int_kind) :: dimid
integer(kind=int_kind) :: status
integer(kind=int_kind) :: begin (2)  ! indices where to start and 
integer(kind=int_kind) :: finish (2) !   finish reading from file
integer(kind=int_kind) :: yyid       ! ndvi file id#
integer(kind=int_kind) :: ndvi_id    ! ndvi variable id#
integer(kind=int_kind) :: phys_id    ! physfrac variable id#
integer(kind=int_kind) :: d13_id     ! d13cresp variable id#
integer(kind=int_kind) :: x,y
real(kind=real_kind), dimension(nsib) :: ndvi ! temp variables for parameters
real(kind=real_kind), dimension(nsib,physmax)  :: frac
real(kind=real_kind), dimension(nsib)  :: d13
character(len=10) :: name


    ! Open the sib_bc file
    status = nf90_open ( trim(filename), nf90_nowrite, yyid )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_dimid ( yyid, 'nsib', dimid )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension ( yyid, dimid, name, ntest1 )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid ( yyid, 'ndvi', ndvi_id )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid ( yyid, 'd13cresp', d13_id )
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inq_varid ( yyid, 'physfrac', phys_id )
    if (status /= nf90_noerr) call handle_err(status)

    ! test if boundary condition files are correct
    if(ntest1 /= nsib) stop ' open: file sib_bc no match with model for nsib'

    ! read in one month's data for all nsib points
    begin (1)  = month 
    begin (2)  = 1 
    finish (1) = 1 
    finish (2) = nsib 
    status = nf90_get_var ( yyid, ndvi_id, ndvi, begin, finish )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( yyid, d13_id, d13 )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( yyid, phys_id, frac )
    if (status /= nf90_noerr) call handle_err (status)
    
    ! read in global attributes that are passed on to output files
    status = nf90_get_att( yyid, nf90_global, 'ndvi_source', ndvi_source )
    status = nf90_get_att( yyid, nf90_global, 'c4_source', c4_source )
    status = nf90_get_att( yyid, nf90_global, 'd13cresp_source',  &
        d13cresp_source )

    ! close file
    status = nf90_close ( yyid )
    if (status /= nf90_noerr) call handle_err (status)

    ! copy only points in subdomain
    do i = 1, subcount
        ndvidata(i) = ndvi(subset(i))
        d13cresp(i) = d13(subset(i))
        do j = 1, physmax
            physfrac(i,j) = frac(subset(i),j)
        enddo
    enddo

end subroutine read_ndvi

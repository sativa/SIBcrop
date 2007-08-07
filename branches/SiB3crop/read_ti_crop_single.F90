!=============================================================
subroutine read_ti_crop_single (sib)
!=============================================================
! Opens sib_bc_TI.nc and reads in the boundary condition data.
!     REFERENCES:
!     CREATED:
!       Owen Leonard   August 3, 2001
!     MODIFICATIONS:
!   Owen Leonard  added code to read tables   10/15/01
!   Owen Leonard  added check on zwind and ztemp  10/24/01
!   Owen Leonard  modified code to be f90 freeform  4/25/02
!   Owen Leonard  modified code to be compatible w/SiB3  4/25/02
!     SUBROUTINES CALLED:
!     FUNCTIONS CALLED:
!       netcdf library
!

use kinds
use sibtype
use sib_const_module, only:  &
    zwind, &
    ztemp, &
    nsib,  &
    subcount, &
    subset
use sib_io_module, only:  &
    param_path,  &
    biome_source,  &
    soil_source,  &
    soref_source
use sib_bc_module
#ifdef PGF
use netcdf
use typeSizes
#endif


! declare input variables
type(sib_t) :: sib(subcount)

!itb...output variable


! declare local variables
integer(kind=int_kind) ::  i,q,j           !  index variables for loops
integer(kind=int_kind) ::  start (2)      !  arrays that define where
integer(kind=int_kind) ::  done (2)       !  to start and stop reading 
integer(kind=int_kind) ::  begin (3)      !  variables from 
integer(kind=int_kind) ::  finish (3)     !  sib_bc_TI.nc
integer(kind=int_kind) ::  tiid           !  file id#
integer(kind=int_kind) ::  status         !  error status check
integer(kind=int_kind) ::  ndims          !  number of dimensions (soiltype)
integer(kind=int_kind) ::  dimid, nvar    !  dimension id#, number of variables
integer(kind=int_kind) ::  phys           !  phys dimension value
!  used for AeroVar structure
integer(kind=int_kind) ::  varid          ! varible ID
integer(kind=int_kind) ::  numsoil, &       ! soil type #s
    biovar, &
    svar,   &
    morphvar
integer(kind=int_kind) :: soilnum(nsib)
real(kind=real_kind), dimension(nsib) :: biome
real(kind=real_kind), dimension(nsib) :: clayfrac
real(kind=real_kind), dimension(nsib) :: sandfrac
real(kind=real_kind), dimension(nsib) :: vcover
real(kind=real_kind), dimension(nsib) :: soref1
real(kind=real_kind), dimension(nsib) :: soref2
real(kind=real_kind) :: fclay
real(kind=real_kind) :: fsand

real(kind=real_kind), dimension(:,:), allocatable :: biovart3   ! 2D arrays for the tables
real(kind=real_kind), dimension(:,:), allocatable :: biovart4   !   "   "
real(kind=real_kind), dimension(:,:), allocatable :: soilvart   !   "   "
real(kind=real_kind), dimension(:,:), allocatable :: morphvart !    "   "
integer(kind=int_kind), dimension(:,:), allocatable :: phystype

real(kind=dbl_kind) testzwind   !  used to test if zwind and ztemp in TI file
real(kind=dbl_kind) testztemp   !  match the values defined in namel file
character(len=10) name

    ! open sib_bc_TI.nc file 
    status = nf90_open ( trim(param_path)//'_TI.nc', nf90_nowrite, tiid )
    if (status /= nf90_noerr) call handle_err (status)
    print *, '\t reading sib_bc_TI.nc'

    ! check zwind and ztemp values
    status = nf90_inq_varid (tiid, 'zwind', varid)
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, testzwind )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inq_varid ( tiid, 'ztemp', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, testztemp )
    if ( status /= nf90_noerr ) call handle_err (status)
    if ( zwind /= testzwind ) then
        print *, '\t\t zwind value in ti file does not match namel file'
        print *, '\t\t ti:  ', testzwind, 'namel:  ', zwind
        !      stop
    endif
    if ( ztemp /= testztemp ) then
        print *, '\t\t ztemp value in ti file does not match namel file'
        print *, '\t\t ti:  ', testztemp, 'namel:  ', ztemp
        !      stop
    endif


    ! read in variables
    status = nf90_inq_dimid ( tiid, 'numvar', dimid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_dimension ( tiid, dimid, name,nvar )
    if (status /= nf90_noerr) call handle_err (status)

    ! biome
    status = nf90_inq_varid ( tiid, 'biome', varid )
    if (status /= nf90_noerr) call handle_err (status)
!    status = nf90_get_var ( tiid, varid, biome )
!    if (status /= nf90_noerr) call handle_err (status)

    ! soilnum
    status = nf90_inq_varid ( tiid, 'soilnum', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_variable ( tiid, varid, name,ndims )
    if (status /= nf90_noerr) call handle_err (status)
    start (1) = 1
    start (2) = 1
    done (1)  = nsib
    done (2)  = 1


    ! laigrid
    status = nf90_inq_varid ( tiid, 'laigrid', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, laigrid )
    if (status /= nf90_noerr) call handle_err (status)

    ! fvcgrid
    status = nf90_inq_varid ( tiid, 'fvcgrid', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, fvcovergrid )
    if (status /= nf90_noerr) call handle_err (status)

    ! aero_zo
    allocate (aerovar(50,50,nvar))
    status = nf90_inq_varid ( tiid, 'aero_zo', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, aerovar%zo )
    if (status /= nf90_noerr) call handle_err (status)

    ! aero_zp
    status = nf90_inq_varid ( tiid, 'aero_zp', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, aerovar%zp_disp )
    if (status /= nf90_noerr) call handle_err (status)

    ! aero_rbc
    status = nf90_inq_varid ( tiid, 'aero_rbc', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, aerovar%rbc )
    if (status /= nf90_noerr) call handle_err (status)

    ! areo_rdc
    status = nf90_inq_varid ( tiid, 'aero_rdc', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, aerovar%rdc )
    if (status /= nf90_noerr) call handle_err (status)

    !---------------read in the tables----------------------------------

    status = nf90_inq_dimid ( tiid, 'numsoil', dimid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_dimension ( tiid, dimid, name, numsoil )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inq_dimid ( tiid, 'biovar', dimid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_dimension ( tiid, dimid, name,biovar )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inq_dimid ( tiid, 'soilvar', dimid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_dimension ( tiid, dimid, name,svar )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inq_dimid ( tiid, 'morphvar', dimid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_dimension ( tiid, dimid,name,morphvar )
    if (status /= nf90_noerr) call handle_err (status)

    allocate(morphtab(nvar))
    allocate(biovart3(nvar,biovar))
    allocate(biovart4(nvar,biovar))
    allocate(soilvart(numsoil,svar))
    allocate(morphvart(nvar,morphvar))

    status = nf90_inq_varid ( tiid, 'biome_tablec3', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, biovart3 )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inq_varid ( tiid, 'biome_tablec4', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, biovart4 )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inq_varid ( tiid, 'soil_table', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, soilvart )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inq_varid ( tiid, 'morph_table', varid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_get_var ( tiid, varid, morphvart )
    if (status /= nf90_noerr) call handle_err (status)

    ! read in global attributes that are passed on to output files
    status = nf90_get_att( tiid, nf90_global, 'biome_source', biome_source )
    status = nf90_get_att( tiid, nf90_global, 'soil_source', soil_source )
    status = nf90_get_att( tiid, nf90_global, 'soref_source', soref_source )

    status = nf90_close(tiid)
    
    !------------Assign values from tables into data structures--------


    ! Morph Table used in mapper
    do q = 1, nvar
        morphtab(q)%zc = morphvart (q,1)
        morphtab(q)%lwidth = morphvart (q,2)
        morphtab(q)%llength = morphvart (q,3)
        morphtab(q)%laimax = morphvart (q,4)
        morphtab(q)%stems = morphvart (q,5)
        morphtab(q)%ndvimax = morphvart (q,6)
        morphtab(q)%ndvimin = morphvart (q,7)
        morphtab(q)%srmax = morphvart (q,8)
        morphtab(q)%srmin = morphvart (q,9)
    enddo

    !---------deallocate 2D arrays---------------------------

    deallocate(biovart3)
    deallocate(biovart4)
    deallocate(soilvart)
    deallocate(morphvart)


end subroutine read_ti_crop_single


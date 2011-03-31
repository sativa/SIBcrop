!====================================================
subroutine read_ti (sib)
!====================================================
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
    subset, &
    !kdcorbin, 02/11
    phys,  &
    phystype,biovart3,biovart4, &
    soilvart,morphvart
use sib_io_module, only:  &
    param_path,  &
    biome_source,  &
    soil_source,  &
    soref_source, &
    drvr_type  !kdcorbin, 01/11
use sib_bc_module

use netcdf
use typeSizes

#include "nc_util.h"

! declare input variables
type(sib_t) :: sib(subcount)

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

!  used for AeroVar structure
integer(kind=int_kind) ::  varid          ! varible ID
integer(kind=int_kind) ::  numsoil, &       ! soil type #s
    biovar, svar, morphvar
integer(kind=int_kind) :: soilnum(nsib)
real(kind=real_kind), dimension(nsib) :: biome
real(kind=real_kind), dimension(nsib) :: clayfrac
real(kind=real_kind), dimension(nsib) :: sandfrac
real(kind=real_kind), dimension(nsib) :: vcover
real(kind=real_kind), dimension(nsib) :: soref1
real(kind=real_kind), dimension(nsib) :: soref2
real(kind=real_kind) :: fclay
real(kind=real_kind) :: fsand

real(kind=dbl_kind) testzwind   !  used to test if zwind and ztemp in TI file
real(kind=dbl_kind) testztemp   !  match the values defined in namel file
character(len=10) name

    ! open sib_bc_TI.nc file 
    print *, '      reading ',trim(param_path)//'_TI.nc'
    status = nf90_open ( trim(param_path)//'_TI.nc', nf90_nowrite, tiid )
    if (status /= nf90_noerr) call handle_err (status)

    ! check zwind and ztemp values
    ENSURE_VAR(tiid, 'zwind', varid)
    status = nf90_get_var ( tiid, varid, testzwind )
    if (status /= nf90_noerr) call handle_err (status)
    ENSURE_VAR( tiid, 'ztemp', varid )
    status = nf90_get_var ( tiid, varid, testztemp )
    if ( status /= nf90_noerr ) call handle_err (status)
    if ( zwind /= testzwind ) then
        print *, '     zwind value in ti file does not match namel file'
        print *, '     ti:  ', testzwind, 'namel:  ', zwind
        !      stop
    endif
    if ( ztemp /= testztemp ) then
        print *, '     ztemp value in ti file does not match namel file'
        print *, '     ti:  ', testztemp, 'namel:  ', ztemp
        !      stop
    endif

    ! read in variables
    status = nf90_inq_dimid ( tiid, 'numvar', dimid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_dimension ( tiid, dimid, name,nvar )
    if (status /= nf90_noerr) call handle_err (status)

    ! biome - not read in for single point runs
    ! kdcorbin, 01/11
   if (drvr_type .ne. 'single') then
       ENSURE_VAR( tiid, 'biome', varid )
       status = nf90_get_var ( tiid, varid, biome )
       if (status /= nf90_noerr) call handle_err (status)
   endif

    ! soilnum
    ENSURE_VAR( tiid, 'soilnum', varid )
    if (status /= nf90_noerr) call handle_err (status)
    !kdcorbin, 02/11 - added name,ndims specifiers
    status = nf90_inquire_variable ( tiid, varid, name=name,ndims=ndims )
    if (status /= nf90_noerr) call handle_err (status)

    start (1) = 1
    start (2) = 1
    done (1)  = nsib
    done (2)  = 1
    if (ndims.eq.1) then
        status = nf90_get_var ( tiid, varid, soilnum )
        if (status /= nf90_noerr) call handle_err (status)
    else 
        status = nf90_get_var ( tiid, varid, clayfrac, start, done )
        if (status /= nf90_noerr) call handle_err (status) 
        start (2) = 2
        status = nf90_get_var ( tiid, varid, sandfrac, start, done )
        if (status /= nf90_noerr) call handle_err (status)
    endif

    ! fvcover
    ENSURE_VAR( tiid, 'fvcover', varid )
    status = nf90_get_var ( tiid, varid, vcover )
    if (status /= nf90_noerr) call handle_err (status)

    ! phystype
    status = nf90_inq_dimid ( tiid, 'phys', dimid )
    if (status /= nf90_noerr) call handle_err (status)
    status = nf90_inquire_dimension ( tiid, dimid, name, phys )
    if (status /= nf90_noerr) call handle_err (status)
    allocate(phystype(nsib,phys))
    ENSURE_VAR( tiid, 'phystype', varid )
    status = nf90_get_var ( tiid, varid, phystype )
    if (status /= nf90_noerr) call handle_err (status)

    ! sorefvis
    ENSURE_VAR( tiid, 'sorefvis', varid )
    status = nf90_get_var ( tiid, varid, soref1 )
    if (status /= nf90_noerr) call handle_err (status)

    ! sorefnir
    ENSURE_VAR( tiid, 'sorefnir', varid )
    status = nf90_get_var ( tiid, varid, soref2 )
    if (status /= nf90_noerr) call handle_err (status)

    ! laigrid
    ENSURE_VAR( tiid, 'laigrid', varid )
    status = nf90_get_var ( tiid, varid, laigrid )
    if (status /= nf90_noerr) call handle_err (status)

    ! fvcgrid
    ENSURE_VAR( tiid, 'fvcgrid', varid )
    status = nf90_get_var ( tiid, varid, fvcovergrid )
    if (status /= nf90_noerr) call handle_err (status)

    ! aero_zo
    allocate (aerovar(50,50,nvar))
    ENSURE_VAR( tiid, 'aero_zo', varid )
    status = nf90_get_var ( tiid, varid, aerovar%zo )
    if (status /= nf90_noerr) call handle_err (status)

    ! aero_zp
    ENSURE_VAR( tiid, 'aero_zp', varid )
    status = nf90_get_var ( tiid, varid, aerovar%zp_disp )
    if (status /= nf90_noerr) call handle_err (status)

    ! aero_rbc
    ENSURE_VAR( tiid, 'aero_rbc', varid )
    status = nf90_get_var ( tiid, varid, aerovar%rbc )
    if (status /= nf90_noerr) call handle_err (status)

    ! areo_rdc
    ENSURE_VAR( tiid, 'aero_rdc', varid )
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

    ENSURE_VAR( tiid, 'biome_tablec3', varid )
    status = nf90_get_var ( tiid, varid, biovart3 )
    if (status /= nf90_noerr) call handle_err (status)
    ENSURE_VAR( tiid, 'biome_tablec4', varid )
    status = nf90_get_var ( tiid, varid, biovart4 )
    if (status /= nf90_noerr) call handle_err (status)
    ENSURE_VAR( tiid, 'soil_table', varid )
    status = nf90_get_var ( tiid, varid, soilvart )
    if (status /= nf90_noerr) call handle_err (status)
    ENSURE_VAR( tiid, 'morph_table', varid )
    status = nf90_get_var ( tiid, varid, morphvart )
    if (status /= nf90_noerr) call handle_err (status)

    ! read in global attributes that are passed on to output files
    status = nf90_get_att( tiid, nf90_global, 'biome_source', biome_source )
    status = nf90_get_att( tiid, nf90_global, 'soil_source', soil_source )
    status = nf90_get_att( tiid, nf90_global, 'soref_source', soref_source )

    status = nf90_close(tiid)
    
    !------------Assign values from tables into data structures--------

    !Do not assign values for single point runs - kdcorbin, 01/11
    if (drvr_type .ne. 'single') then

      do i = 1, subcount

        sib(i)%param%biome = biome(subset(i))
        sib(i)%param%soref(1) = soref1(subset(i))
        sib(i)%param%soref(2) = soref2(subset(i))
        sib(i)%param%vcover = vcover(subset(i))
        if ( ndims /= 1 ) then
            sib(i)%param%clayfrac = clayfrac(subset(i))
            sib(i)%param%sandfrac = sandfrac(subset(i))
        endif
        do j=1,phys  ! phystype loop
            sib(i)%param%phystype(j) = phystype(subset(i),j)
        enddo

        call set_ti(sib(i))

        ! Soil Variables Table
        if (ndims == 1) then
            sib(i)%param%bee   = soilvart(int(soilnum(subset(i))),1)
            sib(i)%param%phsat = soilvart(int(soilnum(subset(i))),2)
            sib(i)%param%satco = soilvart(int(soilnum(subset(i))),3)
            sib(i)%param%poros = soilvart(int(soilnum(subset(i))),4)
            sib(i)%param%slope = soilvart(int(soilnum(subset(i))),5)
            sib(i)%param%wopt  = soilvart(int(soilnum(subset(i))),6)
            sib(i)%param%zm    = soilvart(int(soilnum(subset(i))),7)
            sib(i)%param%wsat  = soilvart(int(soilnum(subset(i))),8)

            ! calculate %sand and $clay from table based on the approximate
            !   centroid of soil texture category within the UDSA texture triangle
            !   centroid % clay/sand estimated by Kevin Schaefer  (3/30/01)
            select case (soilnum(subset(i)))
            case (1)
                sib(i)%param%clayfrac = 3.
                sib(i)%param%sandfrac = 92.
            case (2)
                sib(i)%param%clayfrac = 5.
                sib(i)%param%sandfrac = 82.
            case (3)
                sib(i)%param%clayfrac = 10.
                sib(i)%param%sandfrac = 65.
            case (4)
                sib(i)%param%clayfrac = 13.
                sib(i)%param%sandfrac = 22.
            case (5)
                sib(i)%param%clayfrac = 7.
                sib(i)%param%sandfrac = 7.
            case (6)
                sib(i)%param%clayfrac = 18.
                sib(i)%param%sandfrac = 42.
            case (7)
                sib(i)%param%clayfrac = 28.
                sib(i)%param%sandfrac = 58.
            case (8)
                sib(i)%param%clayfrac = 40.
                sib(i)%param%sandfrac = 52.
            case (9)
                sib(i)%param%clayfrac = 39.
                sib(i)%param%sandfrac = 32.
            case (10)
                sib(i)%param%clayfrac = 39.
                sib(i)%param%sandfrac = 10.
            case (11)
                sib(i)%param%clayfrac = 41.
                sib(i)%param%sandfrac = 7.
            case (12)
                sib(i)%param%clayfrac = 65.
                sib(i)%param%sandfrac = 19.
            case default
                print *, 'Illegal value for soil type:  ', soilnum(subset(i))
                print *, 'nsib point:  ', (subset(i))
                sib(i)%param%clayfrac = 33.
                sib(i)%param%sandfrac = 33.
            end select

        else  ! ndims
            ! calculate soil variables from %sand, %clay
            ! equations taken from SoilProperties subroutine in mapper.F
            ! resp. variable curve fits from Raich et al., 1991 

            ! patch!!!!!!!!!  %sand, %clay maps have bad data points
            ! values of 3.3961514E+38 :: patch changes them to 33 (neutral ground)
            ! approx. 50 bad data points
            if (sib(i)%param%clayfrac > 100)  sib(i)%param%clayfrac = 33.
            if (sib(i)%param%sandfrac > 100)  sib(i)%param%sandfrac = 33.

            fclay = sib(i)%param%clayfrac/100.
            fsand = sib(i)%param%sandfrac/100.
            sib(i)%param%bee = 2.91+0.159*sib(i)%param%clayfrac
            sib(i)%param%phsat = -10.*10**(1.88-0.0131*sib(i)%param%sandfrac)/ &
                1000.
            sib(i)%param%satco = 0.0070556*10**(-0.884+0.0153* &
                sib(i)%param%sandfrac)/1000.
            sib(i)%param%poros = 0.489-0.00126*sib(i)%param%sandfrac
            sib(i)%param%slope = 0.176
            sib(i)%param%wopt = (-0.08*fclay**2+0.22*fclay+0.59)*100.
            sib(i)%param%zm = -2*fclay**3-0.4491*fclay**2+0.2101*fclay+0.3478
            sib(i)%param%wsat = 0.25*fclay+0.5
        endif  ! ndims

    enddo  ! subcount loop
   endif !non-single statements only

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

    !deallocate(biovart3)
    !deallocate(biovart4)
    !deallocate(soilvart)
    !deallocate(morphvart)
    !deallocate(phystype)

end subroutine read_ti


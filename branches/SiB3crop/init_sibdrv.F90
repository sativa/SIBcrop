!---------------------------------------------------------------------
subroutine init_sibdrv( sib, time )
!---------------------------------------------------------------------

use sibtype
use timetype
use sib_const_module
use sib_io_module

!itb_crop...
use sib_bc_module
!itb_crop_end...


implicit none

!---------------------------------------------------------------------
!itb...init_sibdrv reads most initialization information. takes 
!itb...place of several BUGS routines, most notably init_global.
!
!     REFERENCES:
!
! Modifications:
!  - added dtsibbcin, dtsibmetin for possible different intervals
!    of reading in the sibbc and sibdrv met data dd, jk 980209
!  Kevin Schaefer moved read IC and respfactor to after driver data (8/12/04)
!  Kevin Schaefer added calls to read NCEP1 driver data (8/13/04)
!
!     SUBROUTINES CALLED:
!          none
!     FUNCTIONS CALLED:
!          none
!
!     INCLUDED COMMONS:
!
!     ARGUMENT LIST VARIABLES
!---------------------------------------------------------------------

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(inout) :: time

! local variables


integer(kind=int_kind) :: i

! begin time dependant, output variables
type time_dep_var
   real (kind=real_kind) :: fPAR    ! Canopy absorbed fraction of PAR
   real (kind=real_kind) :: LAI     ! Leaf-area index
   real (kind=real_kind) :: Green   ! Canopy greeness fraction of LAI
   real (kind=real_kind) :: zo      ! Canopy roughness coeff 
   real (kind=real_kind) :: zp_disp ! Zero plane displacement
   real (kind=real_kind) :: RbC     ! RB Coefficient (c1)
   real (kind=real_kind) :: RdC     ! RC Coefficient (c2)
   real (kind=real_kind) :: gmudmu  ! Time-mean leaf projection
end type time_dep_var

type(time_dep_var) TimeVar

!itb_crop

type(aero_var),dimension(50,50) :: tempaerovar
real(kind=real_kind),dimension(2,2) :: temptran,tempref


    print *, 'INIT_SIBDRV:'

    !itb---------------------------------------------------------------------
    !itb...initialize some seasonal diagnostic values...

    do i = 1, subcount
        sib(i)%diag%snow_end(1) = 365.0
        sib(i)%diag%snow_end(2) = 365.0
        sib(i)%diag%snow_end(3) = 365.0
        sib(i)%diag%tot_an(:)   = 0.0_dbl_kind
        sib(i)%diag%tot_ss(:,:)   = 0.0_dbl_kind
        sib(i)%stat%pt_num=i
    enddo    



    ! parse sib_qpopts and sib_pbpopts to see which variables are output
    call read_qp_pbp_opts
    
    ! initialize time variables
    print *, '\t initialize time variables'
    call time_init( time )
    sib(:)%stat%julday = time%doy

    ! read in time-invariant boundary conditions for global runs
    if ( drvr_type /= 'single' ) call read_ti(sib)
    
    ! calculate previous month's time-variant boundary conditions
    !  and read in time-invariant boundary conditions
    print *, '\t obtaining previous month time-variant boundary conditions'


!itb_crop...bypass the parameter file time-varying values. will
!itb_crop...call mapper from here, using minimum NDVI value
!itb_crop...specified in sibtype.F90


!itb_crop...still need to read in the time-invariant values.
    call previous_bc( sib, time )

!itb_crop...now replace the time-varying params
    print*,'previous_bc'


    call read_ti_crop_single(sib)

    tempaerovar = aerovar(:,:,int(sib(1)%param%biome))

!itb_crop...assign values to tempaerovar

    temptran(1,1) = sib(1)%param%tran(1,1)
    temptran(1,2) = sib(1)%param%tran(1,2)
    temptran(2,1) = sib(1)%param%tran(2,1)
    temptran(2,2) = sib(1)%param%tran(2,2)

    tempref(1,1) = sib(1)%param%ref(1,1)
    tempref(1,2) = sib(1)%param%ref(1,2)
    tempref(2,1) = sib(1)%param%ref(2,1)
    tempref(2,2) = sib(1)%param%ref(2,2)


    print*,'read_ti_crop'
    call mapper(          &
            latsib(1),        &
            time%mid_month(time%pmonth),           &
            sib(1)%diag%min_ndvi_crop,      &
            sib(1)%diag%min_ndvi_crop,      &
            sib(1)%diag%min_fvcov_crop, &
            sib(1)%param%chil,   &
            temptran,         &
            tempref,          &
            morphtab(int(sib(1)%param%biome)),      &
            tempaerovar,      &
            laigrid,          &
            fvcovergrid,      &
            timevar)



        sib(1)%param%aparc2 = timevar%fpar
        sib(1)%param%zlt2 = timevar%lai
        sib(1)%param%green2 = timevar%green
        sib(1)%param%z0d2 = timevar%zo
        sib(1)%param%zp_disp2 = timevar%zp_disp
        sib(1)%param%rbc2 = timevar%rbc
        sib(1)%param%rdc2 = timevar%rdc
        sib(1)%param%gmudmu2 = timevar%gmudmu

!itb_crop_end...



    call soil_properties( sib )
    

    ! read in initial driver data
    print *, '\t reading in initial time-step driver data'
    if ( drvr_type == 'ecmwf' ) then
        call sibdrv_read_ecmwf( sib, time )
    elseif ( drvr_type == 'ncep1' ) then
        call sibdrv_read_ncep1( sib, time )
    elseif ( drvr_type == 'ncep2' ) then
        call sibdrv_read_ncep2( sib, time )
    elseif ( drvr_type == 'geos4' ) then
        call sibdrv_read_geos4( sib, time )
    elseif ( drvr_type == 'single' ) then
        call sibdrv_read_single( sib, time )
    else
        stop 'Invalid drvr_type specified'
    endif
!
! read in initial conditions
    print *, '\t reading in initial conditions '
    call read_ic(sib,time)

! read in respfactor file
    print *, '\t read in respFactor'
    call read_respfactor(sib)

! calculate initial solar declination
    call init_solar_dec( time )

end subroutine init_sibdrv
!
!===============================================================================
subroutine read_qp_pbp_opts
!===============================================================================

use sib_io_module
use sib_const_module
implicit none

integer(kind=int_kind) :: i,n
logical(kind=log_kind) :: doqptem
integer(kind=int_kind) :: ipbp, ldummy, ndummy
character (len=16) :: nametem
character (len=80) :: listtem
integer(kind=int_kind), dimension(:), allocatable :: imulttem, imulttem2


    !---------------------------------------------------------------------------
    ! read sib_qpopts and count number of variables to be output
    !---------------------------------------------------------------------------
    open(unit=2,file=qp_path,form='formatted') !jk
    nqpsib = 0
    nqp3sib = 0
    do 
        read(2,*, end=922)doqptem,ldummy,nametem,ndummy,listtem
        if(ldummy.eq.1) then
            nqp3sib = nqp3sib + 1
        else if (ldummy.eq.0) then
            nqpsib = nqpsib + 1
        endif
    enddo

    922  continue

    rewind 2
    allocate (doqp3sib(nqp3sib))
    allocate (nameqp3sib(nqp3sib))
    allocate (listqp3sib(nqp3sib))
    allocate (numqp3sib(nqp3sib))
    allocate (doqpsib(nqpsib))
    allocate (nameqpsib(nqpsib))
    allocate (listqpsib(nqpsib))
    allocate (numqpsib(nqpsib))
    iiqp3sib = 0
    iiqpsib = 0
    do i = 1,nqp3sib+nqpsib
        read(2,*)doqptem,ldummy,nametem,ndummy,listtem
        if(ldummy.eq.1) then
            iiqp3sib = iiqp3sib + 1
            doqp3sib(iiqp3sib) = doqptem
            nameqp3sib(iiqp3sib) = nametem
            listqp3sib(iiqp3sib) = listtem
            numqp3sib(iiqp3sib) = ndummy
        else if (ldummy.eq.0) then
            iiqpsib = iiqpsib + 1
            doqpsib(iiqpsib) = doqptem
            nameqpsib(iiqpsib) = nametem
            listqpsib(iiqpsib) = listtem
            numqpsib(iiqpsib) = ndummy
        endif
    enddo 
    close(2)
    allocate (indxqp3sib(nqp3sib))
    allocate (indxqpsib(nqpsib))

    iiqpsib = 0
    do n = 1,nqpsib
        if(doqpsib(n)) then
            iiqpsib = iiqpsib + 1
            indxqpsib(n) = iiqpsib
        endif
    enddo
    iiqp3sib = 0
    do n = 1,nqp3sib
        if(doqp3sib(n)) then
            iiqp3sib = iiqp3sib + 1
            indxqp3sib(n) = iiqp3sib
        endif
    enddo
    do n = 1,nqpsib
        if(.not.doqpsib(n)) then
            indxqpsib(n) = iiqpsib + 1
        endif
    enddo
    do n = 1,nqp3sib
        if(.not.doqp3sib(n)) then
            indxqp3sib(n) = iiqp3sib + 1
        endif
    enddo


    !      initialize diagnostics         
    allocate (qpsib(subcount,iiqpsib+1))   
    allocate( qp2varid(nqpsib) )
    allocate( qp3varid(nqp3sib) )
    allocate (qp3sib(subcount,nsoil,iiqp3sib+1))   
    qp3sib(:,:,:) = 0.0
    qpsib(:,:) = 0.0
    print*,'\t diagnostics initialized'



    !---------------------------------------------------------------------------
    ! read sib_pbpopts and count number of variables to be output
    !---------------------------------------------------------------------------
    open(unit=2,file=pbp_path,form='formatted')   !jk
    npbpsib = 0
    npbp2sib = 0
    
    ! count number of variables listed for pbp and pbp2 data in sib_pbpopts
    do 
        read(2,*, end=932)doqptem,ldummy,nametem,ndummy,listtem
        if(ldummy.eq.1) then
            npbp2sib = npbp2sib + 1
        else if (ldummy.eq.0) then
            npbpsib = npbpsib + 1
        endif
    enddo 
    932  continue
    rewind 2

    allocate (dopbp2sib(npbp2sib))
    allocate (namepbp2sib(npbp2sib))
    allocate (listpbp2sib(npbp2sib))
    allocate (numpbp2sib(npbp2sib))
    allocate (dopbpsib(npbpsib))
    allocate (namepbpsib(npbpsib))
    allocate (listpbpsib(npbpsib))
    allocate (numpbpsib(npbpsib))

    ! count number of variables that are set to be saved to pbp files
    iipbp2sib = 0
    iipbpsib = 0
    do i = 1,npbp2sib+npbpsib
        read(2,*)doqptem,ldummy,nametem,ndummy,listtem
        if(ldummy.eq.1) then
            iipbp2sib = iipbp2sib + 1
            dopbp2sib(iipbp2sib) = doqptem
            namepbp2sib(iipbp2sib) = nametem
            listpbp2sib(iipbp2sib) = listtem
            numpbp2sib(iipbp2sib) = ndummy
        else if (ldummy.eq.0) then
            iipbpsib = iipbpsib + 1
            dopbpsib(iipbpsib) = doqptem
            namepbpsib(iipbpsib) = nametem
            listpbpsib(iipbpsib) = listtem
            numpbpsib(iipbpsib) = ndummy
        endif
    enddo 
    close(2)
    
    allocate (indxpbp2sib(npbp2sib))
    allocate (indxpbpsib(npbpsib))

    iipbpsib = 0
    do n = 1,npbpsib
        if(dopbpsib(n)) then
            iipbpsib = iipbpsib + 1
            indxpbpsib(n) = iipbpsib
        endif
    enddo
    iipbp2sib = 0
    do n = 1,npbp2sib
        if(dopbp2sib(n)) then
            iipbp2sib = iipbp2sib + 1
            indxpbp2sib(n) = iipbp2sib
        endif
    enddo
    do n = 1,npbpsib
        if(.not.dopbpsib(n)) then
            indxpbpsib(n) = iipbpsib + 1
        endif
    enddo
    do n = 1,npbp2sib
        if(.not.dopbp2sib(n)) then
            indxpbp2sib(n) = iipbp2sib + 1
        endif
    enddo

    allocate( pbpsib(iipbpsib+1,ijtlensib) )
    allocate( pbpvarid(npbpsib) )
    allocate( pbp2sib(nsoil,iipbp2sib+1,ijtlensib) )
    allocate( pbp2varid(npbp2sib) )
    pbpsib(:,:) = 0.0
    pbp2sib(:,:,:) = 0.0

end subroutine read_qp_pbp_opts
!
!===============================================================================
subroutine read_ic(sib,time)
!===============================================================================
!  Author:  Ian Baker
!  Modified by:  Owen Leonard
!  Date :  March 30, 2004
!  Purpose:
!    This subroutine reads in the initial conditions file and pulls out
!  only those points in the subdomain
!
! Modifications:
!  Kevin Schaefer moved soil layer calculations to soil_properties (10/27/04)
!===============================================================================

#ifdef PGF
use netcdf 
use typeSizes
#endif
use kinds
use sibtype
use timetype
use sib_const_module
use sib_io_module

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(inout) :: time
! netcdf variables
integer(kind=int_kind) :: ierr
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: varid

! local variables
integer(kind=int_kind) :: i,j,k
integer(kind=int_kind) :: nsibt
integer(kind=int_kind) :: nsoilt
integer(kind=int_kind) :: nsnowt
real(kind=int_kind) :: versiont
integer(kind=int_kind) :: subcountt
integer(kind=int_kind), dimension(2) :: start
integer(kind=int_kind), dimension(2) :: finish
real(kind=dbl_kind), dimension(nsib) :: ta
real(kind=dbl_kind), dimension(nsib) :: tc
!EL...crop vars (kind=int_kind)
integer(kind=int_kind), dimension(nsib) :: ndf_opt,pd,pd7,pd7_est,pdindx7
!EL...crop vars (kind=int_kind) end..
integer(kind=int_kind), dimension(nsib) :: nsl
real(kind=dbl_kind), dimension(nsib) :: pco2ap
real(kind=dbl_kind), dimension(nsib) :: d13cca
real(kind=dbl_kind), dimension(nsib) :: snow_veg
real(kind=dbl_kind), dimension(nsib) :: snow_age
real(kind=dbl_kind), dimension(nsib) :: snow_depth
real(kind=dbl_kind), dimension(nsib) :: snow_mass
real(kind=dbl_kind), dimension(nsib) :: tke
real(kind=dbl_kind), dimension(nsib) :: sha
real(kind=dbl_kind), dimension(nsib) :: capac1
real(kind=dbl_kind), dimension(nsib) :: capac2
real(kind=dbl_kind), dimension(nsib) :: coszbar
real(kind=dbl_kind), dimension(nsib) :: dayflag
!EL...crop vars..real(kind=dbl_kind)
real(kind=dbl_kind), dimension(nsib) :: tempf
real(kind=dbl_kind), dimension(nsib) :: gdd
real(kind=dbl_kind), dimension(nsib) :: w_main
real(kind=dbl_kind), dimension(nsib,4) :: cum_wt_prev
real(kind=dbl_kind), dimension(nsib,4) :: cum_drywt_prev
!EL...crop vars..real(kind=dbl_kind) end..
real(kind=dbl_kind), dimension(12,nsib) :: tot_an
real(kind=dbl_kind), dimension(nsib,6) :: rst
real(kind=dbl_kind), dimension(nsib,-nsnow+1:nsoil) :: deept
real(kind=dbl_kind), dimension(nsib,-nsnow+1:nsoil) :: www_liq
real(kind=dbl_kind), dimension(nsib,-nsnow+1:nsoil) :: www_ice
real(kind=dbl_kind), dimension(nsib,nsnow) :: nz_snow
real(kind=dbl_kind), dimension(nsib,nsnow) :: lz_snow
real(kind=dbl_kind), dimension(nsib,nsnow) :: dz_snow
real(kind=dbl_kind), dimension(12,nsib,nsoil) :: tot_ss

integer(kind=int_kind),dimension(11) :: map_totals
integer(kind=int_kind)               :: jday

DATA map_totals/31,59,90,120,151,181,212,243,273,304,334/

    print*,'ic_path=',trim(ic_path)

    ! read in initial conditions (restart file)
    ierr = nf90_open( trim(ic_path), nf90_nowrite, ncid )
    if( ierr /= nf90_noerr ) call handle_err(ierr)

    print *,'\t opened ic file', trim(ic_path)

    !itb...read some scalars
    ierr = nf90_inq_varid( ncid, 'nsib', varid )
    ierr = nf90_get_var( ncid, varid, nsibt )
    print *, '\t nsib=',nsib, ' total nsib=',nsibt
    if(nsib /= nsibt) stop'INITIAL CONDITIONS: NSIB INCORRECT'

    ierr = nf90_inq_varid( ncid, 'nsoil', varid )
    ierr = nf90_get_var( ncid, varid, nsoilt )
    if(nsoil /= nsoilt) stop'INITIAL CONDITIONS: NSOIL INCORRECT'

    ierr = nf90_inq_varid( ncid, 'nsnow', varid )
    ierr = nf90_get_var( ncid, varid, nsnowt )
    if(nsnow /= nsnowt) stop'INITIAL CONDITIONS: NSNOW INCORRECT'

!    ierr = nf90_inq_varid( ncid, 'subcount', varid )
!    ierr = nf90_get_var( ncid, varid, subcountt )
!    if(subcount /= subcountt) stop'INITIAL CONDITIONS: SUBCOUNT INCORRECT'

    ierr = nf90_inq_varid( ncid, 'version', varid )
    ierr = nf90_get_var( ncid, varid, versiont )

    ierr = nf90_inq_varid( ncid, 'nsecond', varid )
    ierr = nf90_get_var( ncid, varid, nsecond )
    if(nsecond /= time%sec_year) stop 'NSECONDS DOES NOT MATCH STARTTIME'

    !itb...read nsib vectors

    ierr = nf90_inq_varid( ncid, 'ta', varid )
    ierr = nf90_get_var( ncid, varid, ta )

    ierr = nf90_inq_varid( ncid, 'tc', varid )
    ierr = nf90_get_var( ncid, varid, tc )

    ierr = nf90_inq_varid( ncid, 'nsl', varid )
    ierr = nf90_get_var( ncid, varid, nsl )

    ierr = nf90_inq_varid( ncid, 'pco2a', varid )
    ierr = nf90_get_var( ncid, varid, pco2ap )

    ierr = nf90_inq_varid( ncid, 'd13cca', varid )
    ierr = nf90_get_var( ncid, varid, d13cca )

    ierr = nf90_inq_varid( ncid, 'snow_veg', varid )
    ierr = nf90_get_var( ncid, varid, snow_veg )

    ierr = nf90_inq_varid( ncid, 'snow_age', varid )
    ierr = nf90_get_var( ncid, varid, snow_age )

    ierr = nf90_inq_varid( ncid, 'snow_depth', varid )
    ierr = nf90_get_var( ncid, varid, snow_depth )

    ierr = nf90_inq_varid( ncid, 'snow_mass', varid )
    ierr = nf90_get_var( ncid, varid, snow_mass )

    ierr = nf90_inq_varid( ncid, 'tke', varid )
    ierr = nf90_get_var( ncid, varid, tke )

    ierr = nf90_inq_varid( ncid, 'sha', varid )
    ierr = nf90_get_var( ncid, varid, sha )

    !itb...read some 2-d vars
    ierr = nf90_inq_varid( ncid, 'td', varid )
    ierr = nf90_get_var( ncid, varid, deept )

    ierr = nf90_inq_varid( ncid, 'www_liq', varid )
    ierr = nf90_get_var( ncid, varid, www_liq )

    ierr = nf90_inq_varid( ncid, 'www_ice', varid )
    ierr = nf90_get_var( ncid, varid, www_ice )

    !itb...now the rest...
    ierr = nf90_inq_varid( ncid, 'capac1', varid )
    ierr = nf90_get_var( ncid, varid, capac1 )

    ierr = nf90_inq_varid( ncid, 'capac2', varid )
    ierr = nf90_get_var( ncid, varid, capac2 )

    ierr = nf90_inq_varid( ncid, 'coszbar', varid )
    ierr = nf90_get_var( ncid, varid, coszbar )

    ierr = nf90_inq_varid( ncid, 'dayflag', varid )
    ierr = nf90_get_var( ncid, varid, dayflag )

    ierr = nf90_inq_varid( ncid, 'rst', varid )
    ierr = nf90_get_var( ncid, varid, rst )

    !EL...crop variables..

    ierr = nf90_inq_varid( ncid, 'pd', varid )
    ierr = nf90_get_var( ncid, varid, pd )

print*,'init: PD=',pd

    ierr = nf90_inq_varid( ncid, 'pd7', varid )
    ierr = nf90_get_var( ncid, varid, gdd )

    ierr = nf90_inq_varid( ncid, 'pd7_est', varid )
    ierr = nf90_get_var( ncid, varid, w_main )

    ierr = nf90_inq_varid( ncid, 'pdindx7', varid )
    ierr = nf90_get_var( ncid, varid,cum_wt_prev )

    ierr = nf90_inq_varid( ncid, 'ndf_opt', varid )
    ierr = nf90_get_var( ncid, varid, cum_wt_prev )
    
    ierr = nf90_inq_varid( ncid, 'tempf', varid )
    ierr = nf90_get_var( ncid, varid, tempf )

    ierr = nf90_inq_varid( ncid, 'gdd', varid )
    ierr = nf90_get_var( ncid, varid, gdd )

    ierr = nf90_inq_varid( ncid, 'w_main', varid )
    ierr = nf90_get_var( ncid, varid, w_main )

    ierr = nf90_inq_varid( ncid, 'cum_wt_prev', varid )
    ierr = nf90_get_var( ncid, varid,cum_wt_prev )

    ierr = nf90_inq_varid( ncid, 'cum_drywt_prev', varid )
    ierr = nf90_get_var( ncid, varid, cum_wt_prev )

    !El.. end crop vars

    print*,'\t\t read in slabs...'

    !itb...don't know how to read slabs directly into the structure yet...
    start(1) = 1
    start(2) = 1
    finish(1) = nsib
    finish(2) = nsnow
    ierr = nf90_inq_varid( ncid, 'dzsnow', varid )
    ierr = nf90_get_var( ncid, varid, dz_snow, start, finish )

    ierr = nf90_inq_varid( ncid, 'lzsnow', varid )
    ierr = nf90_get_var( ncid, varid, lz_snow, start, finish )

    ierr = nf90_inq_varid( ncid, 'nzsnow', varid )
    ierr = nf90_get_var( ncid, varid, nz_snow, start, finish )

    ! read in tot_an and tot_ss for rolling respfactor
    ierr = nf90_inq_varid( ncid, 'tot_an', varid )
    ierr = nf90_get_var( ncid, varid, tot_an )
    ierr = nf90_inq_varid( ncid, 'tot_ss', varid )
    ierr = nf90_get_var( ncid, varid, tot_ss )

    !itb...close the file
    ierr = nf90_close( ncid )

    print *,'\t\t load data into the structure'

    !itb...need to load these data into sibtype arrays
    do i = 1,subcount
        sib(i)%prog%ta = ta(subset(i))
        sib(i)%prog%tc = tc(subset(i))
        sib(i)%prog%nsl = nsl(subset(i))
        sib(i)%prog%pco2ap = pco2ap(subset(i))
        sib(i)%prog%d13cca = d13cca(subset(i))
        sib(i)%prog%snow_veg = snow_veg(subset(i))
        sib(i)%prog%snow_age = snow_age(subset(i))
        sib(i)%prog%snow_depth = snow_depth(subset(i))
        sib(i)%prog%snow_mass = snow_mass(subset(i))
        sib(i)%prog%tke = max( tkemin, tke(subset(i)) )
        sib(i)%prog%sha = sha(subset(i))
        sib(i)%stat%coszbar = coszbar(subset(i))
        sib(i)%stat%dayflag = dayflag(subset(i))
        
        sib(i)%prog%capac(1) = capac1(subset(i))
        sib(i)%prog%capac(2) = capac2(subset(i))
        
     !EL..crop vars added
        sib(i)%diag%tempf = tempf(subset(i))
        sib(i)%diag%gdd = gdd(subset(i))
        sib(i)%diag%w_main = w_main(subset(i))
        sib(i)%diag%pd = pd(subset(i))
        sib(i)%diag%pd7 = pd7(subset(i))
        sib(i)%diag%pd7_est = pd7_est(subset(i))
        sib(i)%diag%pdindx7 = pdindx7(subset(i))
        sib(i)%diag%ndf_opt = ndf_opt(subset(i))

print*,'PLANTING DATE 0:',sib%diag%pd
   
        do j = 1, 4
            sib(i)%diag%cum_wt_prev(j) = cum_wt_prev(subset(i),j)
        enddo
       
        do j = 1, 4
            sib(i)%diag%cum_drywt_prev(j) = cum_drywt_prev(subset(i),j)
        enddo


      !EL...end crop vars..

        do k = 1, 12
            sib(i)%diag%tot_an(k) = tot_an(k,subset(i))

            do j = 1, nsoil
                sib(i)%diag%tot_ss(k,j) = tot_ss(k,subset(i),j)
            enddo
        enddo

        do j = 1, 6
            sib(i)%prog%rst(j) = rst(subset(i),j)
        enddo

        do j = 1,nsnow
            k = j - 5
            sib(i)%prog%dz(k)      = dz_snow(subset(i),j)
            sib(i)%prog%node_z(k)  = nz_snow(subset(i),j)
            sib(i)%prog%layer_z(k-1) = lz_snow(subset(i),j)
        enddo

        do j=-nsnow+1,nsoil
            sib(i)%prog%td(j)      = deept(subset(i),j)
            sib(i)%prog%www_liq(j) = www_liq(subset(i),j)
            sib(i)%prog%www_ice(j) = www_ice(subset(i),j)
        enddo

    enddo   !subcount loop

    print *, '\t\t read in sib initial conditions'

!itb...need to manipulate tot_an and tot_ss for restart/initial conditions...
    if(time%sec_year /= 0) then
       jday = nsecond/86400
       month_loop: do j = 1, 11
         if(jday == map_totals(j)) then
           do i=1,subcount
             sib(i)%diag%tot_ss(j+1:13,:) = 0.0_dbl_kind
             sib(i)%diag%tot_an(j+1:13)   = 0.0_dbl_kind
           enddo
           exit month_loop
         endif
       enddo month_loop

    else

     do i=1,subcount
       sib(i)%diag%tot_ss(:,:) = 0.0_dbl_kind
       sib(i)%diag%tot_an(:)   = 0.0_dbl_kind
     enddo

    endif

!EL..manipulating crop vars for restart/initial conditions...
   




      do i=1,subcount
!itb...sanity check
       if(sib(i)%diag%pd > 365) then
          sib(i)%diag%pd = 0.0_int_kind
          sib(i)%diag%pd7 = 0.0_int_kind
          sib(i)%diag%pd7_est = 0.0_int_kind
          sib(i)%diag%pdindx7 = 0.0_int_kind
          sib(i)%diag%ndf_opt = 0.0_int_kind
          sib(i)%diag%tempf = 0.0_dbl_kind
          sib(i)%diag%gdd = 0.0_dbl_kind
          sib(i)%diag%w_main = 0.0001_dbl_kind
          sib(i)%diag%cum_wt_prev(:) = 0.0001_dbl_kind
          sib(i)%diag%cum_drywt_prev(:)   = 0.0001_dbl_kind
        endif

      enddo
   
        

end subroutine read_ic

!===============================================================================
subroutine read_respfactor(sib)
!===============================================================================
! reads in a respfactor from an external file
!
! Modifications:
!  Kevin Schaefer filled in respfactor for any error (status/=0 rather than status>0) (11/11/04)
!
use kinds
use sibtype
use sib_const_module
use sib_io_module
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib

! local variables
integer(kind=int_kind) :: i,j
integer(kind=int_kind) :: status
real(kind=dbl_kind), dimension(nsib,nsoil) :: respfactor


    !     Read the SiB-CO2 respiration factor 
    if(drvr_type=='single')then
        open( unit=3, file=co2_path, form='formatted', status='old', iostat=status) !jk
        do i = 1,nsoil
            read( 3,*, iostat = status ) respfactor(1,i)
        enddo
    else
        open( unit=3, file=co2_path, form='unformatted', status='old', iostat=status )
        read( unit=3, iostat=status ) i
        read( unit=3, iostat=status ) j
        read( unit=3, iostat=status ) respfactor(:,:)
    endif
    close(unit=3)

    if ( status /= 0 ) then
        print *, ' INIT_SIB: error reading in respFactor'
        print *, '    respFactor set globally to 3.0e-6'
        do i=1,nsib
            do j=1,nsoil
                respfactor(i,j) = 3.0e-6_dbl_kind
            enddo
        enddo
    endif

    !itb...copy respfactor into the structure...
    do i=1,subcount
        do j=1,nsoil
            sib(i)%param%respfactor(j) = respfactor(subset(i),j)
      print*,i,j,sib(i)%param%respfactor(j)
        enddo
    enddo

end subroutine read_respfactor
!
!===============================================================================
subroutine soil_properties(sib)
!===============================================================================
! calculates various soil parameters that do not change with time
!
! Modifications:
!  Kevin Schaefer moved soil layer calculates here from read_ic (10/27/04)
!===============================================================================
!
use kinds
use sibtype
use sib_const_module
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib

! local variables
integer(kind=int_kind) :: i, j      ! (-) indeces
real(kind=real_kind) :: tkm        ! (W/m K) mineral conductivity
real(kind=real_kind) :: bd         ! (kg/m^3) bulk density of dry soil material
real(kind=real_kind) :: kroot(12)  ! (?) root density extinction coeficient
real(kind=real_kind) :: totalroot  ! (?) total root density in soil column
real(kind=real_kind) :: ztop       ! (-) normalized depth of soil layer top
real(kind=real_kind) :: zbot       ! (-) normalized depth of soil layer bottom

real(kind=real_kind) :: pot_fc     ! water potential at field capacity (J/kg)
real(kind=real_kind) :: pot_wp     ! water potential at wilt point (J/kg)


!
! assign values of root density profiles
DATA KROOT/3.9,3.9,2.0,5.5,5.5,2.0,5.5,2.0,2.0,5.5,2.0,5.5/


    do i = 1,subcount
        !Bio------------------------------------------------------------------
        !Bio   miscellaneous soil properties
        !Bio-------------------------------------------------------------------


!itb...fixing field capacity and wilting point, based on %sand/%clay basis
!itb...the stress performance is directly tied to FC and WP values. We are
!itb...playing with the 'operating range' that gives us the best model 
!itb...performance. 

        pot_fc = -15.0   ! field capacity (J/kg)

        pot_wp = -1500.0 ! wilt point (J/kg)

        sib(i)%param%fieldcap = sib(i)%param%poros*             &
                   ((pot_fc/9.8)/sib(i)%param%phsat) ** (-1.0 / sib(i)%param%bee)

        sib(i)%param%vwcmin = sib(i)%param%poros *               &
                 ((pot_wp/9.8)/sib(i)%param%phsat) ** (-1.0 / sib(i)%param%bee)


        
        tkm = ( 8.80*sib(i)%param%sandfrac + 2.92*sib(i)%param%clayfrac ) /    &
            ( sib(i)%param%sandfrac + sib(i)%param%clayfrac )

        bd = (1.0 - sib(i)%param%poros) * 2.7E3

        do j=1,nsoil

            sib(i)%param%tkmg(j)    = tkm**(1.0 - sib(i)%param%poros)

            sib(i)%param%tksatu(j)  = sib(i)%param%tkmg(j) * 0.57**sib(i)%param%poros 

            sib(i)%param%tkdry(j)   = (0.135*bd + 64.7) / (2.7E3 - 0.947*bd)

            sib(i)%param%csolid(j)  = (2.128*sib(i)%param%sandfrac       &
                + 2.385*sib(i)%param%clayfrac)/      &
                (sib(i)%param%sandfrac + sib(i)%param%clayfrac)*1.0E6
        enddo

        !Bio-------------------------------------------------------------------
        !Bio  compute soil layer values
        !Bio-------------------------------------------------------------------

        !itb...CLM uses a 'scalez' (0.025) factor to determine soil layer depths. 
        !itb...for now i'm going to use it as well...
        do j=1,nsoil
            sib(i)%prog%node_z(j) = 0.025*(exp(0.5*(j-0.5))-1.0)
        enddo

        sib(i)%prog%dz(1) = 0.5*(sib(i)%prog%node_z(1)+sib(i)%prog%node_z(2))

        do j=2,nsoil-1
            sib(i)%prog%dz(j) = 0.5*(sib(i)%prog%node_z(j+1)-  &
                sib(i)%prog%node_z(j-1))
        enddo

        sib(i)%prog%dz(nsoil) = sib(i)%prog%node_z(nsoil) -   &
            sib(i)%prog%node_z(nsoil-1)
        sib(i)%prog%layer_z(0) = 0.0

        do j=1,nsoil-1
            sib(i)%prog%layer_z(j) = 0.5*(sib(i)%prog%node_z(j) +   &
                sib(i)%prog%node_z(j+1))
        enddo
        sib(i)%prog%layer_z(nsoil) = sib(i)%prog%node_z(nsoil) +    &
            0.5*sib(i)%prog%dz(nsoil)


!itb...seems like I'm always wanting info about soil layers...
!        do j=1,nsoil
!          print'(i5,3f15.5)',j,sib(i)%prog%layer_z(j),sib(i)%prog%node_z(j), &
!                sib(i)%prog%dz(j)
!        enddo
!        stop

        !Bio-------------------------------------------------------------
        !Bio   compute root fractions for each layer
        !Bio-------------------------------------------------------------
!
! total roots
        totalroot = (1.0 - exp(-kroot(int(sib(i)%param%biome))*     &
          sib(i)%prog%layer_z(nsoil))) / kroot(int(sib(i)%param%biome))
!
! root fraction per soil layer
        ztop = 0.0
        do j=1,nsoil
            zbot = ztop + sib(i)%prog%dz(j)
            sib(i)%param%rootf(j) = (exp(-kroot(int(sib(i)%param%biome))*ztop) &
                - exp(-kroot(int(sib(i)%param%biome))*zbot))/ &
                (kroot(int(sib(i)%param%biome)) * totalroot)

!          print*, j,sib(i)%param%rootf(j)

            ztop = zbot

        enddo


        !itb...quick patch to cover some underflow problems...
        if(sib(i)%param%vcover < sib(i)%param%zlt/10.0)  then
            sib(i)%param%vcover = sib(i)%param%vcover * 10.0
        endif
    enddo  ! subcount loop


end subroutine soil_properties

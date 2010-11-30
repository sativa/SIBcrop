!-------------------------------------------------------------------------------
subroutine rtape_sib ( sib, time, rank )
!-------------------------------------------------------------------------------
! This subroutine creates a netcdf restart file for all prognostic variables.
!
! Modifications:
!  Kevin Schaefer added calls to netcdf error handling routine (10/26/04)
!-------------------------------------------------------------------------------

use kinds
#ifdef PGF
use netcdf
use typeSizes
#endif
use sibtype
use timetype
use sib_io_module
use sib_const_module


!jlc...parameters
type(sib_t), dimension(subcount), intent(in) :: sib
type(time_struct), intent(in) :: time
integer(kind=int_kind), intent(in) :: rank

!Bio...local variables
integer(kind=int_kind) :: i,j,ierr,ncid,i1,i2
integer(kind=int_kind) :: vdims(3)
integer(kind=int_kind) :: start(2)
integer(kind=int_kind) :: vcount(2)
integer(kind=int_kind) :: nsibdid   ! dimension id - nsib
integer(kind=int_kind) :: nsoildid  ! dimension id - nsoil
integer(kind=int_kind) :: nsnowdid  ! dimension id - nsnow
integer(kind=int_kind) :: nphysdid  ! dimension id - physiology types(6)
!EL...crop pool no.ID  added..
integer(kind=int_kind) :: npooldid  ! dimension id - crop pool types(4)
!EL....crop pool no.ID..end
integer(kind=int_kind) :: ntotdid   ! dimension id - soil + snow 
integer(kind=int_kind) :: monthid   ! dimension id - number of months
integer(kind=int_kind) :: nsibvid   ! variable id - nsib
integer(kind=int_kind) :: nsoilvid  ! variable id - nsoil
integer(kind=int_kind) :: nsnowvid  ! variable id - nsnow
integer(kind=int_kind) :: vervid    ! variable id - version number
integer(kind=int_kind) :: subcountid! variable id - subcount
integer(kind=int_kind) :: tavid     ! variable id - ta (CAS temp)
integer(kind=int_kind) :: tcvid     ! variable id - tcanopy
integer(kind=int_kind) :: nslvid    ! variable id - # of snow layers
integer(kind=int_kind) :: pco2avid  ! variable id - pco2a
integer(kind=int_kind) :: d13ccaid  ! variable id - d13cca
integer(kind=int_kind) :: svegvid   ! variable id - snow_veg
integer(kind=int_kind) :: sagevid   ! variable id - snow_age
integer(kind=int_kind) :: sdepthid  ! variable id - snow_depth
integer(kind=int_kind) :: smassid   ! variable id - snow_mass
integer(kind=int_kind) :: capac1vid ! variable id - capac1  
integer(kind=int_kind) :: capac2vid ! variable id - capac2  
integer(kind=int_kind) :: coszbarid ! variable id - coszbar
integer(kind=int_kind) :: dayflagid ! variable id - dayflag
integer(kind=int_kind) :: totanid   ! variable id - totan
integer(kind=int_kind) :: tkevid    ! variable id - tke     
integer(kind=int_kind) :: tdvid     ! variable id - td (soil/snow temp)
integer(kind=int_kind) :: wwwlvid   ! variable id - www_liq 
integer(kind=int_kind) :: wwwivid   ! variable id - www_ice 
integer(kind=int_kind) :: dzsvid    ! variable id - dz (snow) 
integer(kind=int_kind) :: lzsvid    ! variable id - layer_z (snow)
integer(kind=int_kind) :: nzsvid    ! variable id - node_z (snow)
integer(kind=int_kind) :: rstvid    ! variable id - rst  
integer(kind=int_kind) :: nsecvid   ! variable id - nsecond 
integer(kind=int_kind) :: nsectemp  ! temporary var to hold time%sec_year
integer(kind=int_kind) :: shavid    ! variable id - sha
integer(kind=int_kind) :: totssid   ! variable id - totss
!EL..crop vars added
integer(kind=int_kind) :: tempfid    ! variable id - tempf 
integer(kind=int_kind) :: gddvid    ! variable id - gdd 
integer(kind=int_kind) :: w_mainid  ! variable id - w_main
integer(kind=int_kind) :: w_main_potid  ! variable id - w_main_pot
integer(kind=int_kind) :: pdvid     ! variable id - pd 
integer(kind=int_kind) :: pd7vid    ! variable id - pd7
integer(kind=int_kind) :: pd7_estid ! variable id - pd7_est 
integer(kind=int_kind) :: emerg_did ! variable id - emerg_d
integer(kind=int_kind) :: pdindx7id ! variable id - pdindx7
integer(kind=int_kind) :: ndf_optid ! variable id - ndf_opt
integer(kind=int_kind) :: nd_emergid ! variable id - nd_emerg
integer(kind=int_kind) :: cum_wt_id    ! variable id - cum_wt
integer(kind=int_kind) :: cum_drywt_id ! variable id - cum_drywt
!EL..crop vars end..


character*10 curmon  !jlc
character*256 rmon  !jk

integer(kind=int_kind), dimension(nsib) :: nsl


real(kind=dbl_kind), dimension(nsib) ::  &
    ta,         &
    tc,         &
    pco2ap,     &
    d13cca,     &
    snow_veg,   &
    snow_age,   &
    snow_depth, &
    snow_mass,  &
    tke,        &
    sha,        &
    capac1,     &
    capac2,     &
    coszbar,    &
    dayflag

!EL..crop variables integer(kind=int_kind)

integer(kind=int_kind), dimension(nsib) :: pd
integer(kind=int_kind), dimension(nsib) :: pd7
integer(kind=int_kind), dimension(nsib) :: pd7_est
integer(kind=int_kind), dimension(nsib) :: emerg_d
integer(kind=int_kind), dimension(nsib) :: pdindx7
integer(kind=int_kind), dimension(nsib) :: ndf_opt
integer(kind=int_kind), dimension(nsib) :: nd_emerg

!EL...crop variables (real(kind=dbl_kind))

real(kind=dbl_kind), dimension(nsib) :: tempf,gdd,w_main,w_main_pot

real(kind=dbl_kind), dimension(nsib,4) :: cum_wt,cum_drywt

!EL...end crop vars (real(kind=dbl_kind))     
    
real(kind=dbl_kind), dimension(12,nsib) :: tot_an

real(kind=dbl_kind), dimension(nsib,6) :: rst

real(kind=dbl_kind), dimension(nsib,-nsnow+1:nsoil) ::  &
    deept,   & ! ->sib%prog%deept
    www_liq, & ! ->sib%prog%www_liq
    www_ice    ! ->sib%prog%www_ice

real(kind=dbl_kind), dimension(nsib,-nsnow+1:0) ::  &
    dz_snow,    &
    nz_snow

real(kind=dbl_kind), dimension(nsib,-nsnow:0) :: lz_snow

real(kind=dbl_kind), dimension(12, nsib, nsoil) :: tot_ss


    !jlc...copy data into temporary arrays initialized to 1.e36
    ta(:) = 1.e36
    tc(:) = 1.e36
    nsl(:) = 0
    pco2ap(:) = 1.e36
    d13cca(:) = 1.e36
    snow_veg(:) = 1.e36
    snow_age(:) = 1.e36
    snow_depth(:) = 1.e36
    snow_mass(:) = 1.e36
    tke(:) = 1.e36
    sha(:) = 1.e36
    capac1(:) = 1.e36
    capac2(:) = 1.e36
    coszbar(:) = 1.e36
    tot_an(:,:) = 1.e36
    rst(:,:) = 1.e36
    deept(:,:) = 1.e36
    www_liq(:,:) = 1.e36
    www_ice(:,:) = 1.e36
    dz_snow(:,:) = 1.e36
    nz_snow(:,:) = 1.e36
    lz_snow(:,:) = 1.e36
    tot_ss(:,:,:) = 1.e36
!EL...crop vars added..
    tempf(:)   = 1.e36
    gdd(:)     = 1.e36
    w_main(:)  = 1.e36
    w_main_pot(:)  = 1.e36
    pd(:)      = 0
    pd7(:)     = 0
    pd7_est(:) = 0
    emerg_d(:) = 0
    pdindx7(:) = 0
    ndf_opt(:) = 0
    nd_emerg(:) = 0
    cum_wt(:,:) = 1.e36
    cum_drywt(:,:) = 1.e36
!EL..end initializing crop vars.
    
    do i = 1, subcount
        ta(subset(i)) = sib(i)%prog%ta
        tc(subset(i)) = sib(i)%prog%tc
        nsl(subset(i)) = sib(i)%prog%nsl
        pco2ap(subset(i)) = sib(i)%prog%pco2ap
        d13cca(subset(i)) = sib(i)%prog%d13cca
        snow_veg(subset(i)) = sib(i)%prog%snow_veg
        snow_age(subset(i)) = sib(i)%prog%snow_age
        snow_depth(subset(i)) = sib(i)%prog%snow_depth
        snow_mass(subset(i)) = sib(i)%prog%snow_mass
        tke(subset(i)) = sib(i)%prog%tke
        sha(subset(i)) = sib(i)%prog%sha
        capac1(subset(i)) = sib(i)%prog%capac(1)
        capac2(subset(i)) = sib(i)%prog%capac(2)
        coszbar(subset(i)) = sib(i)%stat%coszbar
        dayflag(subset(i)) = sib(i)%stat%dayflag
        tot_an(:,subset(i)) = sib(i)%diag%tot_an(1:12)
!EL..crop vars added.
        tempf(subset(i)) = sib(i)%diag%tempf
        gdd(subset(i)) = sib(i)%diag%gdd
        w_main(subset(i)) = sib(i)%diag%w_main
	w_main_pot(subset(i)) = sib(i)%diag%w_main_pot
        pd(subset(i)) = sib(i)%diag%pd
        pd7(subset(i)) = sib(i)%diag%pd7
        pd7_est(subset(i)) = sib(i)%diag%pd7_est
        emerg_d(subset(i)) = sib(i)%diag%emerg_d
        pdindx7(subset(i)) = sib(i)%diag%pdindx7
        ndf_opt(subset(i)) = sib(i)%diag%ndf_opt
        nd_emerg(subset(i)) = sib(i)%diag%nd_emerg
       
!EL..crop vars contd below..
    enddo

    do j = 1, nsoil
        do i = 1, subcount
            tot_ss(:,subset(i),j) = sib(i)%diag%tot_ss(1:12,j)
        enddo
    enddo

    do j=1,6
        do i=1,subcount
            rst(subset(i),j) = sib(i)%prog%rst(j)
        enddo
    enddo

    do j=-nsnow+1,nsoil
        do i=1,subcount
            deept(subset(i),j)   = sib(i)%prog%td(j)
            www_liq(subset(i),j) = sib(i)%prog%www_liq(j)
            www_ice(subset(i),j) = sib(i)%prog%www_ice(j)
            if (j <= 0) then
                dz_snow(subset(i),j) = sib(i)%prog%dz(j)
                nz_snow(subset(i),j) = sib(i)%prog%node_z(j)
            endif
        enddo
    enddo
!EL...crop vars contd..
    do j=1,4
        do i=1,subcount
           cum_wt(subset(i),j) = sib(i)%diag%cum_wt(j)
        enddo
    enddo

    do j=1,4
        do i=1,subcount
            cum_drywt(subset(i),j) = sib(i)%diag%cum_drywt(j)
        enddo
    enddo
!EL..crop vars end..
    
    !jlc...is this a typo in sibtype? does layer_z have the same
    !      indicies as node_z and dz?
    !itb...layer_z has one more value, since it holds the edges
    !      of the layers

    do j=-nsnow,0
        do i=1,subcount
            lz_snow(subset(i),j) = sib(i)%prog%layer_z(j)
        enddo
    enddo


    !-----------------------------------------------------------------
    ! This subroutine writes out a restart file for SiBDRV
    !-----------------------------------------------------------------


    !itb...open netcdf history (restart) file
    
    write(curmon, '(i4.4,i2.2,a,i3.3)') time%year, time%month, 'p', rank
    rmon = trim(out_path)//"sib_r"//trim(curmon)//".nc" !jk
    print*, 'write restart ', trim(rmon)
    ierr = nf90_create( rmon, nf90_clobber, ncid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',1)


    !itb...dimensions...
    ierr = nf90_def_dim( ncid, 'nsib', nsib, nsibdid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',2)
    ierr = nf90_def_dim( ncid, 'nsoil', nsoil, nsoildid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',3)
    ierr = nf90_def_dim( ncid, 'nsnow', nsnow, nsnowdid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',4)
    ierr = nf90_def_dim( ncid, 'nphys', physmax+1, nphysdid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',5)
    ierr = nf90_def_dim( ncid, 'ntot', nsoil+nsnow, ntotdid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',6)
    ierr = nf90_def_dim( ncid, 'nmonths', 12, monthid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',7)
    ierr = nf90_def_dim( ncid, 'npool', 4, npooldid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',8)


    !itb... define scalar variables...
    vdims(1) = monthid
    vdims(2) = nsibdid
    vdims(3) = ntotdid

    ierr = nf90_def_var( ncid, 'nsib', nf90_int, nsibvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',9)
    ierr = nf90_def_var( ncid, 'nsoil', nf90_int, nsoilvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',10)
    ierr = nf90_def_var( ncid, 'nsnow', nf90_int, nsnowvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',11)
    ierr = nf90_def_var( ncid, 'version', nf90_float, vervid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',12) 

!itb...if nsecond = 86400 * 356, set nsecond = 0
    nsectemp = time%sec_year
    if(nsectemp == 31536000) nsectemp = 0


    ierr = nf90_def_var( ncid, 'nsecond', nf90_int, nsecvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',13)
    ierr = nf90_def_var( ncid, 'subcount', nf90_int, subcountid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',14)


    !itb...define vector (length=nsib ) variables...
    ierr = nf90_def_var( ncid, 'ta', nf90_double, vdims(2), tavid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',15)
    ierr = nf90_def_var( ncid, 'tc', nf90_double, vdims(2), tcvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',16)
    ierr = nf90_def_var( ncid, 'nsl', nf90_int, vdims(2), nslvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',17)
    ierr = nf90_def_var( ncid, 'pco2a', nf90_double, vdims(2), pco2avid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',18)
    ierr = nf90_def_var( ncid, 'd13cca', nf90_double, vdims(2), d13ccaid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',19)
    ierr = nf90_def_var( ncid, 'snow_veg', nf90_double, vdims(2), svegvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',20)
    ierr = nf90_def_var( ncid, 'snow_age', nf90_double, vdims(2), sagevid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',21)
    ierr = nf90_def_var( ncid, 'snow_depth', nf90_double, vdims(2), sdepthid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',22)
    ierr = nf90_def_var( ncid, 'snow_mass', nf90_double, vdims(2), smassid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',23)
    ierr = nf90_def_var( ncid, 'capac1', nf90_double, vdims(2), capac1vid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',24)
    ierr = nf90_def_var( ncid, 'capac2', nf90_double, vdims(2), capac2vid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',25)
    ierr = nf90_def_var( ncid, 'coszbar', nf90_double, vdims(2), coszbarid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',26)
    ierr = nf90_def_var( ncid, 'dayflag', nf90_double, vdims(2), dayflagid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',27)
    ierr = nf90_def_var( ncid, 'tke', nf90_double, vdims(2), tkevid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',28)
    ierr = nf90_def_var( ncid, 'sha', nf90_double, vdims(2), shavid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',29)
    ierr = nf90_def_var( ncid, 'tot_an', nf90_double, vdims(1:2), totanid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',30)


    !itb...define 2-D variables...
    ierr = nf90_def_var( ncid, 'td', nf90_double, vdims(2:3), tdvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',31)
    ierr = nf90_def_var( ncid, 'www_liq', nf90_double, vdims(2:3), wwwlvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',32)
    ierr = nf90_def_var( ncid, 'www_ice', nf90_double, vdims(2:3), wwwivid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',33)

    !EL...crop vars added..
    !EL..define crop vars..
    !EL..define..define vector crop variables...
    
    ierr = nf90_def_var( ncid, 'tempf', nf90_double, vdims(2), tempfid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',34)
    ierr = nf90_def_var( ncid, 'gdd', nf90_double, vdims(2), gddvid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',35)
    ierr = nf90_def_var( ncid, 'w_main', nf90_double, vdims(2), w_mainid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',36)
    ierr = nf90_def_var( ncid, 'w_main_pot', nf90_double, vdims(2), w_main_potid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',37)
    ierr = nf90_def_var( ncid, 'pd', nf90_int, vdims(2), pdvid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',37)
    ierr = nf90_def_var( ncid, 'pd7', nf90_int, vdims(2), pd7vid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',38)
    ierr = nf90_def_var( ncid, 'pd7_est', nf90_int, vdims(2), pd7_estid)
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',38)
    ierr = nf90_def_var( ncid, 'emerg_d', nf90_int, vdims(2), emerg_did)
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',39)
    ierr = nf90_def_var( ncid, 'pdindx7', nf90_int, vdims(2), pdindx7id) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',40)
    ierr = nf90_def_var( ncid, 'ndf_opt', nf90_int, vdims(2), ndf_optid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',41)
    ierr = nf90_def_var( ncid, 'nd_emerg', nf90_int, vdims(2), nd_emergid) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',41)
   
    vdims(3) = npooldid
    ierr = nf90_def_var( ncid, 'cum_wt', nf90_double, vdims(2:3),cum_wt_id)  
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',42)
    ierr = nf90_def_var( ncid, 'cum_drywt', nf90_double, vdims(2:3), cum_drywt_id) 
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',43)
     
    !EL...end defining crop vars..

    vdims(3) = nsnowdid
    ierr = nf90_def_var( ncid, 'dzsnow', nf90_double, vdims(2:3), dzsvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',44)
    ierr = nf90_def_var( ncid, 'nzsnow', nf90_double, vdims(2:3), nzsvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',45)
    ierr = nf90_def_var( ncid, 'lzsnow', nf90_double, vdims(2:3), lzsvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',46)

    vdims(3) = nphysdid
    ierr = nf90_def_var( ncid, 'rst', nf90_double, vdims(2:3), rstvid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',47)

    vdims(3) = nsoildid
    ierr = nf90_def_var( ncid, 'tot_ss', nf90_double, vdims, totssid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',48)

    !itb...take file out of define mode, into data mode
    ierr = nf90_enddef( ncid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',49)

    !itb...load the variables...
    ierr = nf90_put_var( ncid, nsibvid, nsib )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',50)
    ierr = nf90_put_var( ncid, nsoilvid, nsoil )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',51)
    ierr = nf90_put_var( ncid, nsnowvid, nsnow )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',52)
    ierr = nf90_put_var( ncid, nsecvid, nsectemp )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',53)
    ierr = nf90_put_var( ncid, vervid, version )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',54)
    ierr = nf90_put_var( ncid, subcountid, subcount )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',55)

    ierr = nf90_put_var( ncid, tavid, ta )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',56)
    ierr = nf90_put_var( ncid, tcvid, tc )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',57)
    ierr = nf90_put_var( ncid, nslvid, nsl )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',58)
    ierr = nf90_put_var( ncid, pco2avid, pco2ap )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',59)
    ierr = nf90_put_var( ncid, d13ccaid, d13cca )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',60)
    ierr = nf90_put_var( ncid, svegvid, snow_veg )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',61)
    ierr = nf90_put_var( ncid, sagevid, snow_age )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',62)
    ierr = nf90_put_var( ncid, sdepthid, snow_depth )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',63)
    ierr = nf90_put_var( ncid, smassid, snow_mass )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',64)
    ierr = nf90_put_var( ncid, tkevid, tke )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',65)
    ierr = nf90_put_var( ncid, shavid, sha )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',66)
    ierr = nf90_put_var( ncid, coszbarid, coszbar )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',67)
    ierr = nf90_put_var( ncid, dayflagid, dayflag )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',68)
    ierr = nf90_put_var( ncid, totanid, tot_an )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',69)

    !jlc...these are the temporary arrays
    ierr = nf90_put_var( ncid, tdvid, deept )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',70)
    ierr = nf90_put_var( ncid, wwwlvid, www_liq )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',71)
    ierr = nf90_put_var( ncid, wwwivid, www_ice )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',72)
    ierr = nf90_put_var( ncid, rstvid, rst )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',73)

    !itb...slabs...
    ierr = nf90_put_var( ncid, capac1vid, capac1 )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',74)
    ierr = nf90_put_var( ncid, capac2vid, capac2 )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',75)

    
    !EL..load the crop vars.
   
    ierr = nf90_put_var( ncid, tempfid, tempf )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',76)
    ierr = nf90_put_var( ncid, gddvid, gdd )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',77)
    ierr = nf90_put_var( ncid, w_mainid, w_main )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',78)
    ierr = nf90_put_var( ncid, w_main_potid, w_main_pot )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',79)
    ierr = nf90_put_var( ncid, pdvid, pd )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',80)
    ierr = nf90_put_var( ncid, pd7vid, pd7 )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',81)
    ierr = nf90_put_var( ncid, pd7_estid, pd7_est )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',82)
    ierr = nf90_put_var( ncid, emerg_did, emerg_d )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',83)
    ierr = nf90_put_var( ncid, pdindx7id, pdindx7 )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',84)
    ierr = nf90_put_var( ncid, ndf_optid, ndf_opt )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',85)
    ierr = nf90_put_var( ncid, nd_emergid, nd_emerg )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',86)
    ierr = nf90_put_var( ncid, cum_wt_id, cum_wt )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',87)
    ierr = nf90_put_var( ncid, cum_drywt_id, cum_drywt )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',88)
    !EL...crop vars end..


    !itb...netcdf does not deal well with negative indices: the 
    !itb...dz/node_z/layer_z arrays appear to go from (1:nsoil+nsnow)
    !itb...rather than (-nsnow+1:nsoil)
    start(1) = 1
    start(2) = 1
    vcount(1) = nsib
    vcount(2) = nsnow


    ierr = nf90_put_var( ncid, dzsvid, dz_snow, start, vcount )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',89)
    ierr = nf90_put_var( ncid, nzsvid, nz_snow, start, vcount )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',90)
    ierr = nf90_put_var( ncid, lzsvid, lz_snow, start,vcount )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',91)


    ierr = nf90_put_var( ncid, totssid, tot_ss )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',92)

    !itb...close the file
    ierr = nf90_close( ncid )
    if(ierr/=nf90_noerr) call handle_err(ierr,'rtape',93)

end subroutine rtape_sib

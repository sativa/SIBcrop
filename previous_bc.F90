subroutine previous_bc( sib, time)

use kinds
use sibtype
use timetype
use sib_io_module
use sib_const_module
use sib_bc_module
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time

! local variables
integer(kind=int_kind) :: ntest1    ! compare nsib value of file to simulation
integer(kind=int_kind) :: i, k
real :: dummy
real(kind=real_kind), dimension(subcount) :: curndvi  ! NDVI variable
character*100 filename  ! filename used to read in ndvi data

    if( drvr_type == 'single' ) then
         write (filename, "(a,a1,i4)") trim(param_path), '_', time%pyear
         print*,'      reading file: ',trim(filename)

        open(unit=32, file=trim(filename), form='formatted')   
        read(32,*)ntest1
        if(ntest1.ne.nsib)  stop &
            ' file sib_bc no match with model for nsib'

        ! read in time-invariant boundary condition variables
        read(32,*)
        read(32,*)
        read(32,*) sib%param%biome
        read(32,*) dummy
        sib%param%phystype(1) = dummy
        read(32,*) dummy
        sib%param%phystype(2) = dummy
        read(32,*) dummy
        sib%param%phystype(3) = dummy
        read(32,*) dummy
        sib%param%phystype(4) = dummy
        read(32,*) dummy
        sib%param%phystype(5) = dummy
        read(32,*) sib%param%z2
        read(32,*) sib%param%z1
        read(32,*) sib%param%vcover
        read(32,*) sib%param%chil
        read(32,*) sib%param%phc
        read(32,*) sib%param%tran(1,1)
        read(32,*) sib%param%tran(1,2)
        read(32,*) sib%param%tran(2,1)
        read(32,*) sib%param%tran(2,2)
        read(32,*) sib%param%ref(1,1)
        read(32,*) sib%param%ref(1,2)
        read(32,*) sib%param%ref(2,1)
        read(32,*) sib%param%ref(2,2)
        read(32,*) sib%param%vmax0(1)
        read(32,*) sib%param%vmax0(2)
        read(32,*) sib%param%vmax0(3)
        read(32,*) sib%param%vmax0(4)
        read(32,*) sib%param%vmax0(5)
        read(32,*) sib%param%effcon(1)
        read(32,*) sib%param%effcon(2)
        read(32,*) sib%param%effcon(3)
        read(32,*) sib%param%effcon(4)
        read(32,*) sib%param%effcon(5)
        read(32,*) sib%param%gradm(1)
        read(32,*) sib%param%gradm(2)
        read(32,*) sib%param%gradm(3)
        read(32,*) sib%param%gradm(4)
        read(32,*) sib%param%gradm(5)
        read(32,*) sib%param%binter(1)
        read(32,*) sib%param%binter(2)
        read(32,*) sib%param%binter(3)
        read(32,*) sib%param%binter(4)
        read(32,*) sib%param%binter(5)
        read(32,*) sib%param%atheta(1)
        read(32,*) sib%param%atheta(2)
        read(32,*) sib%param%atheta(3)
        read(32,*) sib%param%atheta(4)
        read(32,*) sib%param%atheta(5)
        read(32,*) sib%param%btheta(1)
        read(32,*) sib%param%btheta(2)
        read(32,*) sib%param%btheta(3)
        read(32,*) sib%param%btheta(4)
        read(32,*) sib%param%btheta(5)
        read(32,*) sib%param%trda(1)
        read(32,*) sib%param%trda(2)
        read(32,*) sib%param%trda(3)
        read(32,*) sib%param%trda(4)
        read(32,*) sib%param%trda(5)
        read(32,*) sib%param%trdm(1)
        read(32,*) sib%param%trdm(2)
        read(32,*) sib%param%trdm(3)
        read(32,*) sib%param%trdm(4)
        read(32,*) sib%param%trdm(5)
        read(32,*) sib%param%trop(1)
        read(32,*) sib%param%trop(2)
        read(32,*) sib%param%trop(3)
        read(32,*) sib%param%trop(4)
        read(32,*) sib%param%trop(5)
        read(32,*) sib%param%respcp(1)
        read(32,*) sib%param%respcp(2)
        read(32,*) sib%param%respcp(3)
        read(32,*) sib%param%respcp(4)
        read(32,*) sib%param%respcp(5)
        read(32,*) sib%param%slti(1)
        read(32,*) sib%param%slti(2)
        read(32,*) sib%param%slti(3)
        read(32,*) sib%param%slti(4)
        read(32,*) sib%param%slti(5)
        read(32,*) sib%param%shti(1)
        read(32,*) sib%param%shti(2)
        read(32,*) sib%param%shti(3)
        read(32,*) sib%param%shti(4)
        read(32,*) sib%param%shti(5)
        read(32,*) sib%param%hltii(1)
        read(32,*) sib%param%hltii(2)
        read(32,*) sib%param%hltii(3)
        read(32,*) sib%param%hltii(4)
        read(32,*) sib%param%hltii(5)
        read(32,*) sib%param%hhti(1)
        read(32,*) sib%param%hhti(2)
        read(32,*) sib%param%hhti(3)
        read(32,*) sib%param%hhti(4)
        read(32,*) sib%param%hhti(5)
        read(32,*) sib%param%soref(1)
        read(32,*) sib%param%soref(2)
        read(32,*) sib%param%bee
        read(32,*) sib%param%phsat
        read(32,*) sib%param%satco
        read(32,*) sib%param%poros
        read(32,*) sib%param%slope
        read(32,*) sib%param%wopt
        read(32,*) sib%param%wsat
        read(32,*) sib%param%zm
        read(32,*) sib%param%sandfrac
        read(32,*) sib%param%clayfrac

        ! read in previous month's time-variant boundary condition variables
        do k = 1,time%pmonth
            read(32,*)
            read(32,*) sib%param%aparc2
            read(32,*) sib%param%zlt2
            read(32,*) sib%param%green2
            read(32,*) sib%param%z0d2
            read(32,*) sib%param%zp_disp2
            read(32,*) sib%param%rbc2
            read(32,*) sib%param%rdc2
            read(32,*) sib%param%gmudmu2
            read(32,*) sib%param%d13cresp2
            read(32,*) sib%param%physfrac2(1)
            read(32,*) sib%param%physfrac2(2)
            read(32,*) sib%param%physfrac2(3)
            read(32,*) sib%param%physfrac2(4)
            read(32,*) sib%param%physfrac2(5)
        enddo

        close( 32 )

    else    ! global run
    
        ! make sure prevndvi is allocated
        allocate(prevndvi(subcount))

        ! 2 months previous
        !write (filename, "(a,a1,i4,a3)") trim(param_path), '_', &
        !    time%ppyear, '.nc'

        !print*,'      reading file: ',trim(filename)
        !!kdcorbin, 02/11 - modified variables passed to read_ndvi
        !call read_ndvi(filename, sib, time )

        ! 1 month previous
        write (filename, "(a,a1,i4,a3)") trim(param_path), '_', &
            time%pyear, '.nc'
        print*,'      reading file: ',trim(filename)
        call read_ndvi(filename, sib, time )
            
        ! calculate time-variant boundary conditions
        print*,'   calculating time-variant boundary conditions'
        call calculate_td (sib, time%mid_month(time%pmonth), curndvi)

        ! copy current fpar/lai to previous values and fill in initial d13cresp and
        !   physfrac values
        do i=1,subcount
            prevndvi(i) = curndvi(i)
            sib(i)%param%d13cresp2 = sib(i)%param%d13cresp
            do k = 1, physmax
                !kdcorbin, 02/11 - added physfrac1
                sib(i)%param%physfrac1(k) = sib(i)%param%physfrac(k)
                sib(i)%param%physfrac2(k) = sib(i)%param%physfrac(k)
            enddo
        enddo
    endif

end subroutine previous_bc


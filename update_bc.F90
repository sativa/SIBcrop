subroutine update_bc( sib, time)

!Program to update time-dependent boudary conditions
!kdcorbin, 02/11

use sibtype
use timetype
use sib_io_module, only: param_path, drvr_type
use sib_const_module, only : subcount

implicit none

! parameters
type(sib_t), dimension(subcount),intent(inout) :: sib
type(time_struct), intent(in) :: time

!local variables
integer(kind=int_kind) :: i,k
real :: dummy  !used to read in variables that are real in file,
                       !but need to be cast as ints
character*100 filename  !filename used to read in ndvi data

    if( drvr_type == 'single' ) then
        write (filename, "(a,a1,i4)") trim(param_path), '_', time%nyear
        open(unit=32, file=trim(filename), form='formatted')   
 
         do i=1,104
               read(32,*) dummy

              if(i == 4) then
                   !print*,'biome =',dummy
                  sib(1)%param%biome = dummy
              endif
         enddo
       
        ! set the previous variables - kdcorbin, 03/11
         if (sib(1)%param%biome < 20) then
           sib(1)%param%aparc1 = sib(1)%param%aparc2
           sib(1)%param%zlt1 = sib(1)%param%zlt2
           sib(1)%param%green1 = sib(1)%param%green2
           sib(1)%param%z0d1 = sib(1)%param%z0d2
           sib(1)%param%zp_disp1 = sib(1)%param%zp_disp2
           sib(1)%param%rbc1 = sib(1)%param%rbc2
           sib(1)%param%rdc1 = sib(1)%param%rdc2
           sib(1)%param%gmudmu1 = sib(1)%param%gmudmu2
          endif

          sib(1)%param%d13cresp1 = sib(1)%param%d13cresp2
          sib(1)%param%physfrac1 = sib(1)%param%physfrac2

          ! read in the phys fracs from the current month
          do k = 1,time%nmonth
              read (32,*)
        
              if (sib(1)%param%biome < 20) then
                 read(32,*) sib(1)%param%aparc2
                 read(32,*) sib(1)%param%zlt2
                 read(32,*) sib(1)%param%green2
                 read(32,*) sib(1)%param%z0d2
                 read(32,*) sib(1)%param%zp_disp2
                 read(32,*) sib(1)%param%rbc2
                 read(32,*) sib(1)%param%rdc2
                 read(32,*) sib(1)%param%gmudmu2
              else
                 do i=1,8
                      read(32,*) dummy
                 enddo
              endif

              read(32,*) sib(1)%param%d13cresp2
              read(32,*) sib(1)%param%physfrac2(1)
              read(32,*) sib(1)%param%physfrac2(2)
              read(32,*) sib(1)%param%physfrac2(3)
              read(32,*) sib(1)%param%physfrac2(4)
              read(32,*) sib(1)%param%physfrac2(5)
          enddo

       close( 32 )

    endif
    
end subroutine update_bc

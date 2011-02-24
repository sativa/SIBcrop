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
integer(kind=int_kind) :: ntest1  ! compare nsib value of file to simulation
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
       
        ! set the previous fractions
        sib(1)%param%physfrac1 = sib(1)%param%physfrac2

        ! read in the phys fracs from the current month
        do k = 1,time%nmonth
            read (32,*)
        
            do i=1,9
                  read(32,*) dummy
            enddo
            read(32,*) sib(1)%param%physfrac2(1)
            read(32,*) sib(1)%param%physfrac2(2)
            read(32,*) sib(1)%param%physfrac2(3)
            read(32,*) sib(1)%param%physfrac2(4)
            read(32,*) sib(1)%param%physfrac2(5)
        enddo
        close( 32 )

    endif
    
end subroutine update_bc

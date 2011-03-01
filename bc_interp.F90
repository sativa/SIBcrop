!-----------------------------------------------------------------
subroutine bc_interp( sib, time )
!-----------------------------------------------------------------
! interpolates between bc data points
!
! Modifications:
!  Kevin Schaefer created subrotine from bc_update_sib code (3/14/03)
!  Kevin Schaefer removed print statements (3/17/03)
!  Kevin Schaefer moved isotope interpolation to separate routine (2/21/03)
!  Kevin Schaefer deleted LAI underflow patch, already done in laigrn (4/4/03)
!  Kevin Schaefer added ndvi interpolation (7/9/03)
!  Kevin Schaefer changed ndvi variable names to match bc convention (3/15/05)
!         prevndvi to ndvi1  curndvi to ndvi2 
!-----------------------------------------------------------------

use sibtype
use timetype
use sib_const_module
use sib_bc_module
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time

! local variables
integer i, k  ! indices
real tpgf1    ! scaling factor between 1st bc value and current time
real tpgf2    ! scaling factor between 2nd bc value and current time
integer k1    ! 1st month index for interpolation scaling factors
integer k2    ! 2nd month index for interpolation scaling factors

    ! Calculate scaling factors for interpolation of boundary conditions
    if ( time%doy < time%mid_month(time%month) ) then
        k1 = time%month - 1
        k2 = time%month 
        if ( k1 < 1 ) k1 = 12
        tpgf1 = (time%mid_month(k2) - time%doy) * 2.0 /  &
            real(time%days_per_month(k1) + time%days_per_month(k2))
        tpgf2 = 1.0 - tpgf1
    else
        k1 = time%month
        k2 = time%month + 1
        if ( k2 > 12 ) k2 = 1
        tpgf2 = (time%doy - time%mid_month(k1)) * 2.0 /  &
            real(time%days_per_month(k1) + time%days_per_month(k2))
        tpgf1 = 1.0 - tpgf2
    endif

    ! update boundary conditions
    do i = 1, subcount
        sib(i)%param%aparc    = tpgf1*sib(i)%param%aparc1    + tpgf2*sib(i)%param%aparc2
        sib(i)%param%zlt      = tpgf1*sib(i)%param%zlt1      + tpgf2*sib(i)%param%zlt2
        sib(i)%param%green    = tpgf1*sib(i)%param%green1    + tpgf2*sib(i)%param%green2
        sib(i)%param%z0d      = tpgf1*sib(i)%param%z0d1      + tpgf2*sib(i)%param%z0d2
        sib(i)%param%zp_disp  = tpgf1*sib(i)%param%zp_disp1  + &
                      tpgf2*sib(i)%param%zp_disp2

        sib(i)%param%cc1      = tpgf1*sib(i)%param%rbc1      + tpgf2*sib(i)%param%rbc2
        sib(i)%param%cc2      = tpgf1*sib(i)%param%rdc1      + tpgf2*sib(i)%param%rdc2
        sib(i)%param%gmudmu   = tpgf1*sib(i)%param%gmudmu1   + &
                      tpgf2*sib(i)%param%gmudmu2
        sib(i)%param%d13cresp = tpgf1*sib(i)%param%d13cresp1 + &
                      tpgf2*sib(i)%param%d13cresp2

        do k=1,physmax
           sib(i)%param%physfrac(k) = tpgf1*sib(i)%param%physfrac1(k) +    &
                                   tpgf2*sib(i)%param%physfrac2(k)
        enddo

       !kdcorbin, 03/11 - calculate daily scatp, scatg and park
        sib(i)%param%scatp = sib(i)%param%green * &
                                         (sib(i)%param%tran(1,1) + sib(i)%param%ref(1,1)) + &
                                         (1-sib(i)%param%green) * &
                                         (sib(i)%param%tran(1,2) + sib(i)%param%ref(1,2))
        sib(i)%param%scatg = sib(i)%param%tran(1,1) + sib(i)%param%ref(1,1)
        sib(i)%param%park = sqrt(1.-sib(i)%param%scatp)*sib(i)%param%gmudmu

    enddo

end subroutine bc_interp

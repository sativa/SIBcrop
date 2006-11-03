!----------------------------------------------------------------------
subroutine combine_snow(sib)
!----------------------------------------------------------------------
!
!   Based on CLM subroutine CLM_COMBIN
!
!   CLM web info: http://clm.gsfc.nasa.gov
!
!   Description:
!   This subroutine checks for elements which are below the prescribed 
!   minimum for thickness or mass.  If the snow element thickness or 
!   mass is less than a prescribed minimum, then it is combined with a 
!   neighboring element.  The subroutine subdivide_snow then executes 
!   the combination of mass and energy.
!
!   Revision History:
!   15 September 1999: Yongjiu Dai; initial code
!   15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!   30 January   2002: Ian Baker; SiB integration
!
!----------------------------------------------------------------------

use kinds
use sibtype

use sib_const_module, only: &
    nsnow,  &
    denice, &
    denh2o, &
    dti


implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  


!   local variables
integer(kind=int_kind) :: i,i0,j,k,l,m  ! loop variables
integer(kind=int_kind) :: nsnow_old ! copy of number of snow layers
integer(kind=int_kind) :: mssi
integer(kind=int_kind) :: neibor

real(kind=dbl_kind) :: zwice  ! snow-total ice (kg m^-2)
real(kind=dbl_kind) :: zwliq  ! snow-total liquid (kg m^-2)
real(kind=dbl_kind) :: dzmin(5) 
! minimum depth for snow layers (m)
!----------------------------------------------------------------------



data dzmin/0.010,0.015,0.025,0.055,0.115/


!    if ( sib%prog%nsl == 0 ) return

    nsnow_old = sib%prog%nsl

    do j=nsnow_old+1,0


        !itb...Im still not comfortable with the threshold amount to maintain
        !itb...a snow layer. Its been 0.1 (kg/m^2 water) in the past...

        if(sib%prog%www_ice(j) <= 0.05 ) then

            !itb...need to prevent supersaturating the soil column.
            if (  (sib%prog%www_ice(j)/(sib%prog%dz(j)*denice)) + &
                (sib%prog%www_liq(j)/(sib%prog%dz(j)*denh2o)) &
                .gt. sib%param%poros) then
                sib%diag%roffo = sib%diag%roffo + &
                    (sib%prog%www_ice(j) + sib%prog%www_liq(j)) * dti

            else
                sib%prog%www_liq(j+1) = sib%prog%www_liq(j+1) + &
                    sib%prog%www_liq(j)
                sib%prog%www_ice(j+1) = sib%prog%www_ice(j+1) + &
                    sib%prog%www_ice(j)
            endif

            !...shift all elements above this down one.

            if( j > sib%prog%nsl+1 .and. sib%prog%nsl < -1) then
                do k = j,sib%prog%nsl+2, -1
                    sib%prog%td(j)      = sib%prog%td(j-1)
                    sib%prog%www_liq(j) = sib%prog%www_liq(j-1)
                    sib%prog%www_ice(j) = sib%prog%www_ice(j-1)
                    sib%prog%dz(j)      = sib%prog%dz(j-1)
                enddo
            endif

            sib%prog%nsl = sib%prog%nsl + 1

            !itb...zero out the layer that just disappeared
            sib%prog%td(sib%prog%nsl) = 0.0_dbl_kind
            sib%prog%dz(sib%prog%nsl) = 0.0_dbl_kind
            sib%prog%layer_z(sib%prog%nsl-1) = 0.0_dbl_kind
            sib%prog%node_z(sib%prog%nsl) = 0.0_dbl_kind
            sib%prog%www_liq(sib%prog%nsl) = 0.0_dbl_kind
            sib%prog%www_ice(sib%prog%nsl) = 0.0_dbl_kind

        endif
    enddo


    if(sib%prog%nsl == 0) then

        sib%prog%snow_depth = 0.0_dbl_kind
        sib%prog%snow_mass  = 0.0_dbl_kind

        sib%diag%snow_end(1) = min(sib%diag%snow_end(3),   &
                                     (sib%stat%julday))

        !...set layer values to zero
        do j=-nsnow+1,0
            sib%prog%td(sib%prog%nsl) = 0.0_dbl_kind
            sib%prog%dz(j)      = 0.0_dbl_kind
            sib%prog%layer_z(j) = 0.0_dbl_kind
            sib%prog%node_z(j)  = 0.0_dbl_kind
            sib%prog%www_liq(sib%prog%nsl) = 0.0_dbl_kind
            sib%prog%www_ice(sib%prog%nsl) = 0.0_dbl_kind
        enddo

        !...make sure top layer is zero-ed out also
        sib%prog%layer_z(-nsnow) = 0.0_dbl_kind

        return

    else
        sib%prog%snow_depth = 0.0_dbl_kind
        sib%prog%snow_mass  = 0.0_dbl_kind
        zwice      = 0.0_dbl_kind
        zwliq      = 0.0_dbl_kind

        do j=sib%prog%nsl+1,0
            sib%prog%snow_mass  = sib%prog%snow_mass + sib%prog%www_ice(j) + &
                sib%prog%www_liq(j)
            sib%prog%snow_depth = sib%prog%snow_depth + sib%prog%dz(j)
            zwice      = zwice + sib%prog%www_ice(j)
            zwliq      = zwliq + sib%prog%www_liq(j)
        enddo

        if(sib%prog%snow_mass < 1.0) then
           sib%diag%snow_end(3) = min(sib%diag%snow_end(3),   &
                                     (sib%stat%julday))
        endif

    endif


    !...check the snow depth

    if(sib%prog%snow_depth < 0.01 ) then   ! all snow gone!

        sib%prog%nsl = 0
        sib%prog%snow_mass = zwice
        sib%prog%dz(0) = 0.0_dbl_kind
        sib%prog%node_z(0) = 0.0_dbl_kind
        sib%prog%layer_z(-1) = 0.0_dbl_kind
        if(sib%prog%snow_mass <= 0.0 ) sib%prog%snow_depth = 0.0_dbl_kind

        !...the liquid water assumed ponding on soil surface

        sib%prog%capac(2) = sib%prog%capac(2) + zwliq
        return
    else

        !...two or more layers

        if(sib%prog%nsl < -1) then
            nsnow_old = sib%prog%nsl
            mssi = 1
            do k = nsnow_old+1,0

                !...if top node is removed, combine with bottom neighbor

                if(sib%prog%dz(k) < dzmin(mssi)) then
                    if(k == sib%prog%nsl+1)then
                        neibor = k + 1

                        !...if the bottom neighbor is not snow, 
                        !...combine with the top neighbor

                    elseif(k == 0) then
                        neibor = k - 1

                        !...if neither of the above apply, 
                        !...combine with thinnest neighbor

                    else
                        neibor = k + 1
                        if((sib%prog%dz(k-1) + sib%prog%dz(k)) < &
                            (sib%prog%dz(k+1) + sib%prog%dz(k))) neibor = k-1
                    endif

                    !...node l and j are combined and stored as node j

                    if(neibor > k) then
                        j = neibor
                        l = k
                    else
                        j = k
                        l = neibor
                    endif

                    call clm_combo(sib%prog%dz(j), sib%prog%www_liq(j), &
                        sib%prog%www_ice(j), sib%prog%td(j), sib%prog%dz(l), &
                        sib%prog%www_liq(l), sib%prog%www_ice(l), sib%prog%td(l) )


                    !...now shift all elements above this down one
                    if(j-1 > sib%prog%nsl+1) then
                        do m= j-1, sib%prog%nsl+2, -1
                            sib%prog%td(m)      = sib%prog%td(m-1)
                            sib%prog%www_ice(m) = sib%prog%www_ice(m-1)
                            sib%prog%www_liq(m) = sib%prog%www_liq(m-1)
                            sib%prog%dz(m)      = sib%prog%dz(m-1)
                        enddo
                    endif

                    sib%prog%nsl = sib%prog%nsl + 1
                    if(sib%prog%nsl >= 1 ) cycle

                else

                    mssi = mssi + 1

                endif ! if thickness greater than minimum

            enddo  ! (k) snow layer loop

        endif ! if two or more layers condition

        !...reset the node depth and the depth of layer interface
        do j=0,sib%prog%nsl+1,-1

            sib%prog%node_z(j)    = sib%prog%layer_z(j) - 0.5 * sib%prog%dz(j)
            sib%prog%layer_z(j-1) = sib%prog%layer_z(j) - sib%prog%dz(j)

        enddo

    endif  ! what's the total depth? condition     


end subroutine combine_snow

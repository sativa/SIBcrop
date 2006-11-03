!----------------------------------------------------------------------
subroutine subdivide_snow(sib)
!----------------------------------------------------------------------
!
!   Based on CLM subroutine CLM_SUBDIV
!
!   Description
!   Subdivides snow layers if they exceed their prescribed 
!   maximum thickness
!
!   Revision History:
!   15 September 1999: Yongjiu Dai, initial code
!   15 December  1999: Paul Houser and Jon Radakovich, F90 revision
!   30 January   2002: Ian Baker, SiB integration
!----------------------------------------------------------------------

use kinds
use sibtype

implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  


!...local variables
integer(kind=int_kind) :: msno,j
real(kind=dbl_kind) :: dzsno(5)
real(kind=dbl_kind) :: swice(5)
real(kind=dbl_kind) :: swliq(5)
real(kind=dbl_kind) :: tsnow(5)
real(kind=dbl_kind) :: drr
real(kind=dbl_kind) :: propor
real(kind=dbl_kind) :: zwice
real(kind=dbl_kind) :: zwliq

    !if ( sib%prog%nsl == 0 ) return

    msno = abs(sib%prog%nsl)

    do j=1,msno
        dzsno(j) = sib%prog%dz(j+sib%prog%nsl)
        swice(j) = sib%prog%www_ice(j+sib%prog%nsl)
        swliq(j) = sib%prog%www_liq(j+sib%prog%nsl)
        tsnow(j) = sib%prog%td(j+sib%prog%nsl)
    enddo

    if(msno == 1) then
        if(dzsno(1) > 0.03) then
            msno = 2

            !...specify a new snow layer
            dzsno(1) = dzsno(1)/2.0
            swice(1) = swice(1)/2.0
            swliq(1) = swliq(1)/2.0

            dzsno(2) = dzsno(1)
            swice(2) = swice(1)
            swliq(2) = swliq(1)
            tsnow(2) = tsnow(1)
        endif

    endif   ! if msno == 1 condition


    if(msno > 1) then

        if(dzsno(1) > 0.02 ) then
            drr      = dzsno(1) - 0.02
            propor   = drr/dzsno(1)
            zwice    = propor*swice(1)
            zwliq    = propor*swliq(1)
            propor   = 0.02/dzsno(1)
            swice(1) = propor*swice(1)
            swliq(1) = propor*swliq(1)
            dzsno(1) = 0.02


            call clm_combo(dzsno(2),swliq(2),swice(2),tsnow(2),         &
                drr,zwliq,zwice,tsnow(1))


            if(msno <= 2  .AND. dzsno(2) > 0.07 ) then

                !...subdivide a new layer
                msno = 3
                dzsno(2) = dzsno(2)/2.0
                swice(2) = swice(2)/2.0
                swliq(2) = swliq(2)/2.0
                dzsno(3) = dzsno(2)
                swice(3) = swice(2)
                swliq(3) = swliq(2)
                tsnow(3) = tsnow(2)
            endif
        endif     ! if dzsno(1) > 0.02 condition
    endif       ! if msno > 1 condition


    if(msno > 2) then
        if(dzsno(2) > 0.05) then

            drr      = dzsno(2) - 0.05
            propor   = drr/dzsno(2)
            zwice    = propor*swice(2)
            zwliq    = propor*swliq(2)
            propor   = 0.05/dzsno(2)
            swice(2) = propor*swice(2)
            swliq(2) = propor*swliq(2)
            dzsno(2) = 0.05

            call clm_combo(dzsno(3),swliq(3),swice(3),tsnow(3),         &
                drr,zwliq,zwice,tsnow(2))



            if(msno <= 3  .AND.  dzsno(3) > 0.18) then

                !...subdivide a new layer
                msno = 4
                dzsno(3) = dzsno(3)/2.0
                swice(3) = swice(3)/2.0
                swliq(3) = swliq(3)/2.0
                dzsno(4) = dzsno(3)
                swice(4) = swice(3)
                swliq(4) = swliq(3)
                tsnow(4) = tsnow(3) 
            endif
        endif    ! if dzsno(2) > 0.05 condition
    endif      ! if msno > 2 condition


    if(msno > 3) then
        if(dzsno(3) > 0.11) then

            drr      = dzsno(3) - 0.11
            propor   = drr/dzsno(3)
            zwice    = propor*swice(3)
            zwliq    = propor*swliq(3)
            propor   = 0.11/dzsno(3)
            swice(3) = propor*swice(3)
            swliq(3) = propor*swliq(3)
            dzsno(3) = 0.11

            call clm_combo(dzsno(4),swliq(4),swice(4),tsnow(4),         &
                drr,zwliq,zwice,tsnow(3))


            if(msno <= 4  .AND.  dzsno(4) > 0.41) then

                !...subdivide a new layer
                msno = 5
                dzsno(4) = dzsno(4)/2.0
                swice(4) = swice(4)/2.0
                swliq(4) = swliq(4)/2.0
                dzsno(5) = dzsno(4)
                swice(5) = swice(4)
                swliq(5) = swliq(4)
                tsnow(5) = tsnow(4)
            endif
        endif    ! if dzsno(3) > 0.11 condition
    endif      ! if msno > 3 condition


    if(msno > 4) then
        if(dzsno(4) > 0.23) then
            drr      = dzsno(4) - 0.23
            propor   = drr/dzsno(4)
            zwice    = propor*swice(4)
            zwliq    = propor*swliq(4)
            propor   = 0.23/dzsno(4)
            swice(4) = propor*swice(4)
            swliq(4) = propor*swliq(4)
            dzsno(4) = 0.23

            call clm_combo(dzsno(5),swliq(5),swice(5),tsnow(5),         &
                drr,zwliq,zwice,tsnow(4))



        endif    ! if dzsno(4) > 0.23 condition
    endif      ! if msno > 4 condition

    sib%prog%nsl = -msno

    do j=sib%prog%nsl+1,0
        sib%prog%dz(j) = dzsno(j - sib%prog%nsl)
        sib%prog%www_ice(j) = swice(j - sib%prog%nsl)
        sib%prog%www_liq(j) = swliq(j - sib%prog%nsl)
        sib%prog%td(j)      = tsnow(j - sib%prog%nsl)
    enddo

    do j=0,sib%prog%nsl+1,-1

        sib%prog%node_z(j) = sib%prog%layer_z(j) - 0.5 * sib%prog%dz(j)
        sib%prog%layer_z(j-1) = sib%prog%node_z(j) - 0.5*sib%prog%dz(j)

    enddo


end subroutine subdivide_snow

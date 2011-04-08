
!===================SUBROUTINE SORTIN===================================

subroutine sortin( eyy, pco2y, range, gammas, ic)

!=======================================================================
!
!     ARRANGES SUCCESSIVE PCO2/ERROR PAIRS IN ORDER OF INCREASING PCO2.
!       ESTIMATES NEXT GUESS FOR PCO2 USING COMBINATION OF LINEAR AND
!       QUADRATIC FITS.
!
!=======================================================================


    use kinds

    implicit none


    !Bio...INPUT VARIABLES

    integer(kind=int_kind),intent(in) :: ic     ! iteration count

    real(kind=dbl_kind),intent(in) :: range     !
    real(kind=dbl_kind),intent(in) :: gammas    !


    !Bio...OUTPUT VARIABLES...
    real(kind=dbl_kind),intent(inout),dimension(6) :: eyy    !
    real(kind=dbl_kind),intent(inout),dimension(6) :: pco2y  !

    !Bio...LOCAL VARIABLES

    integer(kind=int_kind) ::  i,j,n,i1,i2,i3,isp,is,ix

    logical (kind=log_kind) :: bitx    !

    real(kind=dbl_kind) :: one         !
    real(kind=dbl_kind) :: pmin        !
    real(kind=dbl_kind) :: emin        !
    real(kind=dbl_kind) :: a           !
    real(kind=dbl_kind) :: b           !
    real(kind=dbl_kind) :: pco2yl      !
    real(kind=dbl_kind) :: pco2yq      !
    real(kind=dbl_kind) :: ac1         !
    real(kind=dbl_kind) :: ac2         !
    real(kind=dbl_kind) :: bc1         !
    real(kind=dbl_kind) :: bc2         !
    real(kind=dbl_kind) :: cc1         !
    real(kind=dbl_kind) :: cc2         !
    real(kind=dbl_kind) :: aterm       !
    real(kind=dbl_kind) :: bterm       !
    real(kind=dbl_kind) :: cterm       !
    real(kind=dbl_kind) :: pco2b       !
    real(kind=dbl_kind) :: eyyisp      !
    real(kind=dbl_kind) :: eyyis       !
    real(kind=dbl_kind) :: eyyi1       !
    real(kind=dbl_kind) :: eyyi2       !
    real(kind=dbl_kind) :: eyyi3       !
    real(kind=dbl_kind) :: pco2yisp    !
    real(kind=dbl_kind) :: pco2yis     !
    real(kind=dbl_kind) :: pco2yi1     !
    real(kind=dbl_kind) :: pco2yi2     !
    real(kind=dbl_kind) :: pco2yi3     !


    one = 1.0_dbl_kind

    if( ic < 4 ) then
        pco2y(1) = gammas + 0.5_dbl_kind*range
        pco2y(2) = gammas                                             &
            + range*( 0.5_dbl_kind - 0.3_dbl_kind*sign(one,eyy(1)) )
        pco2y(3) = pco2y(1)- (pco2y(1)-pco2y(2))                      &
            /(eyy(1)-eyy(2)+1.e-10_dbl_kind)*eyy(1)
        pmin = min( pco2y(1), pco2y(2) )
        emin = min(   eyy(1),   eyy(2) )
        if ( emin > 0. .and. pco2y(3) > pmin )                        &
            pco2y(3) = gammas
    else

        n = ic - 1
        bitx = abs(eyy(n)) > 0.1
        if(.not. bitx) pco2y(ic) = pco2y(n)
        if(bitx) then
            do j = 2, n
                a = eyy(j)
                b = pco2y(j)
                do i = j-1,1,-1
                    if(eyy(i) <= a ) go to 100
                    eyy(i+1) = eyy(i)
                    pco2y(i+1) = pco2y(i)
                enddo ! i loop
                i = 0
                100        continue
                eyy(i+1) = a
                pco2y(i+1) = b
            enddo  ! j loop
        endif

!-----------------------------------------------------------------------

        if(bitx) then
            pco2b = 0.
            is    = 1
        endif

        do ix = 1, n
            if(bitx) then
                if( eyy(ix) < 0. )  then
                    pco2b = pco2y(ix)
                    is = ix
                endif
            endif
        enddo

        if(bitx) then
            i1 = is-1
            i1 = MAX(1, i1)
            i1 = min(n-2, i1)
            i2 = i1 + 1
            i3 = i1 + 2
            isp   = is + 1
            isp = min0( isp, n )
            is = isp - 1
            eyyisp = eyy(isp)
            eyyis = eyy(is)
            eyyi1 = eyy(i1)
            eyyi2 = eyy(i2)
            eyyi3 = eyy(i3)
            pco2yisp = pco2y(isp)
            pco2yis = pco2y(is)
            pco2yi1 = pco2y(i1)
            pco2yi2 = pco2y(i2)
            pco2yi3 = pco2y(i3)
        endif

        if(bitx) then

            !itb...Neil Suits' patch to check for zero in the denominator...
            if(eyyis /= eyyisp)then
                pco2yl=pco2yis - (pco2yis-pco2yisp) / (eyyis-eyyisp)*eyyis
            else
                pco2yl = pco2yis * 1.01
            endif

            !   METHOD USING A QUADRATIC FIT

            ac1 = eyyi1*eyyi1 - eyyi2*eyyi2
            ac2 = eyyi2*eyyi2 - eyyi3*eyyi3
            bc1 = eyyi1 - eyyi2
            bc2 = eyyi2 - eyyi3
            cc1 = pco2yi1 - pco2yi2
            cc2 = pco2yi2 - pco2yi3

            !itb...Neil Suits' patch to prevent zero in denominator...
            if(bc1*ac2-ac1*bc2 /= 0.0 .and. ac1 /= 0.0_dbl_kind)then
                bterm = (cc1*ac2-cc2*ac1)/(bc1*ac2-ac1*bc2)
                aterm = (cc1-bc1*bterm)/ac1
                cterm = pco2yi2-aterm*eyyi2*eyyi2-bterm*eyyi2
                pco2yq= cterm
                pco2yq= MAX( pco2yq, pco2b )
                pco2y(ic) = ( pco2yl+pco2yq)/2.0_dbl_kind
            else
                pco2y(ic) = pco2y(ic) * 1.01_dbl_kind
            endif

        endif

    endif
!
! make sure pco2 does not fall below compensation point
    pco2y(ic) = MAX(pco2y(ic),gammas+0.01_dbl_kind)

end

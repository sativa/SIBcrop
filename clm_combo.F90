!----------------------------------------------------------------------
subroutine CLM_COMBO(dz1,liq1,ice1,temp1,dz2,liq2,ice2,temp2)

    use kinds
    use eau_params, only: &
        lfus
    use physical_parameters, only: &
        tice
    use sib_const_module, only: &
        cpice, &
        cpliq 

    implicit none
    real(kind=dbl_kind),intent(inout) :: dz1
    real(kind=dbl_kind),intent(inout) :: liq1
    real(kind=dbl_kind),intent(inout) :: ice1
    real(kind=dbl_kind),intent(inout) :: temp1

    real(kind=dbl_kind),intent(in) :: dz2
    real(kind=dbl_kind),intent(in) :: liq2
    real(kind=dbl_kind),intent(in) :: ice2
    real(kind=dbl_kind),intent(in) :: temp2


    !...local variables...
    real(kind=dbl_kind) :: dzc
    real(kind=dbl_kind) :: wicec
    real(kind=dbl_kind) :: wliqc
    real(kind=dbl_kind) :: tc
    real(kind=dbl_kind) :: h1
    real(kind=dbl_kind) :: h2
    real(kind=dbl_kind) :: hc

    !
    !   Code based on CLM subroutine CLM_COMBO, modified for use
    !   with SiB
    !
    !   CLM Web Info:  http://clm.gsfc.nasa.gov
    !
    !   Description
    !   Combines two elements, and returns dz (thickness) temperature
    !   www_liq and www_ice.
    !
    !   Revision History:
    !   15 September 1999; Yongjiu Dai, original code
    !   15 December 1999;  Paul Houser and Jon Radakovich, F90 revision
    !   01 March 2002;     Ian Baker, SiB integration

    dzc   = dz1 + dz2
    wicec = ice1 + ice2
    wliqc = liq1 + liq2

    h1    = (cpice*ice1+cpliq*liq1) &
        *(temp1-tice)+lfus*liq1
    h2    = (cpice*ice2+cpliq*liq2) &
        *(temp2-tice)+lfus*liq2
    hc    = h1 + h2

    if(hc < 0.0) then
        tc = tice + hc/(cpice*wicec+cpliq*wliqc)
    elseif(hc <= lfus*wliqc) then
        tc = tice
    else
        tc = tice + (hc - lfus*wliqc)/(cpice*wicec &
            +cpliq*wliqc)
    endif

    dz1   = dzc
    ice1  = wicec
    liq1  = wliqc
    temp1 = tc

    return
end subroutine clm_combo

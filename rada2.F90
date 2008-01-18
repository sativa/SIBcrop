!=================SUBROUTINE RADA2======================================
subroutine rada2(sib,sib_loc)



    !=======================================================================
    !
    !     CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION( DIRECT
    !     AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY
    !
    !-----------------------------------------------------------------------


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       SALB(2,2)      SURFACE ALBEDOS 
    !       TGEFF4         EFFECTIVE SURFACE RADIATIVE TEMPERATURE (K) 
    !       RADFAC(2,2,2)  RADIATION ABSORPTION FACTORS 
    !       THERMK         CANOPY GAP FRACTION FOR TIR RADIATION 
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       ALBEDO(2,2,2)  COMPONENT REFLECTANCES 
    !
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !      REFERENCE
    !
    !      Sellers, P.J., 1985: Canopy Reflectance, Photosynthesis and
    !                     Respiration. Int. J. Remote Sensing, 6(8)
    !                     1335-1372.



    use kinds
    use sibtype
    use sib_const_module, only:  &
        zlnd
    use physical_parameters, only :  &
        tice
    implicit none   

    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib
    type(sib_local_vars)     ,intent(inout) :: sib_loc
    ! variables local to SiB

    !----------------------------------------------------------------------  

    !...LOCAL VARIABLES...
    integer(kind=int_kind) :: iwave, irad   ! loop variables

    real(kind=dbl_kind) :: fff        ! looks like a minimum cosz value
    real(kind=dbl_kind) :: facs       ! 1/20th of Td(1) - Tice, constrained
    !  to be between 0.0 and 0.4  
    real(kind=dbl_kind) :: fmelt      ! 1-facs
    real(kind=dbl_kind) :: scov
    real(kind=dbl_kind) :: reff1
    real(kind=dbl_kind) :: reff2
    real(kind=dbl_kind) :: tran1
    real(kind=dbl_kind) :: tran2
    real(kind=dbl_kind) :: scat
    real(kind=dbl_kind) :: chiv      ! copy of sib_param%chil
    real(kind=dbl_kind) :: aa
    real(kind=dbl_kind) :: bb
    real(kind=dbl_kind) :: proj
    real(kind=dbl_kind) :: extkb
    real(kind=dbl_kind) :: zmew
    real(kind=dbl_kind) :: acss
    real(kind=dbl_kind) :: upscat
    real(kind=dbl_kind) :: betao
    real(kind=dbl_kind) :: be
    real(kind=dbl_kind) :: ce
    real(kind=dbl_kind) :: bot
    real(kind=dbl_kind) :: de
    real(kind=dbl_kind) :: fe
    real(kind=dbl_kind) :: ge
    real(kind=dbl_kind) :: hh1
    real(kind=dbl_kind) :: hh2
    real(kind=dbl_kind) :: hh3
    real(kind=dbl_kind) :: hh4
    real(kind=dbl_kind) :: hh5
    real(kind=dbl_kind) :: hh6
    real(kind=dbl_kind) :: hh7
    real(kind=dbl_kind) :: hh8
    real(kind=dbl_kind) :: hh9
    real(kind=dbl_kind) :: hh10
    real(kind=dbl_kind) :: psi
    real(kind=dbl_kind) :: zat
    real(kind=dbl_kind) :: power1
    real(kind=dbl_kind) :: power2
    real(kind=dbl_kind) :: epsi
    real(kind=dbl_kind) :: ek
    real(kind=dbl_kind) :: albedo(2,2,2)
    real(kind=dbl_kind) :: f1
    real(kind=dbl_kind) :: zp
    real(kind=dbl_kind) :: den
    real(kind=dbl_kind) :: zmk
    real(kind=dbl_kind) :: tranc1(2)
    real(kind=dbl_kind) :: tranc2(2)
    real(kind=dbl_kind) :: tranc3(2)
    real(kind=dbl_kind) :: tsurf
    real(kind=dbl_kind) :: tg4
    real(kind=dbl_kind) :: tc4
    real(kind=dbl_kind) :: fac2
    real(kind=dbl_kind) :: zkat


    !
    !----------------------------------------------------------------------
    !
    !
    !     MODIFICATION FOR EFFECT OF SNOW ON UPPER STOREY ALBEDO
    !         SNOW REFLECTANCE   = 0.80, 0.40 . MULTIPLY BY 0.6 IF MELTING
    !         SNOW TRANSMITTANCE = 0.20, 0.54
    !
    !
    !-----------------------------------------------------------------------
    !
    !

    sib%diag%canex         = 1.-( sib%prog%snow_veg*5.-sib%param%z1)/  &
        (sib%param%z2-sib%param%z1)
    sib%diag%canex         = max( 0.1_dbl_kind, sib%diag%canex )
    sib%diag%canex         = min( 1.0_dbl_kind, sib%diag%canex )


    !...both satcap values are in meters: multiply by density to get
    !...kg/m^2...
    sib%param%satcap(1) = sib%param%zlt * 0.0001 * sib%diag%canex
!    sib%param%satcap(2) = 0.0002           ! lahouari
    sib%param%satcap(2) = 0.01           ! Baker

    sib%diag%areas = sib%prog%snow_depth / (zlnd*10.0 + sib%prog%snow_depth)

    !itb...areas criteria for vanishing snow
    if(sib%diag%areas < 0.25 .and. sib%stat%julday > 3 ) then
      sib%diag%snow_end(2) = MIN(sib%diag%snow_end(2),(sib%stat%julday))
    endif

    fff = max(0.01746_dbl_kind,sib%stat%cosz)

    !itb...facs only accounts for ground sfc temperature-no snow influence 
    !      facs  = ( sib%prog%td(1) - tice ) * 0.04

    !itb...so we'll change it...
    facs = (sib%prog%td(sib%prog%nsl+1) - tice) * 0.04

    facs  = max( 0.0_dbl_kind , facs)
    facs  = min( 0.4_dbl_kind, facs)
    fmelt = 1.0 - facs

    !-----------------------------------------------------------------------
    do iwave = 1, 2
!print*,sib%prog%snow_veg/1000.0,sib%param%satcap(1)
        scov =  min( 0.5_dbl_kind, (sib%prog%snow_veg/1000.0)/sib%param%satcap(1) )

        reff1 = ( 1. - scov ) * sib%param%ref(iwave,1) + scov * ( 1.2 -     &
            iwave * 0.4 ) * fmelt
        reff2 = ( 1. - scov ) * sib%param%ref(iwave,2) + scov * ( 1.2 -     &
            iwave * 0.4 ) * fmelt
        tran1 = sib%param%tran(iwave,1) * ( 1. - scov )                     &
            + scov * ( 1.- ( 1.2 - iwave * 0.4 ) * fmelt )         &
            * sib%param%tran(iwave,1)
        tran2 = sib%param%tran(iwave,2) * ( 1. - scov )                     &
            + scov * ( 1.- ( 1.2 - iwave * 0.4 ) * fmelt ) * 0.9   &
            * sib%param%tran(iwave,2)
        !-----------------------------------------------------------------------
        !
        !     CALCULATE AVERAGE SCATTERING COEFFICIENT, LEAF PROJECTION AND
        !     OTHER COEFFICIENTS FOR TWO-STREAM MODEL.
        !
        !      SCAT  (OMEGA)         : EQUATION (1,2) , SE-85
        !      PROJ  (G(MU))         : EQUATION (13)  , SE-85
        !      EXTKB (K, G(MU)/MU)   : EQUATION (1,2) , SE-85
        !      ZMEW  (INT(MU/G(MU))  : EQUATION (1,2) , SE-85
        !      ACSS  (A-S(MU))       : EQUATION (5)   , SE-85
        !      EXTK  (K, VARIOUS)    : EQUATION (13)  , SE-85
        !      UPSCAT(OMEGA*BETA)    : EQUATION (3)   , SE-85
        !      BETAO (1BETA-0)       : EQUATION (4)   , SE-85 
        !
        !-----------------------------------------------------------------------

        scat = sib%param%green*( tran1 + reff1 ) +( 1.-sib%param%green ) *  &
            ( tran2 + reff2)
        chiv = sib%param%chil

        if ( abs(chiv) .LE. 0.01 ) chiv = 0.01

        aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
        bb = 0.877 * ( 1. - 2. * aa )

        proj = aa + bb * fff
        extkb = ( aa + bb * fff ) / fff
        zmew = 1. / bb * ( 1. - aa / bb   &
            * log ( ( aa + bb ) / aa ) )
        acss = scat / 2. * proj / ( proj + fff * bb )
        acss = acss * ( 1. - fff * aa     &
            / ( proj + fff * bb ) * log ( ( proj   &
            +   fff * bb + fff * aa ) / ( fff * aa ) ) )

        upscat = sib%param%green * tran1 + ( 1.- sib%param%green ) * tran2
        upscat = 0.5 * ( scat + ( scat - 2. * upscat ) *   &
            (( 1. - chiv ) / 2. ) ** 2 )
        betao = ( 1. + zmew * extkb )   &
            / ( scat * zmew * extkb ) * acss

        !-----------------------------------------------------------------------
        !
        !     Intermediate variables identified in appendix of SE-85.
        !
        !      BE          (B)     : APPENDIX      , SE-85
        !      CE          (C)     : APPENDIX      , SE-85
        !      BOT         (SIGMA) : APPENDIX      , SE-85
        !      HH1         (H1)    : APPENDIX      , SE-85
        !      HH2         (H2)    : APPENDIX      , SE-85
        !      HH3         (H3)    : APPENDIX      , SE-85
        !      HH4         (H4)    : APPENDIX      , SE-85
        !      HH5         (H5)    : APPENDIX      , SE-85
        !      HH6         (H6)    : APPENDIX      , SE-85
        !      HH7         (H7)    : APPENDIX      , SE-85
        !      HH8         (H8)    : APPENDIX      , SE-85
        !      HH9         (H9)    : APPENDIX      , SE-85
        !      HH10        (H10)   : APPENDIX      , SE-85
        !      PSI         (H)     : APPENDIX      , SE-85
        !      ZAT         (L-T)   : APPENDIX      , SE-85
        !      EPSI        (S1)    : APPENDIX      , SE-85
        !      EK          (S2)    : APPENDIX      , SE-85
        !-----------------------------------------------------------------------

        be = 1. - scat + upscat
        ce = upscat
        bot = ( zmew * extkb ) ** 2 + ( ce**2 - be**2 )

        if ( abs(bot) <= 1.e-10) then
            scat = scat* 0.98
            be = 1. - scat + upscat
            bot = ( zmew * extkb ) ** 2 + ( ce**2 - be**2 )
        endif

        de = scat * zmew * extkb * betao
        fe = scat * zmew * extkb * ( 1. - betao )
        hh1 = -de * be + zmew * de * extkb - ce * fe
        hh4 = -be * fe - zmew * fe * extkb - ce * de

        psi = sqrt(be**2 - ce**2)/zmew

        zat = sib%param%zlt/sib%param%vcover*sib%diag%canex

        power1 = min( psi*zat, 50.0_dbl_kind )
        power2 = min( extkb*zat, 50.0_dbl_kind )
        epsi = exp( - power1 )
        ek = exp ( - power2 )

        albedo(2,iwave,1) = sib%param%soref(iwave)*(1.-sib%diag%areas)  &
            + ( 1.2-iwave*0.4 )*fmelt * sib%diag%areas
        albedo(2,iwave,2) = sib%param%soref(iwave)*(1.-sib%diag%areas)  &
            + ( 1.2-iwave*0.4 )*fmelt * sib%diag%areas
        ge = albedo(2,iwave,1)/albedo(2,iwave,2)

        !----------------------------------------------------------------
        !     CALCULATION OF DIFFUSE ALBEDOS
        !
        !     ALBEDO(1,IWAVE,2) ( I-UP ) : APPENDIX , SE-85
        !----------------------------------------------------------------

        f1 = be - ce / albedo(2,iwave,2)
        zp = zmew * psi

        den = ( be + zp ) * ( f1 - zp ) / epsi -   &
            ( be - zp ) * ( f1 + zp ) * epsi
        hh7 = ce * ( f1 - zp ) / epsi / den
        hh8 = -ce * ( f1 + zp ) * epsi / den
        f1 = be - ce * albedo(2,iwave,2)
        den = ( f1 + zp ) / epsi - ( f1 - zp ) * epsi

        hh9 = ( f1 + zp ) / epsi / den
        hh10 = - ( f1 - zp ) * epsi / den
        tranc2(iwave) = hh9 * epsi + hh10 / epsi

        albedo(1,iwave,2) =  hh7 + hh8

        !-----------------------------------------------------------------
        !     CALCULATION OF DIRECT ALBEDOS AND CANOPY TRANSMITTANCES.
        !
        !      ALBEDO(1,IWAVE,1) ( I-UP )   : EQUATION(11)   , SE-85
        !      TRANC(IWAVE)      ( I-DOWN ) : EQUATION(10)   , SE-85
        !
        !-----------------------------------------------------------------
        f1 = be - ce / albedo(2,iwave,2)
        zmk = zmew * extkb

        den = ( be + zp ) * ( f1 - zp ) / epsi -     &
            ( be - zp ) * ( f1 + zp ) * epsi
        hh2 = ( de - hh1 / bot * ( be + zmk ) )      &
            * ( f1 - zp ) / epsi -                   &
            ( be - zp ) * ( de - ce*ge - hh1 / bot   &
            * ( f1 + zmk ) ) * ek
        hh2 = hh2 / den
        hh3 = ( be + zp ) * (de - ce * ge -          &
            hh1 / bot * ( f1 + zmk )) * ek -         &
            ( de - hh1 / bot * ( be + zmk ) ) *      &
            ( f1 + zp ) * epsi
        hh3 = hh3 / den
        f1 = be - ce * albedo(2,iwave,2)
        den = ( f1 + zp ) / epsi - ( f1 - zp ) * epsi
        hh5 = - hh4 / bot * ( f1 + zp ) / epsi -     &
            ( fe + ce*ge*albedo(2,iwave,2) +         &
            hh4 / bot * ( zmk - f1 ) ) * ek
        hh5 = hh5 / den
        hh6 =   hh4 / bot * ( f1 - zp ) * epsi +     &
            ( fe + ce * ge * albedo(2,iwave,2) +     &
            hh4 / bot*( zmk - f1 ) ) * ek
        hh6 = hh6 / den
        tranc1(iwave) = ek
        tranc3(iwave) = hh4 / bot * ek + hh5 * epsi + hh6 / epsi

        albedo(1,iwave,1) = hh1 / bot + hh2 + hh3
        !
        !----------------------------------------------------------------------
        !
        !
        !----------------------------------------------------------------------
        !     CALCULATION OF TERMS WHICH MULTIPLY INCOMING SHORT WAVE FLUXES
        !     TO GIVE ABSORPTION OF RADIATION BY CANOPY AND GROUND
        !
        !      RADFAC      (F(IL,IMU,IV)) : EQUATION (19,20) , SE-86
        !----------------------------------------------------------------------
        !
        sib%diag%radfac(2,iwave,1) = ( 1.-sib%param%vcover )   &
            * ( 1.-albedo(2,iwave,1) ) + sib%param%vcover      &
            * ( tranc1(iwave) * ( 1.-albedo(2,iwave,1) )    &
            + tranc3(iwave) * ( 1.-albedo(2,iwave,2) ) )

        sib%diag%radfac(2,iwave,2) = ( 1.-sib%param%vcover )   &
            * ( 1.-albedo(2,iwave,2) ) + sib%param%vcover      &
            *  tranc2(iwave) * ( 1.-albedo(2,iwave,2) )

        sib%diag%radfac(1,iwave,1) = sib%param%vcover          &
            * ( ( 1.-albedo(1,iwave,1) )                    &
            - tranc1(iwave) * ( 1.-albedo(2,iwave,1) )      &
            - tranc3(iwave) * ( 1.-albedo(2,iwave,2) ) )

        sib%diag%radfac(1,iwave,2) = sib%param%vcover          &
            * ( ( 1.-albedo(1,iwave,2) )                    &
            - tranc2(iwave) * ( 1.-albedo(2,iwave,2) ) )
        !
        !----------------------------------------------------------------------
        !     CALCULATION OF TOTAL SURFACE ALBEDOS ( SALB ) WITH WEIGHTING
        !     FOR COVER FRACTIONS.
        !----------------------------------------------------------------------
        !
        do irad = 1, 2
            sib%diag%salb(iwave,irad) = ( 1.-sib%param%vcover )   & 
                * albedo(2,iwave,irad) +                       &
                sib%param%vcover * albedo(1,iwave,irad)
        enddo
        !
        !----------------------------------------------------------------------
        !
    enddo  ! iwave loop
    !
    !----------------------------------------------------------------------
    !
    !     CALCULATION OF LONG-WAVE FLUX TERMS FROM CANOPY AND GROUND
    !
    !----------------------------------------------------------------------
    !
    tsurf = min(tice,sib%prog%td(sib%prog%nsl+1))*sib%diag%areas   &
        + sib%prog%td(1)*(1.-sib%diag%areas)
    tc4 = sib%prog%tc*sib%prog%tc*sib%prog%tc*sib%prog%tc
    tg4 = tsurf*tsurf*tsurf*tsurf   

    zkat = 1./zmew * sib%param%zlt / sib%param%vcover
    zkat = min( 50.0_dbl_kind , zkat )
    zkat = max( 1.E-5_dbl_kind, zkat )
    sib%diag%thermk = exp(-zkat)

    sib_loc%fac1 =  sib%param%vcover * ( 1.-sib%diag%thermk )

end subroutine rada2

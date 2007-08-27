!==================SUBROUTINE PHOSIB===================================
subroutine phosib(sib,sib_loc)

    use kinds
    use sibtype
    use cfrax
    use sib_const_module, only: &
        po2m, &
        pdb,  &
        dtt,  &
        dti
    use physical_parameters, only: &
        p0 => p0_sfc, &
        tice

    implicit none


    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    type(sib_local_vars)     ,intent(inout) :: sib_loc
    ! variables local to SiB

    !----------------------------------------------------------------------  



    !=======================================================================
    !
    !     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED   
    !     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
    !     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
    !     BASE REFERENCE IS SE-92A
    !
    !                          UNITS
    !                         -------
    !
    !      PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
    !      CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
    !      VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
    !      EFFCON                                   : MOL CO2 MOL QUANTA-1
    !      GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
    !      EVAPKG                                   : KG M-2 S-1
    !      Q                                        : KG KG-1
    !
    !                       CONVERSIONS
    !                      -------------
    !
    !      1 MOL H2O           = 0.018 KG
    !      1 MOL CO2           = 0.044 KG
    !      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
    !      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
    !pl the next line applies to the Ci to Cs pathway
    !      GS  (CO2)           = GS (H2O) * 1./1.6
    !pl 44.6 is the number of moles of air per cubic meter
    !      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO
    !      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
    !      MM  (MOLAIR/MOLH2O) = 1.611
    !
    !
    !                         OUTPUT
    !                      -------------
    !
    !      ASSIMN              = CANOPY NET ASSIMILATION RATE
    !      EA                  = CANOPY AIR SPACE VAPOR PRESSURE
    !      1/RST               = CANOPY CONDUCTANCE
    !      PCO2I               = INTERNAL CO2 CONCENTRATION
    !      RESPC               = CANOPY RESPIRATION
    !      RESPG               = GROUND RESPIRATION
    !
    !----------------------------------------------------------------------
    !
    !         RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
    !         RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
    !         RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
    !         RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
    !
    !-----------------------------------------------------------------------


    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1) 
    !       RST            CANOPY RESISTANCE (S M-1)
    !       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS 
    !
    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
    !
    !       RESPC          CANOPY RESPIRATION (MOL M-2 S-1)
    !       RESPG          GROUND RESPIRATION (MOL M-2 S-1)
    !       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
    !       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
    !       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !     Modifications:
    !       - gs (stomatal conductance reduced for freezing soils per Jim Collatz
    !         dd 950221      
    !
    !      Modified for multitasking - introduced gather/scatter indices
    !          - DD 951206
    !
    !itb   Added in pco2c (chloroplast partial co2) for neil's fractionation
    !itb   calculations
    !itb       - IB Sep99


    !Bio...LOCAL VARIABLES
    integer ::  icconv    !

    integer :: &
        ic,  & ! iteration loop
        ic1, & !
        i      ! loop variable

    real(kind=dbl_kind) :: co2cap    ! conversion factor from ppm to mole/m2 in the CAS
    real(kind=dbl_kind) :: c3        ! C3 flag
    real(kind=dbl_kind) :: c4        ! C4 flag
    real(kind=dbl_kind) :: scatp     ! 
    real(kind=dbl_kind) :: scatg     !
    real(kind=dbl_kind) :: park      !
    real(kind=dbl_kind) :: qt        !
    real(kind=dbl_kind) :: respn     !
    real(kind=dbl_kind) :: vm        !
    real(kind=dbl_kind) :: templ     !
    real(kind=dbl_kind) :: temph     !
    real(kind=dbl_kind) :: zkc       !
    real(kind=dbl_kind) :: zko       !
    real(kind=dbl_kind) :: spfy      !
    real(kind=dbl_kind) :: gammas    !
    real(kind=dbl_kind) :: gah2o     !
    real(kind=dbl_kind) :: gbh2o     !
    real(kind=dbl_kind) :: gsh2o     !
    real(kind=dbl_kind) :: xgco2m    !
    real(kind=dbl_kind) :: rrkk      !
    real(kind=dbl_kind) :: soilfrz   !
    real(kind=dbl_kind) :: soilfrztg !
    real(kind=dbl_kind) :: soilfrztd !
    real(kind=dbl_kind) :: bintc     !
    real(kind=dbl_kind) :: omss      !
    real(kind=dbl_kind) :: range1     !
    real(kind=dbl_kind) :: par       !
    real(kind=dbl_kind) :: co2a      ! CAS CO2 concentration (mol C/mol air)
    real(kind=dbl_kind) :: co2m      ! reference level CO2 concentration 
    !   (mol C/mol air)
    real(kind=dbl_kind) :: co2s      ! leaf surface CO2 concentration 
    !   (mol C/mol air)
    real(kind=dbl_kind) :: pco2a     ! intermediate CAS CO2 concentration 
    !   (Pa)
    real(kind=dbl_kind) :: pco2s     ! intermediate leaf surface CO2 
    !   partial pressure (Pa)
    real(kind=dbl_kind) :: pco2i     ! intermediate stomatal (internal) 
    !   CO2 partial pressure (Pa)
    real(kind=dbl_kind) :: pco2c     ! intermediate leaf chloroplast CO2 
    !   partial pressure (Pa)
    real(kind=dbl_kind) :: h2oi      !
    real(kind=dbl_kind) :: h2oa      !
    real(kind=dbl_kind) :: h2os      !
    real(kind=dbl_kind) :: h2osrh    !
    real(kind=dbl_kind) :: ecmole    !
    real(kind=dbl_kind) :: pco2y(6)  !
    real(kind=dbl_kind) :: eyy(6)    !
    real(kind=dbl_kind) :: assimny(6)!
    real(kind=dbl_kind) :: assimy(6) !
    real(kind=dbl_kind) :: gsh2oinf  !
    real(kind=dbl_kind) :: drst(5)   ! delta of stomatal resistance (sec/m)
    real(kind=dbl_kind) :: pco2ipot  !
    real(kind=dbl_kind) :: omcpot    !
    real(kind=dbl_kind) :: sqrtin    !
    real(kind=dbl_kind) :: omppot    !
    real(kind=dbl_kind) :: omspot    !
    real(kind=dbl_kind) :: omcci     !
    real(kind=dbl_kind) :: ompci     !
    real(kind=dbl_kind) :: omsci     !
    real(kind=dbl_kind) :: dompdomc  !
    real(kind=dbl_kind) :: ascitemp  !
    real(kind=dbl_kind) :: ccomc     !
    real(kind=dbl_kind) :: ccoms     !
    real(kind=dbl_kind) :: cwsfws    !
    real(kind=dbl_kind) :: cwsfht    !
    real(kind=dbl_kind) :: cwsflt    !
    real(kind=dbl_kind) :: rstfac3(6)! intermediate temperature stress factor
    real(kind=dbl_kind) :: zln2      ! used in calculating pdamp,qdamp
    real(kind=dbl_kind) :: ghalf     ! used in calculating pdamp,qdamp
    real(kind=dbl_kind) :: dttin     ! used in calculating pdamp,qdamp
    real(kind=dbl_kind) :: dmin      ! used in calculating pdamp,qdamp
    real(kind=dbl_kind) :: pdamp
    real(kind=dbl_kind) :: qdamp
    real(kind=dbl_kind) :: tprcor    ! temperature correction (K)


!itb...playing with carbon mass balance...
    real(kind=dbl_kind) :: co2a_star ! intermediate value of co2a to be used with
                                      ! time filter when conditions warrant

    real(kind=dbl_kind) :: gah2o_crit ! critical value of CAS-ref level conductance. 
                                       ! Conductance above this value invokes time filter.
    real(kind=dbl_kind) :: rstar    ! universal gas constant (N m mole^-1 K^-1)
    
    rstar = 8.3143


    !pl introduce a co2 capacity 
    !pl this will basically be the mass of air under the top of the canopy (in
    !pl this case (CHEAS-RAMS) O(10-30m), that is, ground to displacemnt height.

    !pl all the carbon fluxes are expresse as Mol C / m2 s and resistances for
    !pl carbon are in m2 s / mol air

    !pl one mole of gas occupies 22.4 cubic dm
    !pl 1 cubic meter contains therefore 1000./22.4  = 44.6 moles of gas
    !pl the units for the carbon capacity are mol air /m2. 
    !pl (e.g. here 893 moles if thickness of the layer is 20m)
    !pl this means that the units for pc02ap should be mol co2 / mol air, but
    !pl it is also possible to keep just co2 pressure and convert

    !
    !   calculate damping factors
    !
    zln2 = 6.9314718e-1
    ghalf = 1.0257068e1
    dttin = 3.6e3 
    dmin = 6.0e1
    pdamp = exp (-1.0 * zln2*(dtt*dmin)/(dttin*ghalf))
    qdamp = 1.0 - pdamp
    tprcor = tice*sib%prog%ps*100.0/p0
    co2cap = sib%diag%cas_cap_co2 * 44.6 * tprcor/sib%prog%ta  ! moles air / m2
    co2cap = sib%diag%cas_cap_co2 * sib%prog%ps*100.0 /rstar/sib%prog%ta


    !-----------------------------------------------------------------------
    !
    !     CALCULATION OF CANOPY PAR USE PARAMETER.
    !
    !      APARKK      (PI)     : EQUATION (31) , SE-92A
    !-----------------------------------------------------------------------

    scatp =     sib%param%green  *   (  sib%param%tran(1,1) +  sib%param%ref(1,1) )   &
   +  ( 1.- sib%param%green ) *  (  sib%param%tran(1,2) +  sib%param%ref(1,2) )

    scatg =  sib%param%tran(1,1) +  sib%param%ref(1,1)

    park = sqrt(1.-scatp) *  sib%param%gmudmu

    !itb...Niall integrates physfrac into aparkk. I'm not sure I like
    !itb...doing it that way--SO I WON'T, FOR NOW...

    sib%diag%aparkk   = sib%param%aparc / park * sib%param%green






    !itb...start PHYSIOLOGY LOOP here...
    !itb...potentially 5 different physiologies can share the same
    !itb...soil and CAS (making it different from a normal tile).
    !itb...loop will cycle out of an unused physiology type.

    !
    ! zero out physiology-specific values
    !
    sib%prog%rst(6)     = 0.0

    
        sib%diag%pco2c(6)   = 0.0
        sib%diag%pco2i(6)   = 0.0
        sib%diag%pco2s(6)   = 0.0
        sib%diag%assimn(6)  = 0.0
        sib%diag%assim(6)   = 0.0
        sib%diag%ggl(6)     = 0.0
        sib%diag%antemp(6)  = 0.0
        sib%diag%omepot(6)  = 0.0
        sib%diag%ansqr(6)   = 0.0
        sib%diag%wsfws(6)   = 0.0
        sib%diag%wsflt(6)   = 0.0
        sib%diag%wsfht(6)   = 0.0
        sib%diag%wci(6)     = 0.0
        sib%diag%whs(6)     = 0.0
        sib%diag%wags(6)    = 0.0
        sib%diag%wegs(6)    = 0.0
        sib%diag%assimnp(6) = 0.0
        sib%diag%assimci(6) = 0.0
        sib%diag%assimpot(6)= 0.0
        sib%diag%respc(6)   = 0.0
        
        rstfac3(6)          = 0.0
    

    phys_loop : do i=1,5

        if ( sib%param%physfrac(i) == 0.0 ) cycle phys_loop


        if( sib%param%phystype(i) == 3) then
            c3 = 1.
            c4 = 0.
        elseif(sib%param%phystype(i) == 4) then
            c3 = 0.
            c4 = 1.
        else
            print*,'loop index=',i,' phystype=',sib%param%phystype(i)
            stop'ERROR:UNKNOWN PHYSIOLOGY TYPE IN PHOSIB'

        endif


        !-----------------------------------------------------------------------
        !
        !     Q-10 AND STRESS TEMPERATURE EFFECTS
        !
        !      QT          (QT)    : TABLE (2)     , SE-92A
        !-----------------------------------------------------------------------

        qt = 0.1 * ( sib%prog%tc - sib%param%trop(i) )

        respn = sib%param%respcp(i) * sib%param%vmax0(i) * sib%diag%rstfac(2)


        !itb...patch to prevent underflow if temp is too cool...
        if(sib%prog%tc >= sib%param%trdm(i) )then
            sib%diag%respc(i) = respn * 2.0**qt                              &
                /( 1. + EXP( sib%param%trda(i)*(sib%prog%tc-sib%param%trdm(i) )))
        else
            sib%diag%respc(i) = respn * 2.0**qt
        endif

        vm = sib%param%vmax0(i) * 2.1**qt

        templ = 1. + EXP(sib%param%slti(i) * (sib%param%hltii(i) - sib%prog%tc))

        temph = 1. + EXP(sib%param%shti(i) * (sib%prog%tc - sib%param%hhti(i) ))

        rstfac3(i) = 1./( templ*temph)
        
        vm    = vm/temph * sib%diag%rstfac(2)*c3 &
            + vm * sib%diag%rstfac(2)*rstfac3(i) * c4

!print *,sib%param%phystype(i),sib%param%physfrac(i),vm,sib%param%vmax0(i)
        !-----------------------------------------------------------------------
        !
        !     MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY,
        !     COMPENSATION POINT       
        !
        !      ZKC          (KC)     : TABLE (2)     , SE-92A
        !      ZKO          (KO)     : TABLE (2)     , SE-92A
        !      SPFY         (S)      : TABLE (2)     , SE-92A
        !      GAMMAS       (GAMMA-*): TABLE (2)     , SE-92A
        !      OMSS         (OMEGA-S): EQUATION (13) , SE-92A
        !      BINTC        (B*ZLT)  : EQUATION (35) , SE-92A
        !-----------------------------------------------------------------------

        zkc     = 30. * 2.1**qt
        zko     = 30000. * 1.2**qt
        spfy    = 2600. * 0.57**qt
        gammas  = 0.5 * po2m/spfy * c3
        !itb...check for underflow...


        if(sib%prog%radvbc < 1.E-6) then
            sib%diag%pfd=0.0
        else
            sib%diag%pfd = 4.6E-6 * sib%param%gmudmu * &
                ( sib%prog%radvbc + sib%prog%radvdc )
        endif


        !...convert resistance to conductance  to  mol/ (m2 sec)
        !...(44.6 mol m^-3 conversion factor)
        if ( sib%prog%rst(i) == 0. ) sib%prog%rst(i) = sib%prog%rst(1)
        gsh2o  = 1.0/sib%prog%rst(i) * 44.032476*tprcor/sib%prog%tc
        gbh2o  = 0.5/sib%diag%rb     * 44.032476*tprcor/sib%prog%tc
        gah2o  = 1.0/sib%diag%ra     * 44.032476*tprcor/sib%prog%tm

        xgco2m = 4000.0 * sib%param%vmax0(i) * sib%diag%aparkk * sib%diag%rstfac(2)

        rrkk   = zkc*( 1. + po2m/zko ) * c3 &
            + sib%param%vmax0(i)/5.* ( 1.8**qt) * c4

        par    = sib%diag%pfd * sib%param%effcon(i) * ( 1.-scatg )

        soilfrztg = 1.+exp(-1.5 * &
            (max(270.0_dbl_kind,sib%prog%td(1))-273.16))
        soilfrztd = 1.+exp(-1.5 * &
            (max(270.0_dbl_kind,sib%prog%td(2))-273.16))
        soilfrz   = max(1./soilfrztg, 1./soilfrztd)
        soilfrz   = max( soilfrz, 0.05_dbl_kind)

        bintc  = sib%param%binter(i) * sib%param%zlt * sib%param%green * &
            sib%diag%rstfac(2) * soilfrz

        omss   = ( sib%param%vmax0(i) / 2.0 ) * ( 1.8**qt ) &
            /templ * sib%diag%rstfac(2) * c3 &
            + rrkk * sib%diag%rstfac(2) * c4

        !-----------------------------------------------------------------------
        !
        !     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
        !     ASSIMILATION RATE.
        !
        !-----------------------------------------------------------------------


        range1  = sib%prog%pco2m * ( 1. - 1.6/sib%param%gradm(i) ) - gammas
        icconv = 1

        do ic = 1, 6
            pco2y(ic) = 0.
            eyy(ic)   = 0.
        enddo


        !Bio...HERE IS THE ITERATION LOOP.
        !Bio...
        !Bio...We iterate on PCO2C-sortin makes a 'first guess' at
        !Bio...then orders PCO2C/Error pairs on increasing error size,
        !Bio...then uses a combination of linear/quadratic fit to obtain 
        !Bio...the 'next best guess' as iteration count increases.
        !Bio...CYCALC uses that value of PCO2C to get the next value 
        !Bio...of ASSIMN. CO2A and CO2S follow.

        do ic = 1, 6

            call sortin( eyy, pco2y, range1, gammas, ic)

            call cycalc( sib%diag%aparkk, VM, sib%param%atheta(i), &
                sib%param%btheta(i),par,gammas, sib%diag%respc(i),    &
                rrkk, omss, c3, c4,pco2y(ic), assimny(ic),      &
                assimy(ic))



            !pl now prognose the new CAS CO2 according to flux divergence
            !pl we are going to do this in mol C / mol air (same as PaC/PaAir)

            co2a    = sib%prog%pco2ap /   (sib%prog%ps*100.)
            co2m    = sib%prog%pco2m  /   (sib%prog%ps*100.)   

            co2a   = (  co2a + (dtt/co2cap) *  & 
                (sib%diag%respg - assimny(ic)  &
                +co2m*gah2o        ) )         &
                / (1+dtt*gah2o/ co2cap ) 


            pco2a = co2a * sib%prog%ps * 100.

            !itb...intermediate leaf surface CO2
            pco2s = pco2a - (1.4/gbh2o * assimny(ic) &
                * sib%prog%ps*100.)

            !itb...intermediate leaf internal CO2
            pco2i = pco2s - assimny(ic) * sib%prog%ps * 100.0 * 1.6/gsh2o

            !itb...intermediate leaf chloroplast CO2-this is what we iterate on 
            pco2c = pco2i - assimny(ic) * sib%prog%ps * 100.0 * 1.0/xgco2m

            eyy(ic) = pco2y(ic) - pco2c

            if(ic.ge.2) then
                ic1 = ic-1
                if(abs(eyy(ic1)).ge.0.1)then
                    icconv = ic
                else
                    eyy(ic) = eyy(ic1)
                    pco2y(ic) = pco2y(ic1)
                endif
            endif

        enddo !iteration loop


!print'(3g16.8)',sib%param%atheta(i),sib%diag%pfd,sib%diag%assimn(i)


        !itb...have iterated for physiology-specific values, now save them

        sib%diag%pco2c(i)  = pco2y(icconv)
        sib%diag%assimn(i) = assimny(icconv)
        sib%diag%assim(i)  = assimy(icconv)
        sib%diag%pco2i(i)  = sib%diag%pco2c(i) +  &
            sib%diag%assimn(i)/xgco2m*sib%prog%ps*100.0
        sib%diag%pco2s(i)  = sib%diag%pco2i(i) +  &
            sib%diag%assimn(i)/gsh2o *sib%prog%ps*100.0
        !
        !  update stomatal resistance...
        !
        h2oi   = sib_loc%etc / sib%prog%ps
        h2oa   =  sib%prog%ea / sib%prog%ps

        ecmole = 55.56 * sib%diag%ecmass * dti 

        h2os = h2oa + ecmole / gbh2o
        h2os  = min( h2os, h2oi )
        h2os  = max( h2os, 1.0e-7_dbl_kind)
        h2osrh = h2os / h2oi

        sib%diag%rstfac(1) = h2os/h2oi

        !Bio relaxed this condition to 1/10 of previous (.05 vs .5). The old way made
        !Bio the CO2 on top of the leaves always at least 1/2 of the value at the
        !Bio reference level.
        co2s = MAX(sib%diag%pco2s(i),sib%prog%pco2m*0.05) / (sib%prog%ps*100.)

        !Bio Ball-Berry equation right here !     
        gsh2oinf = (sib%param%gradm(i) *  &
            MAX(1.0e-12_dbl_kind,sib%diag%assimn(i))  &
            * h2osrh * soilfrz / co2s) + bintc

        !Bio this is the change in stomatal resistance
        !itb...this has been brought here from ADDINC 

        drst(i) = sib%prog%rst(i) * qdamp * ((gsh2o-gsh2oinf)/  &
            (pdamp*gsh2o+qdamp*gsh2oinf))

        bintc = bintc * sib%prog%tc / ( 44.032476 * tprcor)

        sib%prog%rst(i) = sib%prog%rst(i) + drst(i)

        ! bintc(i)- smallest canopy stomatal conductance needs to be passed in here.
        ! ---- c.zhang, 2/3/93

        sib%prog%rst(i)=MIN( 1./bintc, sib%prog%rst(i) )

        !...leaf conductance...
        sib%diag%ggl(i) = 1.0 / (sib%prog%rst(i)* sib%diag%rc)

        !-----------------------------------------------------------------------
        ! CALCULATION OF POTENTIAL ASSIMILATION
        !-----------------------------------------------------------------------

        ! Make assimn a top leaf, not the canopy.
        sib%diag%assimnp(i) = sib%diag%assimn(i) / sib%diag%aparkk

        ! Bottom stopped assim.
        sib%diag%antemp(i) = MAX(0.0_dbl_kind,sib%diag%assimnp(i))

        ! Potential intercellular co2.
        pco2ipot = sib%prog%ps*100.* (co2s-(1.6 * sib%diag%assimnp(i)/   &
            ((sib%param%gradm(i) * sib%diag%antemp(i) / co2s) + bintc)))

        ! Potential rubisco limitation.
        omcpot = sib%param%vmax0(i)*2.1**qt*((pco2ipot-gammas)/             &
            (pco2ipot+rrkk)*c3 + c4)

        ! Potential light limitation.
        sib%diag%omepot(i) = par*((pco2ipot-gammas)/                     &
            (pco2ipot+2.*gammas)*c3 + c4)

        ! Quad 1.
        sqrtin = MAX(0.0_dbl_kind,((sib%diag%omepot(i) + omcpot)**2-     &
            4.* sib%param%atheta(i) * sib%diag%omepot(i) * omcpot))

        ! Quad 1. Intermediate  top leaf photosynthesis.
        omppot = ((sib%diag%omepot(i)+omcpot)-SQRT(sqrtin))/  &
            (2.*sib%param%atheta(i) )

        ! Potential sink or pep limitation.
        omspot = (sib%param%vmax0(i) / 2.0)*(1.8**qt)*c3  &
            + rrkk*pco2ipot*c4

        ! Quad 2.
        sqrtin=MAX(0.0_dbl_kind,((omppot+omspot)**2-4.*sib%param%btheta(i)*   &
            omppot*omspot))

        ! Quad 2. Final Potential top leaf photosynthesis.
        if ( omppot < 1.0e-14 ) omppot = 0.0_dbl_kind

        sib%diag%assimpot(i) = ((omspot + omppot)-SQRT(sqrtin))/           &
            (2. * sib%param%btheta(i))
        !-----------------------------------------------------------------------
        ! CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
        !-----------------------------------------------------------------------

        ! Stressed rubisco limitation.
        omcci = vm*((pco2ipot-gammas)/(pco2ipot+rrkk)*c3+ c4)

        ! Quad 1.
        sqrtin = MAX(0.0_dbl_kind,(sib%diag%omepot(i) + omcci) **2 -       &
            4.*sib%param%atheta(i) * sib%diag%omepot(i) * omcci)

        ! Quad 1. Intermediate stress limited top leaf photosynthesis.
        ompci = ((sib%diag%omepot(i) + omcci) - SQRT(sqrtin))              &
            /(2.*sib%param%atheta(i))

        ! Stressed sink or pep limitation.
        omsci = omss*(c3 + pco2ipot*c4)

        ! Quad 2.
        sqrtin = MAX(0.0_dbl_kind,(ompci + omsci)**2-4. * sib%param%btheta(i)  & 
            * ompci * omsci)

        ! Quad 2. Final stress limited top leaf photosynthesis.
        sib%diag%assimci(i) = ((omsci+ompci)-SQRT(sqrtin))/(2.*sib%param%btheta(i))

        !-----------------------------------------------------------------------
        ! CALCULATION OF CONTROL COEFFICIENTS
        !-----------------------------------------------------------------------

        ! Intermediate.
        dompdomc = (ompci-sib%diag%omepot(i) )/                            &
            (2.*sib%param%atheta(i)*ompci-omcci-sib%diag%omepot(i))

        ! Bottom stopped final stress limited top leaf photosynthesis.
        ascitemp = MAX(sib%diag%assimci(i),1.0e-12_dbl_kind)

        ! Rubisco control coefficient.
        ccomc = (dompdomc*(sib%diag%assimci(i)-omsci)/  &
            (2.*sib%param%btheta(i)*sib%diag%assimci(i)-ompci-omsci))*  &
            omcci/ascitemp

        ! Sink or pep control coefficient.
        ccoms = ((sib%diag%assimci(i)-ompci)/  &
            (2.*sib%param%btheta(i)*sib%diag%assimci(i)-ompci-omsci))*  &
            omsci/ascitemp

        !-----------------------------------------------------------------------
        !  OUTPUT:  POTENTIAL ASSIMILATION RATES TO BE SUMMED
        !-----------------------------------------------------------------------
        ! Canopy values (overwrites top leaf).

        sib%diag%omepot(i)   = sib%diag%omepot(i)   * sib%diag%aparkk
        sib%diag%assimpot(i) = sib%diag%assimpot(i) * sib%diag%aparkk
        sib%diag%assimci(i)  = sib%diag%assimci(i)  * sib%diag%aparkk
        sib%diag%assim(i)    = sib%diag%assim(i)    * sib%diag%aparkk
        sib%diag%antemp(i)   = sib%diag%antemp(i)   * sib%diag%aparkk
        sib%diag%ansqr(i)    = sib%diag%antemp(i)   * sib%diag%antemp(i)
        sib%diag%assimnp(i)  = sib%diag%assimnp(i)  * sib%diag%aparkk

        !-----------------------------------------------------------------------
        ! OUTPUT:  WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE SUMMED
        !-----------------------------------------------------------------------

        ! Water stress.
        sib%diag%wsfws(i) = sib%diag%assimpot(i)*(1.-sib%diag%rstfac(2))*  &
            (ccomc+ccoms)

        ! High temperature stress.
        sib%diag%wsfht(i) = sib%diag%assimpot(i)*(1.-1./temph)*ccomc

        ! Low temperature stress.
        sib%diag%wsflt(i) = sib%diag%assimpot(i)*(1.-1./templ)*  &
            (ccoms*c3+ccomc*c4)

        !  protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
        cwsfws = (1.-sib%diag%rstfac(2))*(ccomc+ccoms)
        if(cwsfws .gt.1. .or. cwsfws .lt. 0.) sib%diag%wsfws(i)=0.

        cwsfht = (1.-1./temph)*ccomc
        if(cwsfht > 1.0_dbl_kind .or. cwsfht < 0.0_dbl_kind)  &
            sib%diag%wsfht(i)=0.

        cwsflt = (1.-1./templ)*(ccoms*c3+ccomc*c4)
        if(cwsflt > 1.0_dbl_kind .or. cwsflt < 0.0_dbl_kind)  &
            sib%diag%wsflt(i)=0.


        ! Intermediate assimilation weighted Ci.
        sib%diag%wci(i) = sib%diag%antemp(i) * sib%diag%pco2i(i)

        ! Intermediate assimilation weighted relative humidty stress factor.
        sib%diag%whs(i) = sib%diag%antemp(i) * sib%diag%rstfac(1)

        ! Intermediate assimilation weighted stomatal conductance.
        sib%diag%wags(i) = gsh2o * sib%diag%antemp(i)

        ! Intermediate evaporation weighted stomatal conductance.(Step 1.
        !   Step 2 after subroutine update)
        sib%diag%wegs(i) = gsh2o

        !
        !itb...determine the weighted mean value for the gridcell
        !
        sib%diag%pco2c(6)   = sib%diag%pco2c(6)  + sib%diag%pco2c(i)  *  &
            sib%param%physfrac(i)
        sib%diag%pco2i(6)   = sib%diag%pco2i(6)  + sib%diag%pco2i(i)  *  &
            sib%param%physfrac(i)
        sib%diag%pco2s(6)   = sib%diag%pco2s(6)  + sib%diag%pco2s(i)  *  &
            sib%param%physfrac(i)
        sib%diag%assimn(6)  = sib%diag%assimn(6) + sib%diag%assimn(i) *  &
            sib%param%physfrac(i)
        sib%diag%assim(6)   = sib%diag%assim(6)  + sib%diag%assim(i)  *  &
            sib%param%physfrac(i)
        rstfac3(6)     = rstfac3(6)    + rstfac3(i)    * sib%param%physfrac(i) 
        sib%prog%rst(6)     = sib%prog%rst(6)    + sib%prog%rst(i)    *  &
            sib%param%physfrac(i)
        sib%diag%ggl(6)     = sib%diag%ggl(6)    + sib%diag%ggl(i)    *  &
            sib%param%physfrac(i)
        sib%diag%antemp(6)  = sib%diag%antemp(6) + sib%diag%antemp(i) *  &
            sib%param%physfrac(i)
        sib%diag%omepot(6)  = sib%diag%omepot(6) + sib%diag%omepot(i) *  &
            sib%param%physfrac(i)
        sib%diag%ansqr(6)   = sib%diag%ansqr(6)  + sib%diag%ansqr(i)  *  &
            sib%param%physfrac(i)
        sib%diag%wsfws(6)   = sib%diag%wsfws(6)  + sib%diag%wsfws(i)  *  &
            sib%param%physfrac(i)
        sib%diag%wsflt(6)   = sib%diag%wsflt(6)  + sib%diag%wsflt(i)  *  &
            sib%param%physfrac(i)
        sib%diag%wsfht(6)   = sib%diag%wsfht(6)  + sib%diag%wsfht(i)  *  &
            sib%param%physfrac(i)
        sib%diag%wci(6)     = sib%diag%wci(6)    + sib%diag%wci(i)    *  &
            sib%param%physfrac(i)
        sib%diag%whs(6)     = sib%diag%whs(6)    + sib%diag%whs(i)    *  &
            sib%param%physfrac(i)
        sib%diag%wags(6)    = sib%diag%wags(6)   + sib%diag%wags(i)   *  &
            sib%param%physfrac(i)
        sib%diag%wegs(6)    = sib%diag%wegs(6)   + sib%diag%wegs(i)   *  &
            sib%param%physfrac(i)
        sib%diag%assimnp(6) = sib%diag%assimnp(6)+ sib%diag%assimnp(i)*  &
            sib%param%physfrac(i)
        sib%diag%assimci(6) = sib%diag%assimci(6)+ sib%diag%assimci(i)*  &
            sib%param%physfrac(i)
        sib%diag%assimpot(6)= sib%diag%assimpot(6)+sib%diag%assimpot(i)* &
            sib%param%physfrac(i)
        sib%diag%respc(6)   = sib%diag%respc(6)  + sib%diag%respc(i)  *  &
            sib%param%physfrac(i)


        !...CFRAX...
        !...at the end of each physiology loop, call the CFRAX code

        call cfrax_physloop(sib,i,c3)

        !...CFRAX...

    enddo phys_loop   ! PHYSIOLOGY LOOP

    !pl now do the real C_A forecast with the iterated fluxes.

    co2a    = sib%prog%pco2ap /   (sib%prog%ps*100.)
    co2m    = sib%prog%pco2m  /   (sib%prog%ps*100.) 

    !itb...carbon flux between CAS and reference level (mol C m^-2 sec^-1)
    sib%diag%cflux = gah2o*(co2a-co2m)

    sib%prog%expand=sib%prog%pco2ap*sib%diag%cas_cap_co2/rstar/sib%prog%ta- sib%prog%cas_old
  
! original semi-implicit time differencing
    co2a = (co2a + (dtt/co2cap) * (sib%diag%respg - sib%diag%assimn(6)    & 
        +co2m*gah2o ) ) / (1+dtt*gah2o/co2cap)
    sib%prog%pco2ap = co2a * sib%prog%ps * 100.
!
! basic forward time differencing
!    sib%prog%pco2ap=sib%prog%pco2ap+(sib%diag%respg - sib%diag%assimn(6)  &
!    -sib%diag%cflux)*dtt*rstar*sib%prog%ta/sib%diag%cas_cap_co2

! semi-implicit time differencing with cflux damping
!    if (sib%diag%cflux*dtt > 0.01*sib%prog%cas_old) then 

!      co2a_star = (co2a + (dtt/co2cap) * (sib%diag%respg - sib%diag%assimn(6)    & 
!        +co2m*gah2o ) ) / (1+dtt*gah2o/co2cap)

!      co2a = 0.4*co2a_star + 0.6* co2a
!      sib%prog%pco2ap = co2a * sib%prog%ps * 100.
!      print*, sib%prog%pco2ap,sib%prog%pco2m
!    else
!      co2a = (co2a + (dtt/co2cap) * (sib%diag%respg - sib%diag%assimn(6)    & 
!        +co2m*gah2o ) ) / (1+dtt*gah2o/co2cap)
!      sib%prog%pco2ap = co2a * sib%prog%ps * 100.
!    endif
    sib%diag%cflux = gah2o*(co2a-co2m)

!
! moles per m2 CO2 in canopy air space
    sib%prog%cas=sib%prog%pco2ap*sib%diag%cas_cap_co2/rstar/sib%prog%ta

    sib%diag%rstfac(3) = rstfac3(6)
    sib%diag%rstfac(4) = sib%diag%rstfac(1) * sib%diag%rstfac(2) *  &
        sib%diag%rstfac(3)

    !...CFRAX...
    !...one last call to calculate canopy-mean discrimination values

    call cfrax_final(sib)

    !...CFRAX...


end subroutine phosib

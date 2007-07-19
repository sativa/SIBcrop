
!==================SUBROUTINE DELEF======================================

subroutine delef(sib,sib_loc)    

    use kinds
    use sibtype

    use physical_parameters, only: &
        cp => spec_heat_cp, &
        hltm
    use sib_const_module, only: &
        snofac, &
        dtt

    !========================================================================
    !
    !     Calculation of partial derivatives of canopy and ground latent
    !        heat fluxes with respect to Tc, Tg, Theta-m, and Qm.
    !     Calculation of initial latent heat fluxes.
    !
    !pl the ETC, ETG and so on are the vapor pressure at temps TC, TG and so on
    !pl the BTC, BTG are the derivatives of ETC, ETG with relation to TC, TG etc.
    !
    !======================================================================== 

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       EC             ECT + ECI
    !       EG             EGS + EGI
    !       ECDTC          dEC/dTC
    !       ECDTG          dEC/dTG
    !       ECDQM          dEC/dQM
    !       EGDTC          dEG/dTC
    !       EGDTG          dEG/dTG
    !       EGDQM          dEG/dQM
    !       BBC            dE/dTC
    !       BBG            dE/dTG
    !       BBM            dE/dQM
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    implicit none

    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    type(sib_local_vars)     ,intent(inout) :: sib_loc
    ! variables local to SiB

    !----------------------------------------------------------------------  


    ! local variables
    real(kind=dbl_kind) :: cpdpsy    ! cp/psy

    !     MODIFICATION FOR SOIL DRYNESS : HR=REL. HUMIDITY IN TOP LAYER     

    cpdpsy = cp / sib%diag%psy

    !-----------------------------------------------------------------------
    !                                                                       
    !     CALCULATION OF SURFACE RESISTANCE COMPONENTS, SEE EQUATIONS (64,66)
    !       OF SE-86                                   
    !pl     gect(i)  =      (1. -wc(i)) /  rc(i)
    !pl     geci(i)  = epsc(i) * wc(i)  / (RB(I) + RB(I))
    !pl     gegs(i)  =        (1-wg(i)) / (rds(i))
    !pl     gegi(i)  = epsg(i) * wg(i)  /  rd(i)
    !                                                                       
    !-----------------------------------------------------------------------

    sib_loc%cog1 = (sib_loc%gegi + sib_loc%gegs*sib%diag%hrr)
    sib_loc%cog2 = (sib_loc%gegi + sib_loc%gegs        )


    !            D2(I)   = 1.0 / RA(I) + COC(I) + COG2(I)
    !-----------------------------------------------------------------------
    !                                                                       
    !     FLUXES EXPRESSED IN JOULES M-2   CPL WHY ?????
    !                                                                       
    !      ec         (EC)    : EQUATION (64) , SE-86
    !      eg         (EG)    : EQUATION (66) , SE-86
    !      es         (ES)    : EQUATION (66) , SE-86
    !      ea         (EA)    : EQUATION ????
    !-----------------------------------------------------------------------

    !pl these are the current time step fluxes in J/m2  

    !pl notice that the fluxes are already limited by the altered e*(T) values

print*,'hi there'

    sib%diag%ec  = (sib_loc%etc - sib%prog%ea) * sib_loc%coc * &
        sib%prog%ros  * dtt * cpdpsy

    sib%diag%eg  = ( sib_loc%etg * sib_loc%cog1 &
        - sib%prog%ea * sib_loc%cog2) * &
        sib%prog%ros * dtt * cpdpsy

    sib%diag%es  = ((sib_loc%ets - sib%prog%ea)/sib%diag%rd )* &
        sib%prog%ros * dtt * cpdpsy/snofac

    sib%diag%fws = ((sib%prog%ea  - sib%prog%em ) / sib%diag%ra) &
        * sib%prog%ros * dtt * cpdpsy

    !pl now we do the partial derivatives  these assume W/m2

    !pl for the canopy leaves vapor pressure: W/ (m2* K)
    sib_loc%ecdtc =    sib_loc%getc * sib_loc%coc * sib%prog%ros * cpdpsy

    sib_loc%ecdea = - sib_loc%coc * sib%prog%ros * cpdpsy

    !pl for ground latent heat fluxes: W/ (m2* K)
    sib_loc%egdtg =   sib_loc%getg * sib_loc%cog1 * sib%prog%ros * cpdpsy

    sib_loc%egdea = - sib_loc%cog2 * sib%prog%ros * cpdpsy               

    !pl for snow latent heat fluxes: W/ (m2* K)
    sib_loc%esdts =   sib_loc%gets * sib%prog%ros * cpdpsy / sib%diag%rd

    !pl for snow latent heat fluxes: W/ (m2 * Pa)
    sib_loc%esdea = - sib%prog%ros * cpdpsy / sib%diag%rd              

    !pl for CAS latent heat fluxes: W/ (m2* Pa)
    sib_loc%eadea = sib%prog%ros * cpdpsy / sib%diag%ra 

    sib_loc%eadem = - sib_loc%eadea            

    !PL ATTENTION !!!! DANGER !!! do not use without sibdrv = true
    !pl these all need to be re-done for the GCM (no sibdrv)
    !-----------------------------------------------------------------------
    !      BBC       (dE/dTC)  : EQUATION (13) , SA-89B
    !      BBG       (dE/dTG) : EQUATION (13) , SA-89B
    !      BBM       (dE/dQM)  : EQUATION (13) , SA-89B
    !-----------------------------------------------------------------------
    !        BBG(I) = (COG1(I) / D2(i))
    !     *          * btg(I) * 0.622 * ps(i)           
    !     *       / ((ps(i) - etg(I)) * (ps(i) - etg(I)))                
    !        BBC(I) = (COC(I)  / D2(i))
    !     *            * getc(I) * 0.622 * ps(i)           
    !     *       / ((ps(i) - etc(I)) * (ps(i) - etc(I)))                
    !        BBM(I) = 1.0   / (ra(I)  * D2(i))                             


end subroutine delef                                              

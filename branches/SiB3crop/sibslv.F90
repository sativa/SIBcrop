
!======================SUBROUTINE SIBSLV=================================

subroutine sibslv(sib,sib_loc)   

!========================================================================
!
!     Calculation of time increments in Tc, Tgs, Theta-m and Qm using an
!        implicit backwards method with explicit coefficients.  
!pl   Similar to equations (10-15), SA-92B. 
!
!     Longwave feedbacks are now really included
!
!======================================================================== 

!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       DTC            CANOPY TEMPERATURE INCREMENT (K)
!       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
!       DTH            MIXED LAYER POTENTIAL TEMPERATURE INCREMENT (K)
!       DQM            MIXED LAYER MIXING RATIO INCREMENT (KG KG-1)
!       ETMASS (FWS)   EVAPOTRANSPIRATION (MM)
!       HFLUX (FSS)    SENSIBLE HEAT FLUX (W M-2)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use kinds
use sibtype

!itb...used for diagnostic print...
use sib_const_module, only:  &
    nsoil,  &
    dtt,    &
    dti

use physical_parameters, only:  &
    grav,                &
    cp => spec_heat_cp,  &
    hltm
IMPLICIT none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

type(sib_local_vars)     ,intent(inout) :: sib_loc
                                   ! variables local to SiB

!----------------------------------------------------------------------  



real(kind=dbl_kind),parameter :: grav2 = grav * 0.01_dbl_kind

integer(kind=int_kind) :: error_flag,i,j,k

real(kind=dbl_kind):: aag
real(kind=dbl_kind):: aac
real(kind=dbl_kind):: aam
real(kind=dbl_kind):: bbg
real(kind=dbl_kind):: bbc
real(kind=dbl_kind):: bbm
real(kind=dbl_kind):: temv 

!...matrix arrays - variable, due to potential snow layers
real(kind=dbl_kind),dimension(15-sib%prog%nsl,15-sib%prog%nsl) :: avec
real(kind=dbl_kind),dimension(15-sib%prog%nsl)   :: bvec,bv_copy
integer(kind=int_kind),dimension(15-sib%prog%nsl)   :: cvec   

!...this routine sets up the coupled system of partial 
!...differential equations described in Sato et al., 
!...with the exception that now Ta and ea are prognostic
!...variables, and so we have two more equations, reflected
!...in two more lines and two more columns as compared 
!...to the old sibslv.F (used for no prognistic CAS calculations)

!...this matrix is expandable, size determined by the number of 
!...snow layers (0 to 5). All energy is passed through the snow layer
!...if snow exists. Partial snowcover is not dealt with. We know
!...that this is physically incorrect, it is an issue that we hope
!...to deal with in the future...


!itb...VARIABLES
!...1 : TREF (lowest model layer/driver data temperature)
!...2 : EREF (lowest model layer/driver data water vapor)
!...3 : TA   (CAS temperature) 
!...4 : EA   (CAS water vapor)
!...5 : TC   (vegetation temperature)
!...6 : TG   (Ground or snow surface temperature)
!...7-14 to 19 : TD (interior soil layers temperature)
!...last member: TD (bottom soil level)


!*****************************************************
!
!   THESE FIRST TWO EQUATIONS (TREF AND EREF) ARE
!   FOR THE 'OFFLINE SiB' SITUATION, WHERE SiB IS
!   DRIVEN BY TOWER DATA. FOR COUPLING TO A MODEL
!   (GCM OR MESOSCALE MODEL), USE THE SECOND SET
!   OF EQUATIONS. REFERENCE IS TO KALNAY AND KANAMITSU
!   (1988) RIGHT NOW, BUT I'LL WRITE IT UP LATER...
!
!*****************************************************

!itb...set up a_areas
      sib%diag%a_areas = sib%diag%areas

!     if(sib%diag%areas > 0.0) then
!       sib%diag%a_areas = 1.0
!     else
!       sib%diag%a_areas = 0.0
!     endif

    !1     TREF EQUATION       

    avec(1,1)             =  1.0                                 
    avec(1,2:15-sib%prog%nsl)  =  0.0            
    bvec(1)               =  0.0  

    !2     EREF EQUATION  

    avec(2,1)             =  0.0                                 
    avec(2,2)             =  1.0                                 
    avec(2,3:15-sib%prog%nsl)  =  0.0            
    bvec(2)               =  0.0        

!********************************************************************
!     TREF EQUATION - COUPLED
!        avec(1,1)             =  (sib%cp * sib%prog%psb)/(grav*    &
!                                  (0.5_dbl_kind * dtt)) + sib_loc%hadta
!        avec(1,2)             = 0.0
!        avec(1,3)             = -sib_loc%hadta
!        avec(1,4:15-sib%prog%nsl)  = 0.0
!        bvec(2)               = sib%diag%fss * dti
!
!     EREF EQUATION - COUPLED
!        avec(2,1)             = 0.0
!        avec(2,2)             = (sib%cp * sib%prog%psb)/    &
!                                (grav * (0.5_dbl_kind * dtt)    &
!                                 * sib%diag%psy) + sib_loc%eadea
!        avec(2,3)             = sib_loc%eadem
!        avec(2,4:15-sib%prog%nsl)  = 0.0
!        bvec(3)               = sib%diag%fws * dti
!*********************************************************************        

    !3     TA EQUATION 
    !itb...start with the 2delta_t stuff here...

!itb...zero all of 'em out first, then fill 'em in...
    avec(3,:) = 0.0

    avec(3,1)  =  sib_loc%hadth                                    
    avec(3,2)  =  0.0                                 
    avec(3,3)  =  sib%diag%cas_cap_heat   * (0.5_dbl_kind *dti )  & 
        + sib_loc%hadta  - sib_loc%hcdta  &
        - (1.-sib%diag%areas)*sib_loc%hgdta - sib%diag%areas*sib_loc%hsdta
    avec(3,4)  =  0.0         
    avec(3,5)  =  - sib_loc%hcdtc     


!itb...my latest crack at this...(15 august)
    if(sib%diag%areas > 0.0 )then           
      avec(3,6)  =   - sib_loc%hsdts * sib%diag%areas 
      avec(3,6-sib%prog%nsl) = -sib_loc%hgdtg * (1.0 - sib%diag%areas)
    else
      avec(3,6) = - sib_loc%hgdtg
    endif

    avec(3,7:15-sib%prog%nsl)  =  0.0   

    bvec(3)    =  sib%diag%hc * dti - sib%diag%fss * dti  &
        + (1.-sib%diag%a_areas)*sib%diag%hg * dti  &
        +     sib%diag%a_areas *sib%diag%hs * dti    

    !4    EA EQUATION  
    avec(4,:) = 0.0


    avec(4,1)  =  0.0                                  
    avec(4,2)  =  sib_loc%eadem                                  
    avec(4,3)  =  0.0                                  
    avec(4,4)  =  sib%diag%cas_cap_vap   * (0.5_dbl_kind * dti)  &
        + sib_loc%eadea  - sib_loc%ecdea        &
        -  (1.-sib%diag%areas)*sib_loc%egdea  &
        -      sib%diag%areas *sib_loc%esdea  
    avec(4,5)  =  - sib_loc%ecdtc               
    avec(4,6)  =  - sib_loc%egdtg * (1.-sib%diag%areas)  

    avec(4,7:15-sib%prog%nsl)  =  0.0             
    bvec(4)    =  (sib%diag%ec * dti  -  sib%diag%fws * dti  &
        +  (1.-sib%diag%a_areas)*sib%diag%eg * dti   &
        +      sib%diag%a_areas *sib%diag%es * dti)    

    !5    TC EQUATION 

    avec(5,1)  =  0.0                                   
    avec(5,2)  =  0.0                               
    avec(5,3)  =  sib_loc%hcdta                                    
    avec(5,4)  =  sib_loc%ecdea          
    avec(5,5)  =  sib%param%czc * (0.5_dbl_kind * dti) + sib_loc%hcdtc    &
        + sib_loc%ecdtc + sib_loc%lcdtc
    avec(5,6)  =  sib_loc%lcdtg

    avec(5,7:15-sib%prog%nsl)  =  0.0             

    bvec(5)    =  sib%diag%radt(1) - sib%diag%hc * dti - sib%diag%ec * dti      

    !6    TOP SOIL LAYER (TD1 OR TG) EQUATION

    !itb...will need an 'if' statement here for snow/no snow case

    if(sib%prog%nsl == 0) then !NO SNOW CASE

        avec(6,1)  =  0.0                                  
        avec(6,2)  =  0.0
        avec(6,3)  =  sib_loc%hgdta                                  
        avec(6,4)  =  sib_loc%egdea          
        avec(6,5)  =  sib_loc%lgdtc               
        avec(6,6)  =  sib%param%slamda(1) +  sib%param%shcap(1) *     &
            (0.5_dbl_kind * dti)          &
            + sib_loc%hgdtg + sib_loc%egdtg + sib_loc%lgdtg          
        avec(6,7)  =  sib%param%slamda(1)
 
        avec(6,8:15-sib%prog%nsl)  =  0.0                    
         
        bvec(6)    =  sib%diag%radt(2) - sib%diag%hg * dti - sib%diag%eg * dti &
            + sib%param%slamda(1) * (sib%prog%td(2) - sib%prog%td(1))

    else   ! SNOW CASE

        avec(6,1)  =  0.0                                  
        avec(6,2)  =  0.0
        avec(6,3)  =  sib_loc%hsdta*sib%diag%areas                                  
        avec(6,4)  =  sib_loc%esdea *sib%diag%areas         
        avec(6,5)  =  sib_loc%lsdtc               
        avec(6,6)  =  sib%param%slamda(sib%prog%nsl+1)*sib%diag%areas +  &
            sib%param%shcap(sib%prog%nsl+1) * (0.5_dbl_kind * dti) +   &
            sib_loc%hsdts + sib_loc%esdts + sib_loc%lsdts          
        avec(6,7)  =  sib%param%slamda(sib%prog%nsl+1)

        avec(6,8:15-sib%prog%nsl)  =  0.0              

        bvec(6)    =  sib%diag%radt(3) - sib%diag%hs * dti*sib%diag%areas -     &
            sib%diag%es * dti*sib%diag%areas +                           &
            sib%param%slamda(sib%prog%nsl+1) * (sib%prog%td(sib%prog%nsl+2) - &
            sib%prog%td(sib%prog%nsl+1))*sib%diag%areas
    endif

    !7-??    INTERIOR SOIL LAYERS

    !itb...need to scale this for snow/no snow (sib%prog%nsl)

    do i = 7, 14 - sib%prog%nsl   ! matrix indices
        k = i - 5 + sib%prog%nsl    ! soil layer indices

    !itb...fill all matrix components with zero, re-fill the 3 necessary
    !   ...positions 
        do j=1, 15  - sib%prog%nsl
            avec(i,j) = 0.0
        enddo

        avec(i,i-1) = -sib%param%slamda(k-1)

        avec(i,i)   = sib%param%shcap(k)*dti + sib%param%slamda(k)  &
            + sib%param%slamda(k-1)

        avec(i,i+1) = -sib%param%slamda(k)

        bvec(i)     = sib%param%slamda(k)*(sib%prog%td(k+1) - sib%prog%td(k))  &
            - sib%param%slamda(k-1) * (sib%prog%td(k) - sib%prog%td(k-1))


    enddo




    !    BOTTOM SOIL LAYER
    i = 15-sib%prog%nsl
    do j=1,i-2
        avec(i,j)  =  0.0 
    enddo                                
    avec(i,i-1) =  sib%param%slamda(9) * (sib%prog%td(10) - sib%prog%td(9))             
    avec(i,i) =  sib%param%shcap(10) * dti  &
        + sib%param%slamda(9) 
    bvec(i)    =  - sib%param%slamda(9) * (sib%prog%td(10) - sib%prog%td(9))


    !     SOLVE MATRIX EQUATION   
    !call dsimul(avec,bvec,cvec,1,15-sib%prog%nsl,15-sib%prog%nsl,error_flag)
    !jlc...use the lapack version of this call, instead of our home-brewed
    !FIXME: cvec is temporarily in the place of IPIV (the pivot array)

    call dgesv( 15-sib%prog%nsl, 1, avec, 15-sib%prog%nsl, cvec, bvec, &
        15-sib%prog%nsl, error_flag )
 
    sib_loc%dth = bvec(1)           
    sib_loc%dqm = bvec(2)
    sib_loc%dta = bvec(3)
    sib_loc%dea = bvec(4)
    sib_loc%dtc = bvec(5)        
    sib_loc%dtg = bvec(6)   
    sib_loc%dts = 0.0 

    sib_loc%dtd(sib%prog%nsl+1) = sib_loc%dtg
    do i=7,15-sib%prog%nsl
        sib_loc%dtd(i-5+sib%prog%nsl) = bvec(i)
 
    enddo


end subroutine sibslv

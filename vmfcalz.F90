subroutine vmfcalz(sib,zzwind,zztemp,cuni)

    use kinds
    use sibtype
    use sib_const_module, only:  &
        vkrmn
    use physical_parameters, only:  &
        grav,  &
        delta


    implicit none

    !----------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    !----------------------------------------------------------------  

!*****************************************************************************                                                                       
!     VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum 
!     and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
!     and Boville(1993), and by Beljaars and Holtslag (1991).              
!  
!     Rerences:
!       Beljars and Holtslag (1991): Flux parameterization over land surfaces
!              for atmospheric models. J. Appl. Meteo., 30, 327-341.
!       Holtslag and Boville (1993): Local versus nonlocal boundary-layer 
!              diffusion in a global climate model. J. of Climate, 6, 1825-
!              1842.
!       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
!              atmosphere. Boundary-Layer Meteo., 17, 187-202.
!       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
!              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
!              Layer parameterization, ECMWF, 59-79.
!
!     General formulation:
!        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.) 
!     Transfer coefficients for mommentum and heat fluxes are:
!        CU = CUN*Fm, and
!        CT = CTN*Fh
!        where  CUN and CTN are nutral values of momentum and heat transfers,
!           and Fm and Fh are stability functions derived from surface
!           similarity relationships.     
!*****************************************************************************



    real(kind=dbl_kind) ,intent(in) ::  &
        zzwind, & ! adjusted wind measurement height (m)       
        zztemp    ! adjusted temp measurement height (m)

    !...OUTPUT VARIABLES
    real(kind=dbl_kind) ,intent(out) :: &
        cuni    ! 1/ momentum transfer coefficient

    !...PATCH WARNING: an unjustified patch has been put in the code, 
    !...whereupon when cuni=1/cun is calculated, the square root is taken.
    !...this is a patch that makes the results better, but it is 
    !...unjustified scientifically.


    !...LOCAL VARIABLES
    integer ::   i

    !  constants for surface flux functions, according to Holtslag and
    !      Boville (1993, J. Climate)
    real(kind=dbl_kind), parameter ::  &
        bunstablM = 10.0, & !  
        bunstablT = 15.0, & !
        cunstablM = 75.0, & ! 
        cunstablT = 75.0, & !
        bstabl =  8.0,    & !
        cstabl = 10.0       ! 

    real(kind=dbl_kind) ::   &
        wgm,     & ! moisture mixing ratio deficit, 
                   !  CAS to reference layer (kg/kg)
        thgm,    & ! temperature difference (theta) CAS-ref level (K)
        z1z0u,   & ! ratio of reference height to roughness length
        z1z0urt, & ! square root of z1z0u
        z1z0t,   & ! ratio of reference height to roughness length
        z1z0trt, & ! square root of z1z0t
        !...currently, z1z0u and z1z0t are identical. theoretically, they 
        !...can be changed for different wind/temp measurement heights. 
        !...they both use zzwind right now.
        cun,     & ! momentum transfer coefficient (?)
        ctn,     & ! thermal transfer coefficient (?)
        temv,    & ! part of Richardson No. calculation
        zrib,    & ! part of Richardson No. calculation
        rib,     & ! Richardson Number
        fmomn,   & !
        fheat      !

    real(kind=dbl_kind) ::  &
        ribtemp, & !
        dm,      & !
        dh         !

    zrib = zzwind **2 / zztemp                                                    
    wgm  = sib%prog%sha - sib%prog%sh   
                                   
    !                                                                       
    ! SFC-AIR DEFICITS OF MOISTURE AND POTENTIAL TEMPERATURE         
    ! WGM IS THE EFFECTIVE SFC-AIR TOTAL MIXING RATIO DIFFERENCE.    
    !                                                                       
    thgm  = sib%prog%tha  - sib%prog%thm                                   
    sib%diag%thvgm = thgm + sib%prog%tha * delta * wgm       

    !   Ratio of reference height (zwind/ztemp) and roughness length:
    z1z0u = zzwind/sib%param%z0
    z1z0urt = sqrt( z1z0U )
    z1z0u = log( z1z0U )
    z1z0t = zzwind/sib%param%z0
    z1z0trt = sqrt( z1z0T )
    z1z0t = log( z1z0T )

    !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:

    cun = vkrmn*vkrmn / (z1z0u*z1z0u )   !neutral Cm & Ct
    ctn = vkrmn*vkrmn / (z1z0t*z1z0t )

    !...PATCH-when 1/cun is calculated, the square root is taken.
    cuni = z1z0u / vkrmn

    !                                                                       
    !   SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.            
    !   RIB IS THE BULK RICHARDSON NUMBER, between reference
    !   height and surface.

    temv = sib%prog%tha * sib%prog%spdm * sib%prog%spdm   
    temv = max(0.000001_dbl_kind,temv)
    rib = -sib%diag%thvgm * grav * zrib / temv 

    !   The stability functions for momentum and heat/moisture fluxes as
    !   derived from the surface-similarity theory by Luis (1079, 1982), and
    !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag 
    !   (1991).

    if(rib .ge. 0.0) then                                           

        !  THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT              

        rib   = min( rib, 0.5_dbl_kind)                   
        fmomn = (1. + cstabl * rib * (1.+ bstabl * rib))
        fmomn = 1. / fmomn
        fmomn = max(0.0001_dbl_kind,fmomn)
        fheat = fmomn

    else                                  

        !  THE UNSTABLE CASE.    

        ribtemp = abs(rib)
        ribtemp = sqrt( ribtemp )
        dm      = 1. + cunstablM * cun * z1z0Urt * ribtemp
        dh      = 1. + cunstablT * ctn * z1z0Trt * ribtemp
        fmomn   = 1. - (bunstablM * rib ) / dm
        fheat   = 1. - (bunstablT * rib ) / dh

    endif    

    !   surface-air transfer coefficients for momentum CU, for heat and 
    !   moisture CT. The CUI and CTI are inversion of CU and CT respectively.

    sib%diag%cu = cun * fmomn 
    sib%diag%ct = ctn * fheat

    !   Ustar and ventlation mass flux: note that the ustar and ventlation 
    !   are calculated differently from the Deardoff's methods due to their
    !   differences in define the CU and CT.

    sib%diag%ustar  = sib%prog%spdm * sib%prog%spdm * sib%diag%cu 
    sib%diag%ustar  = sqrt( sib%diag%ustar ) 
    sib%diag%ventmf = sib%prog%ros * sib%diag%ct * sib%prog%spdm  

    !                                                                       
    ! Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) 
    ! FREE CONVECTION VALUE, like DEARDORFF EQ(40B), because the 
    ! above CU and CT included free convection conditions.                                            

end subroutine vmfcalz

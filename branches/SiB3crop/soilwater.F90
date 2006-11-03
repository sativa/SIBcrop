
subroutine soilwater( sib, infil, hk, dhkdw, dw_liq, dzmm, zmm)

!----------------------------------------------------------------------
!
!     subroutine based on  CLM_SOILWATER
!
!     CLM web info: http://clm.gsfc.nasa.gov
!
!     Description: (taken directly from CLM_SOILWATER.F90)
!     Soil moisture is predicted from a 10-layer model (as with soil
!     temperature), in which the vertical soil moisture transport is 
!     governed by infiltration, runoff, gradient diffusion, gravity 
!     and root extraction through canopy transpiration.  The net water
!     applied to the surface layer is the snowmelt plus precipitation 
!     plus the throughfall of canopy dew minus surface runoff and 
!     evaporation. 
!
!     The vertical water flow in an unsaturated porous media is 
!     described by Darcy's Law, and the hydraulic conductivity and the
!     soil negative potential vary with soil water content and soil 
!     texture based on the work of Clapp and Hornberger (1978) and Cosby 
!     et al. (1984). The equation is integrated over the layer thickness,
!     in which the time rate of change in water mass must equal the net 
!     flow across the bounding interface, plus the rate of internal source
!     or sink.  The terms of water flow across the layer interfaces are 
!     linearly expanded by using first-order Taylor expansion.  The 
!     equations result in a tridiagonal system of equations.

!     Note: lentgh units here are all millimeter
!
!     Richards Equation
!
!     d wat     d      d wat   d psi
!     ----- = - -- [ k(-----   ----- - 1) ] + S
!      dt       dz       dz    d wat
!
!     where:
!     wat = volumetric water content (m^3/m^3)
!     psi = soil matric potential (mm)
!     dt  = time step (sec)
!     z   = depth (mm)
!     qin = inflow at top of layer (mm H2O/sec)
!     qout= outflow at bottom of layer (mm H2O/sec)
!     s   = source/sink flux (mm H2O/sec)
!     k   = hydraulic conductivity (mm H2O/sec)
!
!     Solution:
!     linearize k and psi about d wat and use tridiagonal system of 
!     equations to solve for d wat, where for layer j
!
!     r_j = a_j[d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
!     Revision History
!     15 September 1999: Yongjiu Dai; Initial code
!     15 December  1999: Paul Houser and Jon Radakovich; F90 Revision
!     25 January   2002: Ian Baker; SiB integration
!----------------------------------------------------------------------

use sibtype

use physical_parameters, only:  &
    grav,  &
    tice,  &
    hltm

use sib_const_module, only:  &
    nsoil,  &
    wimp,   &
    phmin,  &
    dti,    &
    dtt

implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  

!...INPUT VARIABLES
real(kind=dbl_kind),intent(in)  :: infil   ! infiltration into soil
                                             !  (kg m^-2 sec^-1) 
real(kind=dbl_kind),intent(in)  ::dzmm(1:nsoil)     
real(kind=dbl_kind),intent(in)  ::zmm(1:nsoil)


!...OUTPUT VARIABLES
real(kind=dbl_kind),intent(out) :: hk(1:nsoil)
                                             ! hydraulic conductivity
                                             !  (mm H2O sec^-1)    
real(kind=dbl_kind),intent(out) :: dhkdw(1:nsoil)   
                                             ! d(hk)/d(water)     
real(kind=dbl_kind),intent(out) :: dw_liq(1:nsoil)
                                             ! change in layer liquid
                                             !  (m^3/m^3)(volumetric)



!...local variables...
integer ::   j     ! loop variables

real(kind=dbl_kind)  :: hltmi     ! 1/hltm
real(kind=dbl_kind)  :: s_node    ! volumetric wetness of node
real(kind=dbl_kind)  :: s1        ! wetness at interface of layer
real(kind=dbl_kind)  :: s2        ! conductivity*wetness**(2b+2)
real(kind=dbl_kind)  :: smp(1:nsoil)
                                    ! soil matric potential (mm)
real(kind=dbl_kind)  :: dsmpdw(1:nsoil)
                                    ! d(smp)/d(wetness)
real(kind=dbl_kind)  :: qin       ! flux of water into layer 
                                    !  (mm H2O sec^-1)
real(kind=dbl_kind)  :: qout      ! flux of water out of layer 
                                    !  (mm H2O sec^-1)  
real(kind=dbl_kind)  :: den       ! used in calculating qin,qout   
real(kind=dbl_kind)  :: num       ! used in calculating qin,qout
real(kind=dbl_kind)  :: dqidw0    ! d(qin)/d(vol_liq(j-1)) 
real(kind=dbl_kind)  :: dqidw1    ! d(qin)/d(vol_liq(j))
real(kind=dbl_kind)  :: dqodw1    ! d(qout)/d(vol_liq(j))
real(kind=dbl_kind)  :: dqodw2    ! d(qout)/d(vol_liq(j+1))
real(kind=dbl_kind)  :: rmx(1:nsoil)
                                    ! "r" forcing term of tridiag matrix
real(kind=dbl_kind)  :: amx(1:nsoil)
                                    ! "a" left off diagonal of tridiag mat
real(kind=dbl_kind)  :: bmx(1:nsoil)
                                    ! "b" diagonal column for tridiag mat
real(kind=dbl_kind)  :: cmx(1:nsoil)
                                    ! "c" right off diagonal of tridiag mat



!----------------------------------------------------------------------

    hltmi = 1.0/hltm

    !...set hydraulic conductivity to zero if effective porosity 5% in 
    !...any two adjoining layers, or if volumetric water (liquid) content
    !...less than 0.001

    do j=1,nsoil
        if( ((sib%diag%eff_poros(j) < wimp)            .or.   &
            (sib%diag%eff_poros(min(nsoil,j+1)) < wimp)) .or. &
            (sib%prog%vol_liq(j) <= 1.E-3)) then
            hk(j)    = 0.0
            dhkdw(j) = 0.0
        else
            s1 = 0.5*(sib%prog%vol_liq(j)+sib%prog%vol_liq(min(nsoil,j+1)))/ &
                sib%param%poros

            s2 = sib%param%satco*1000.0*s1**(2.0*sib%param%bee+2.0)
            hk(j) = s1*s2

            dhkdw(j) = (2.0*sib%param%bee+3.0)*s2*0.5/sib%param%poros
            !itb...I don't like the 0.5 factor in the CLM code-don't understand.
            !itb...it's in the CLM manual-CLM does not follow Bonan exactly.


            if(j==nsoil) dhkdw(j) = dhkdw(j)*2.0
        endif

        !...evaluate hydraulic conductivity, soil matric potential,
        !...d(smp)/d(vol_liq) and d(hk)/d(vol_liq)
        if(sib%prog%td(j) > tice) then

            s_node = max(sib%prog%vol_liq(j)/sib%param%poros,0.01_dbl_kind)
            s_node = min(1.0_dbl_kind,s_node)

            smp(j) = sib%param%phsat*1000.0*s_node**(-sib%param%bee)

            smp(j) = max(phmin,smp(j))

            dsmpdw(j) = -sib%param%bee*smp(j)/(s_node*sib%param%poros)

        else

        !...when ice is present, the matric potential is only related to
        !...temperature by (Fuchs et. al., 1978; Soil Sci. Soc. of Amer. J.
        !...42(3); 379-385) 
        !...Unit 1 joule = 1kg m2/sec2 j/kg/(m/s2) ==> m ==>1e3mm

            smp(j) = 1.e3_dbl_kind * 0.3336e6_dbl_kind/grav *  &
                (sib%prog%td(j) - tice)/sib%prog%td(j)
            smp(j) = max(phmin,smp(j))
            dsmpdw(j) = 0.0

        endif
    enddo


    !...set up a, b, c and r vectors for tridiagonal solver
    !...node j=1

    j = 1
    qin    = infil
    den    = zmm(j+1) - zmm(j)
    num    = (smp(j+1)-smp(j)) - den
    qout   = -hk(j)*num/den
    dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
    dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
    rmx(j) = qin - qout  - (((sib%diag%ect * sib%param%rootr(j)) + &
        sib%diag%egs) * dti * hltmi)
    amx(j) = 0.0
    bmx(j) = dzmm(j) * (dti) + dqodw1
    cmx(j) = dqodw2

    !...nodes 2 through nsoil-1
    do j = 2,nsoil-1 
        den    = zmm(j) - zmm(j-1)
        num    = (smp(j)-smp(j-1)) - den
        qin    = -hk(j-1)*num/den
        dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
        dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den  
        den    = zmm(j+1)-zmm(j)
        num    = smp(j+1)-smp(j) - den
        qout   = -hk(j)*num/den
        dqodw1 = -(-hk(j)*dsmpdw(j)  +  num*dhkdw(j))/den
        dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
        rmx(j) = qin - qout  - (sib%diag%ect * dti * sib%param%rootr(j) * hltmi)
        amx(j) = -dqidw0
        bmx(j) = dzmm(j)*dti - dqidw1 + dqodw1
        cmx(j) = dqodw2
    enddo

    !...node j=nsoil
    j = nsoil
    den    = zmm(j) - zmm(j-1)
    num    = smp(j) - smp(j-1) - den
    qin    = -hk(j-1) * num/den
    dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
    dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den
    qout   = hk(j)
    dqodw1 = dhkdw(j)
    rmx(j) = qin - qout  - (sib%diag%ect * dti * sib%param%rootr(j) * hltmi)
    amx(j) = -dqidw0
    bmx(j) = dzmm(j)*dti - dqidw1 + dqodw1
    cmx(j) = 0.0

    ! Add qout out of the bottom layer to runoff
    sib%diag%roff = sib%diag%roff + qout*dtt
    
    !...solve
    call  clm_tridia (nsoil, amx, bmx, cmx, rmx, dw_liq)

end subroutine soilwater

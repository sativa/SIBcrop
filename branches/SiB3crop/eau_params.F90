!------------------------------------------------------------------------------

module eau_params
      
use kinds
use physical_parameters,  cpelq   => spec_heat_cp, &
                          gravelq => grav,         &
                          lcond   => hltm

implicit none
save
      
!------------------------------------------------------------------------------
!contains the constants needed for the eau_cup and eau_liq parameterizations of
!convection and cloud microphysics.
!Laura D. Fowler/slikrock (10-01-99).

!send comments to laura@atmos.colostate.edu.

!references:
!Fowler, L.D, D.A. Randall, and S.A. Rutledge, 1996: Liquid and Ice Cloud Micro
!physics in the CSU General Circulation Model: Model description and simulated
!cloud microphysical processes.
!J. Climate, 9, 489-529.

!Randall, D.A., and L.D. Fowler, 1999: EAUliq: The next generation. Dept. of
!Atmospheric Science Paper 673, Dept. of Atmospheric Science. Colorado State
!University, Fort Collins, Colorado, 65 pp.
!------------------------------------------------------------------------------

!eau_cup parameters:

integer (kind=int_kind), parameter:: &
   ltpcup = 1,                      &!index of highest convective cloud top.
   ncb    = 2                        !index of highest convective cloud base.

real (kind=dbl_kind), parameter:: &
   alpham = 1.e+08_dbl_kind,        &!
   ckemin = 5.0_dbl_kind,           &!minimum value of cke.
   taudis = 600.0_dbl_kind,         &!cke dissipation time scale.
   amiu   = 1.0_dbl_kind             !cloud-top entrainment.
     
!eau_liq parameters:

integer (kind=int_kind), parameter:: &
   nc1elq = 3                        !number of microphysics time-steps per
                                     !dynamical time-step.

real (kind=dbl_kind), parameter:: &
   a0elq    = -.267_dbl_kind,       &!
   a1elq    = 5.15e+03_dbl_kind,    &!
   a2elq    = -1.0225e+06_dbl_kind, &!
   a3elq    = 7.55e+07_dbl_kind,    &!
   alphaelq = 0.001_dbl_kind,       &!
   aprime   = 3.e+03_dbl_kind,      &!
   asecond  = 1.139_dbl_kind,       &!
   belq     = 0.11_dbl_kind,        &!
   betaelq  = 0.001_dbl_kind,       &!
   diffelq  = 2.26e-05_dbl_kind,    &!
   erccoef  = 1._dbl_kind,          &!
   esccoef  = 1._dbl_kind,          &!
   esicoef  = 0.1_dbl_kind,         &!
   gam3     = 2._dbl_kind,          &!
   gams1    = 2.21891_dbl_kind,     &!
   gams2    = 1.38784_dbl_kind,     &!
   gams3    = 6.90080_dbl_kind,     &!
   kap      = .2861328125_dbl_kind, &!
   muelq    = 1.718e-05_dbl_kind,   &!
   nzeror   = 8.e+06_dbl_kind,      &!
   nzeros   = 2.e+07_dbl_kind,      &!
   pielq    = 3.14159265,           &!
   pzero    = 1.e+05_dbl_kind,      &!
   qci0     = 0.01e-03_dbl_kind,    &!
!  qci0     = 0.1e-03_dbl_kind,     &!
   qcw0     = 0.25e-03_dbl_kind,    &!
!  qcw0     = 0.7e-03_dbl_kind,     &!
   rhor     = 1.e03_dbl_kind,       &!
   rhos     = 1.e02_dbl_kind,       &!
   taul     = 100._dbl_kind,        &!
   tauf     = 100._dbl_kind,        &!
   therco   = 2.43e-02_dbl_kind      !
      
!shared eau_liq and eau_cup parameters:

real (kind=dbl_kind), parameter:: &
   eauc0    = 0.0_dbl_kind,         &!
   eauc1    = 1.0_dbl_kind           !

real (kind=dbl_kind), parameter:: &
   lfus     = 0.3336e+06_dbl_kind,  &!
   lsub     = lcond+lfus,           &!
   t00      = 273.15_dbl_kind,      &!
   tbgmin   = 253.15_dbl_kind,      &!
   tbgmax   = 273.15_dbl_kind        !

!#ifdef eau_ng
!eau_ng parameters:
!real (kind=dbl_kind), parameter:: &
!   cClrCld  = eauc1,                &!exchange terms between cld and clr.
!   cCldClr  = eauc1                  !exchange terms between cld and clr.
!#endif


end module eau_params

!------------------------------------------------------------------------------

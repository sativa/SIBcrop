!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	
 	module physical_parameters

!-----------------------------------------------------------------
!	this modules specifices physical parameters and the units
!	  of those parameters (MKS is standard)
!-----------------------------------------------------------------

	use kinds

 	implicit none
 	save

	real (kind=dbl_kind), parameter :: &
     		pi           = 3.14159265358979323846_dbl_kind, &
            pidaypy      = 0.0172142_dbl_kind, &
     		a            = 6.37122E+06_dbl_kind ,&
     		grav         = 9.8100_dbl_kind,&
     		gravi        = 1.0_dbl_kind/grav,&
            omega        = 2.0_dbl_kind*pi/86400.0_dbl_kind,&
     		earth_area   = 5.100996990707616E+14_dbl_kind,&
     		gas_const_R  = 287.000_dbl_kind,&
     		spec_heat_cp = 1005.000_dbl_kind,&
     		kappa        = gas_const_R/spec_heat_cp,&
     		inv_kappa    = 1.0_dbl_kind / kappa,&
     		p0_sfc       = 1.0E+05_dbl_kind,&
     		inv_p0_sfc   = 1.0_dbl_kind / p0_sfc,&
            tice         = 273.15_dbl_kind,&
            hltm         = 2.25E+06_dbl_kind,&
            gamfac       = hltm*5417.9827_dbl_kind/spec_heat_cp,&
            delta        = 0.608_dbl_kind,&
            stefan       = 5.67e-08_dbl_kind,&
            rv           = 4.61e+02_dbl_kind

	real (kind=dbl_kind), parameter ::&
     		c0     = 0.00000_dbl_kind,&
     		c1     = 1.00000_dbl_kind,&
     		c2     = 2.00000_dbl_kind,&
     		c3     = 3.00000_dbl_kind,&
     		p5     = 0.50000_dbl_kind,&
     		p25    = 0.25000_dbl_kind,&
     		pi2    = 2*pi,&
     		dtr    = pi/180.0_dbl_kind,&
     		rtd    = 180.0_dbl_kind/pi,&
     		alpha2 = 0.0

!-------------------------------------------------------------------
!   VARIABLE DEFINITION
!-------------------------------------------------------------------
!	pi = pi
!	a = earth radius (m)
!	omega = earth angular velocity (1/s)
!	earth_area = surface area of earth (m2)
!	gas_const_R = gas constant for dry air (J/ (kg K))
!	spec_heat_cp = specific heat at constant pressure (J/ (kg K))
!	kappa = gas_const_R/spec_heat_cp
!	p0_sfc = surface pressure (Pa)
!	alpha2 = rotation of axis of rotation from NP/SP
!       tice = freezing temperature of water at 1 atm (K)
!       hltm = latent heat of vaporization (J/kg)
!       gamfac = a moist thermodynamic variable
!       delta = molecular_weight_air/molecular_weight_water - 1
!       inv_kappa = inverse kappa
!       stefan = stefan boltzmann constant (W/m^2/K^4)
!       rv = gas constant for water vapor
!
!	c0 = zero
!       c1 = one
!	c2 = two
!       p5 = 0.5
!       p25 = 0.25
!       pi2 = 2*pi
!       dtr = Degrees To Radians conversion constant
!	rtd = Radians To Degrees conversion constant
!-------------------------------------------------------------------

 	end module physical_parameters

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

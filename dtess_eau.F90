!------------------------------------------------------------------------------
subroutine dtess_eau (len, pl, tl, ess, dtess)
!------------------------------------------------------------------------------
!eau_sat computes the saturation mixing ratio, vapor pressure, and saturation,
!and their derivatives with respect to temperature over water, ice and mixed-
!phase clouds. the parameterization is that used in the Community Climate Com-
!munity Climate Model at NCAR.
!Laura D. Fowler /slikrock (07-01-01).

!send comments to laura@atmos.colostate.edu.

!subroutines called:
!none.

!argument list variables:
!input arguments:
!----------------

use kinds
use physical_parameters, only: hltm, rv
      
implicit none
      

integer (kind=int_kind), intent(in) :: len                !length of vector.

real(kind=dbl_kind), intent(in), dimension(len):: &
   pl,               &!pressure                                           (Pa).
   tl                 !temperature                                         (K).

!output arguments:
!-----------------

real (kind=dbl_kind), intent(out), dimension(len) :: &
   ess,              &!saturation vapor pressure                          (Pa).
   dtess              !derivative of es with respect to temperature     (Pa/K).

!local variables:

integer (kind=int_kind):: i

real (kind=dbl_kind):: &
   twmin=173.16,     &!lowest allowed temperature boundary for water       (K).
   twmax=373.16,     &!highest allowed temperature boundary for water      (K).     
   timin=173.16,     &!lowest allowed temperature boundary for ice         (K).
   timax=273.16,     &!highest allowed temperature boundary for ice        (K).
   tnull=273.16       !freezing temperature                                (K).

real (kind=dbl_kind) :: tstl , t0tl

real (kind=dbl_kind), dimension(len):: &
      esw ,     dtesw ,      esi ,     dtesi ,&
      esm ,     dtesm ,      tl0 ,            &
    wghtm 

!ccm parameterization of saturation vapor pressure over water and ice:

real (kind=dbl_kind), parameter:: &
   ps = 1013.246,                   &!reference pressure             (hPa).
   ts = 373.16,                     &!reference temperature            (K).
   t0 = 273.16,                     &!freezing temperature             (K)
   lsub = hltm+0.3336e+06_dbl_kind, &!
   tbgmin   = 253.15_dbl_kind,      &!
   tbgmax   = 273.15_dbl_kind        !
  
real (kind=dbl_kind):: &
       e1 ,   e2 ,     f ,    f1 ,&
       f2 ,   f3 ,    f4 ,    f5 ,&
   lphase , term , term1 , term2 ,&
   term3     

!------------------------------------------------------------------------------

!initialization of different arrays:

tl0    = tl
esw    = 0.0_dbl_kind
esi    = 0.0_dbl_kind
esm    = 0.0_dbl_kind
dtesw  = 0.0_dbl_kind
dtesi  = 0.0_dbl_kind
dtesm  = 0.0_dbl_kind

ess    = 0.0_dbl_kind
dtess  = 0.0_dbl_kind

!saturation over water:

do i = 1, len

   tl0(i)    = max(twmin,tl0(i))
   tl0(i)    = min(twmax,tl0(i))
   tstl      = ts / tl0(i)
   e1        = 11.344*(1.0 - tl0(i)/ts)
   e2        = -3.49149*(tstl - 1.0)
   f1        = -7.90298*(tstl - 1.0)
   f2        = 5.02808*log10(tstl)
   f3        = -1.3816*(10.0**e1-1.0)/10000000.0
   f4        = 8.1328*(10.0**e2-1.0)/1000.0
   f5        = log10(ps)
   f         = f1 + f2 + f3 + f4 + f5

   esw(i)    = (10.0**f)*1.e+02
   esw(i)    = min(esw(i),pl(i)*0.9)
   dtesw(i)  = hltm*esw(i)/(rv*tl0(i)*tl0(i))

   ess(i)    = esw(i)
   dtess(i)  = dtesw(i)


!saturation over ice:

   if(tl0(i).lt.timax) then

      tl0(i)    = max(tl0(i),timin)
      t0tl      = t0 / tl0(i)
      term1     = 2.01889049/(t0tl)
      term2     = 3.56654*log(t0tl)
      term3     = 20.947031*(t0tl)

      esi(i)    = 575.185606e10*exp(-(term1 + term2 + term3))
      esi(i)    = min(esi(i),pl(i)*0.9)
      dtesi(i)  = lsub*esi(i)/(rv*tl0(i)*tl0(i))

      ess(i)    = esi(i)
      dtess(i)  = dtesi(i)

   endif

!interpolated saturation variables:

   if(tl0(i).lt.tbgmax .and. tl0(i).ge.tbgmin) then

      wghtm(i)  = (tl0(i)-tbgmin)/(tbgmax-tbgmin)
      lphase    = hltm*wghtm(i)+lsub*(1.-wghtm(i))
      esm(i)    = wghtm(i)*esw(i) + (1.-wghtm(i))*esi(i)
      esm(i)    = min(esm(i),pl(i)*0.9)
      dtesm(i)  = lphase*esm(i)/(rv*tl0(i)*tl0(i))

      ess(i)    = esm(i)
      dtess(i)  = dtesm(i)

   endif
enddo

end subroutine dtess_eau

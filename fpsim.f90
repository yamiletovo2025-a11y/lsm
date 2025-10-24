      FUNCTION fpsim(zeta)

      IMPLICIT NONE
      REAL*8     zeta, fpsim, x
      REAL*8     cu, cs, pi
      PARAMETER (cu=16.d0, cs=5.d0, pi=3.141592654d0)

      IF(zeta.gt.0.d0) THEN
        fpsim=-cs*zeta
      ELSE
        x = (1.d0-cu*zeta)**0.25d0
        fpsim = dlog((1.+x**2.d0)*(1.d0+x)**2.d0/8.d0) -  &
     &          2.d0*datan(x) + pi/2.d0
      ENDIF
!
      END FUNCTION fpsim


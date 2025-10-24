      FUNCTION fpsis(zeta)

      IMPLICIT NONE
      
      REAL(8) :: zeta, fpsis, x
      REAL(8) :: cu, cs, pi
      PARAMETER (cu=16.d0, cs=5.d0, pi=3.141592654d0)

      IF(zeta.gt.0.d0) THEN
        fpsis=-cs*zeta
      ELSE
        x = sqrt(1.d0-cu*zeta)
        fpsis = 2.d0 * dlog((1.d0 + x) / 2.d0)

      ENDIF
!
      END FUNCTION fpsis
!




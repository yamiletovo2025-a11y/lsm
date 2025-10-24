      FUNCTION fave(x,ix)
	  
      REAL*8  x(ix), xmean, fave
      INTEGER ix, i

      xmean = 0.d0
      DO i = 1, ix 
         xmean = xmean + x(i)
      ENDDO
      fave = xmean/ix
!
      END FUNCTION fave

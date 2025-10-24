      FUNCTION FickianD(T, Theta, Thetas)

      REAL*8  FickianD, T, Theta, Thetas
      REAL*8  xi

      REAL*8, PARAMETER :: Tst=200.00d0, Dst=3.7d-5, xnd = 0.861d0
                        

      IF ( (Theta .lt. 0.) .or. (Theta .gt. Thetas) ) THEN
        print *, "Fickian diffusivity; Theta = ", Theta, "Thetas = ", Thetas 
        STOP
      ELSEIF (Thetas - Theta < 1.d-10) THEN 
        FickianD = 1.0d-12
      ELSE
        xi = 1/(Thetas- Theta)**0.5d0                                  
        FickianD =(Dst/xi)*(T/Tst)**xnd
      ENDIF
	  IF (.NOT. (FickianD == FickianD)) THEN
        PRINT *, "FickianD NaN"
        PRINT *, "T=", T, "Theta=", Theta, "Thetas=", Thetas
        PRINT *, "porosity_diff=", porosity_diff, "xi=", xi
        FickianD = 1.0d-12  
    ENDIF
	  
      END FUNCTION FickianD 


SUBROUTINE advel(xnu, xK, P0, T0, Tsj, Cvj, Cwj, wS, rz, rdz, m, wel, opt,iw,dz,hiwi,dP_inte,zm,dT_inte,dV_inte)

USE module_random
      IMPLICIT none
!
      INTEGER, parameter :: MXD = 50
      INTEGER :: i, m
      REAL(8) :: xK                   
      REAL(8) :: P0                 
      REAL(8) :: T0, Tp0
      REAL(8) :: wS
      REAL(8), dimension(MXD) :: xnu  
      REAL(8), dimension(MXD) :: Cvj  
	  REAL(8), dimension(MXD) :: dTdz,dVdz,dT_inte,dV_inte,dP_inte  
      REAL(8), dimension(MXD) :: Cwj  
      REAL(8), dimension(MXD) :: Tsj  
      REAL(8), dimension(MXD) :: dCO2 
      REAL(8) :: dCO20             
      REAL(8) :: delVol             
	  REAL(8) :: hydrostatic_term,density_term,temp_term              ! Air volume
      REAL(8), dimension(MXD) :: rz, rdz,dz,zm
      REAL(8), dimension(MXD) :: wel 
      REAL(8), PARAMETER :: Cf = 0.5  
      REAL(8), PARAMETER :: Rluft = 191.0d0      
      REAL(8), PARAMETER :: gacc  = 3.73d0        
      REAL(8), PARAMETER :: rhoCO2ice = 1625.d0    
	  REAL*8, PARAMETER :: Rgas = 8.3143d0  
      REAL(8) :: xwa, xwb, xwc, xwd, w1, w2, xwk
      REAL(8) :: xwc1, xwc2, xwc3
      INTEGER :: opt               
      INTEGER :: iw,ioerr                    
	  character*80 :: hiwi
	  REAL(8) :: p_hydrostatic(MXD)
		DO i = 2, m
			p_hydrostatic(i) =  Cvj(i) * gacc * (zm(i)-zm(i-1))
		ENDDO
		p_hydrostatic(1) =  Cvj(1) * gacc * (zm(1)-0)

      DO i = 1, m
         xwk = Cwj(i)/rhoCO2ice
         delVol  = wS - xwk                    
         dCO2(i) = log( Cvj(i)/delVol )
      ENDDO

	  dCO20 = log( exp( dCO2(1) ) )
      xnu(1) = 5.17d-4
      xwb = xnu(1)/(xK*Cvj(1))

	  xwc2 = -gacc

      xwc1 = Rluft*rdz(1)*( Tsj(1)-T0 )
      xwc3 = Rluft*Tsj(1)*rdz(1)*( dCO2(1) + (dCO2(2)- dCO2(1))*rz(1) -dCO20 )
	  dTdz(1) = xwc1*Cvj(1)
	  dVdz(1) = xwc3*Cvj(1)
    
      xwc  = xwc1 + xwc2 + xwc3
      IF ( opt == 1 ) THEN 
         wel(1) = -xwc/xwb                                                       
      ELSEIF ( opt == 2 ) THEN
         IF ( xwc <= 0.d0 ) THEN
            xwa = Cf/dsqrt(xK)
            wel(1) = ( -xwb + dsqrt( xwb*xwb - 4.d0*xwa*xwc ) )/( 2.d0*xwa )      
         ELSE
            xwa = -Cf/dsqrt(xk)
            wel(1) = ( -xwb + dsqrt( xwb*xwb - 4.d0*xwa*xwc ) )/( 2.d0*xwa )      
         ENDIF
      ENDIF
!
      DO i = 2, m-1

		 xwb = xnu(i)/(xK*Cvj(i))
  
         xwc1 = Rluft*rdz(i)*( Tsj(i)-Tsj(i-1) + (Tsj(i+1)-Tsj(i))*rz(i) - (Tsj(i)-Tsj(i-1))*rz(i-1) )
         xwc3 = Rluft*Tsj(i)*rdz(i)*(dCO2(i)-dCO2(i-1) + (dCO2(i+1)-dCO2(i))*rz(i) - &
              &                                          (dCO2(i)-dCO2(i-1))*rz(i-1) )
		 dTdz(i) = xwc1*Cvj(i)
	     dVdz(i) = xwc3*Cvj(i)
!------------------------------------------------------------------------------         
         xwc = xwc1 + xwc2 + xwc3
         IF ( opt == 1 ) THEN 
            wel(i) = -xwc/xwb                                                  
         ELSEIF ( opt == 2 ) THEN
            IF ( xwc <= 0.d0 ) THEN
               xwa = Cf/dsqrt(xK)
               wel(i) = ( -xwb + dsqrt( xwb*xwb - 4.d0*xwa*xwc ) )/( 2.d0*xwa )     
            ELSE
               xwa = -Cf/dsqrt(xk)
               wel(i) = ( -xwb + dsqrt( xwb*xwb - 4.d0*xwa*xwc ) )/( 2.d0*xwa )     
            ENDIF
         ENDIF
      ENDDO
	  xwb = xnu(m)/(xK*Cvj(m))
      xwc1 = Rluft*rdz(m)*( Tsj(m)-Tsj(m-1) - (Tsj(m)-Tsj(m-1))*rz(m-1) )
      xwc3 = Rluft*Tsj(m)*rdz(m)*(dCO2(m)-dCO2(m-1) - (dCO2(m)-dCO2(m-1))*rz(m-1) )
      xwc  = xwc1 + xwc2 + xwc3
	  dTdz(m) = xwc1*Cvj(m)
	  dVdz(m) = xwc3*Cvj(m)

      IF ( opt == 1 ) THEN 
         wel(m) = -xwc/xwb                                                     
      ELSEIF ( opt == 2 ) THEN
         IF ( xwc <= 0.d0 ) THEN
            xwa = Cf/dsqrt(xK)
            wel(m) = ( -xwb + dsqrt( xwb*xwb - 4.d0*xwa*xwc ) )/( 2.d0*xwa )       
         ELSE
            xwa = -Cf/dsqrt(xk)
            wel(m) = ( -xwb + dsqrt( xwb*xwb - 4.d0*xwa*xwc ) )/( 2.d0*xwa )       
         ENDIF
      ENDIF

	  
	  DO i = 1, m-1
		dT_inte(i) = 0.5 * (zm(i+1)-zm(i))*(dTdz(i)+dTdz(i+1))
		dV_inte(i) = 0.5 * (zm(i+1)-zm(i))*(dVdz(i)+dVdz(i+1))
	  ENDDO
	  
	  dP_inte(1) = P0
	  DO i = 2, m-1
			temp_term = (dT_inte(i-1) + 4.0 * dT_inte(i) + dT_inte(i+1)) / 6.0
			density_term = (dV_inte(i-1) + 4.0 * dV_inte(i) + dV_inte(i+1)) / 6.0
			dP_inte(i) = dP_inte(i-1) + (temp_term + density_term) * (zm(i+1) - zm(i-1)) / 2.0
			hydrostatic_term = p_hydrostatic(i)

	  ENDDO
		dT_inte(m) = 0.5 * (zm(m) - zm(m-1)) * (dTdz(m) + dTdz(m-1))
		dV_inte(m) = 0.5 * (zm(m) - zm(m-1)) * (dVdz(m) + dVdz(m-1))

		dP_inte(m) = dP_inte(m-1) + (dT_inte(m) + dV_inte(m)) * (zm(m) - zm(m-1))

	  DO i = 1, m
		dP_inte(i) = p_hydrostatic(i) + dP_inte(i)
	  ENDDO
      RETURN
    END SUBROUTINE advel

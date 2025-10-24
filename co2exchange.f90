SUBROUTINE co2exchange(Tsj, C_vj, C_wj, wS, rhoA, CO2A, z0CO2, gA, xK, fCO2, &
                       ECO2, Dfickian, mmax,qCO2_poreTop,sdelt,dz,temp)

      use module_random
      IMPLICIT none

      REAL(8) :: rhoA, CO2A, z0CO2, gA, xK
      REAL(8) :: FickianD        
      REAL(8) :: Vdiff           
      REAL(8) :: gc                
      REAL(8) :: deltTop            
      REAL(8) :: T_sub,LatentH,Cp
	  REAL(8) :: sdelt			
      REAL(8) :: fCO2              
      REAL(8) :: ECO2(mmax)      
      REAL(8) :: Tsj(mmax)       
      REAL(8) :: C_vj(mmax),dz(mmax)         
      REAL(8) :: C_wj(mmax)        
      REAL(8) :: Dfickian(mmax)     
      REAL(8) :: delta_Hm           
      REAL(8) :: alpha            
      REAL(8) :: Vth                
      REAL(8) :: NE              
      REAL(8), PARAMETER :: pi = 3.141592654d0
      REAL(8), PARAMETER :: mbar = 7.34d-26       
      REAL(8), PARAMETER :: rhoCO2ice = 1625.d0   
      REAL(8), PARAMETER :: Tsblim = 145.d0        
      REAL(8), PARAMETER :: T0 = 273.15d0        
      REAL(8), PARAMETER :: Rgas = 8.3143d0       
      REAL(8), PARAMETER :: delta_H = 591000.d0    
      REAL(8), PARAMETER :: P0 = 611.d0            
      REAL(8), PARAMETER :: NA = 6.022d+23         
      REAL(8), PARAMETER :: KB = 1.380649d-23      
      REAL(8) :: tau, wS, xwK, ss, eta, hcSl, L_char, L_s
      REAL(8) :: qCO2_z0, qCO2_poreTop, wf, f
      INTEGER :: jj, mmax, temp
      
      external :: FickianD
!
      xwk = C_wj(1)/rhoCO2ice
      Vdiff = FickianD(Tsj(1), xwk, wS)        

      ss = dsqrt( xK )                        
      gc = Vdiff/ss
      gc = ( gA**0.5 )*( gc**0.5 )

      deltTop  = wS - C_wj(1)/rhoCO2ice          
      qCO2_z0  = CO2A                           
      
      qCO2_poreTop = (1/rhoA)*(C_vj(1)/deltTop) 
!
      fCO2 = -gc*rhoA*(qCO2_z0 - qCO2_poreTop)    
      fCO2 = fCO2*(deltTop)       
!
      delta_Hm = delta_H*mbar*NA
      
	  T_sub = 145;
	  LatentH = 5.9d5;
	  Cp = 839.d0;
      DO jj = 1, mmax
		 xwk = C_wj(jj)/rhoCO2ice
         Vdiff = FickianD(Tsj(jj),xwk,wS)              
         Dfickian(jj) = Vdiff
         tau = xK/Vdiff                                   
         alpha = 5.d-11*(Tsj(jj)-Tsblim)/Tsblim

         ECO2(jj) = alpha*C_wj(jj)/tau                                           
		 IF (Tsj(jj)-Tsblim < 0.d0 .AND. C_vj(jj) > 1.d-10) THEN
             
              eta = 0.02d0
              hcSl = 0.035d0
              L_s = delta_H
              L_char = dz(jj)
              
              IF (L_s > 0.d0 .AND. L_char > 0.d0) THEN
                  ECO2(jj) = -eta * hcSl * ABS(Tsj(jj)-Tsblim) / (L_s * L_char)
                  
              ELSE
                  ECO2(jj) = 0.d0
              ENDIF
		ENDIF
      ENDDO
	  
     
      RETURN
END SUBROUTINE co2exchange

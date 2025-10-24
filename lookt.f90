SUBROUTINE lookt(wCO2ice, wS, hCpSl, hcSl, hDSl)
      implicit  none
      REAL*8    hCpSl,hcSl,hDSl,wCO2ice,wS
      REAL*8    CpQ, CpCO2ice, CpA, rhoQ, rhoCO2ice, rhoA,ksoil,kice,kco2

      PARAMETER (CpQ    = 732.d0)       
      PARAMETER (CpCO2ice = 839.d0)    
      PARAMETER (CpA    = 735.d0)       

      PARAMETER (rhoQ   = 1450.d0)       
      PARAMETER (rhoCO2ice = 1625.d0)    
      PARAMETER (rhoA   = 0.020)      
	  PARAMETER (ksoil   = 0.035)  
	  PARAMETER (kice   = 0.65)   
	  PARAMETER (kco2   = 0.0059)  

      hCpSl = (1.d0-wS)*rhoQ*CpQ + wCO2ice*rhoCO2ice*CpCO2ice &
            & +(1.d0-wS-wCO2ice)*rhoA*CpA
	  hcSl = (1.d0-wS)*ksoil + wCO2ice*kice + (1.d0-wS-wCO2ice)*kco2

      hDSl = hcSl/hCpSl  

      END SUBROUTINE lookt


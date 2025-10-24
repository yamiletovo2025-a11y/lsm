
      SUBROUTINE BareSoil(Ebs, T, Ts, RsDnS, RlDnS,       &
     &     Tair, rho, k_sj, delz, emisss, RlUpS, Rnets,   &
     &     SHs,  fhG, gas)


      REAL*8 Ebs
      REAL*8 T, Ts, Tair
      REAL*8 rho, k_sj, delz, emisss
      REAL*8 Rnets, SHs, fhG
      REAL*8 RsDnS, RlDnS, RlUpV, RlUpS
      REAL*8 gas
!
      REAL*8, PARAMETER :: Tc0 = 273.15d0
      REAL*8, PARAMETER :: SBoltz = 5.67d-8    
      REAL*8, PARAMETER :: CpAir  = 735.d0        
!
      SHs  = rho * CpAir * (T - Tair) * gas      
      fhG = -2.d0 * k_sj * (Ts - T) / delz
!
      RlUpS = emisss * SBoltz * T**4.d0
      RnetS = RsDnS + RlDnS - RlUpS
      Ebs = RnetS - SHs - fhG
      
      end subroutine BareSoil

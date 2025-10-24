     PROGRAM mars_alsis
      use module_random
	  USE IEEE_ARITHMETIC

      IMPLICIT none

      INTEGER, parameter :: MXD = 50             
      REAL(8)            :: wS                  

      REAL(8)            :: zC
      REAL(8), dimension(MXD) :: dz, rdz, rz, cz, zz, zm
      REAL(8)            :: r2dt, rdt, hdt
      REAL(8)            :: zeta, rnd, Amp
	  REAL(8)            :: temp               

      INTEGER, parameter :: nmax=200000           
      INTEGER            :: m,intemp	      
      REAL(8), dimension(nmax) :: SolDn, LwDn, Uair, Tair, Pres, QCO2
      
      REAL(8) :: Tr, delt
      INTEGER :: irdr, iprt, idata, debug, itemp  

      REAL(8), PARAMETER :: CpAir  = 735.d0         ! Mars air spec heat at const P [J/kg/K]
      REAL(8), PARAMETER :: CvAir  = 544.d0         ! Mars air spec heat at const V [J/kg/K]
      REAL(8), PARAMETER :: CpCO2  = 839.d0         ! CO2 spec heat                 (J/kg/K), for CO2 gas, not found for CO2 ice
      REAL(8), PARAMETER :: CpW    = 4180.d0        ! water spec heat               [J/kg/K]
      REAL(8), PARAMETER :: CpQ    = 732.d0         ! Mars quartz spec heat         [J/kg/K]
      REAL(8), PARAMETER :: Lsblim = 591000.d0      ! Enthalpy of sublimation at 180K [J/kg]
      REAL(8), PARAMETER :: hcSoil = 0.04d0         ! Mars soil thermal conductivity (W/m/K), Grott etal. doi.org/10.1029/2021JE006861

      REAL(8), PARAMETER :: Xfreepath = 1.d-6       ! Mars air molecular free path, 10 um (DOI: 10.5772/intechopen.90912)
      REAL(8), PARAMETER :: airdynvis = 9.82d-6     ! Mars air dynamic viscosity, N s/m2  (DOI: 10.5772/intechopen.90912)
      REAL(8), PARAMETER :: airkinvis = 5.17d-4     ! Mars air kinematic viscosity, m2/s  (DOI: 10.5772/intechopen.90912)
!
      REAL(8), PARAMETER :: Dsoil = 0.000037d0      ! Hudson et al. 2007: Diffusivity range: 2.8 - 5.4cm2/s at 600 Pa & 260 K,
                                                    ! 1.9–4.7cm2/s to a Mars temp of 200K. Preferred value for JSC Mars–1 at 600 Pa
                                                    ! and 200K is 3.7±0.5cm2/s. Packed dust displays lower diffusivity
                                                    ! 0.38±0.26 cm2/s;
      REAL(8), PARAMETER :: Turt = 1.8d0            ! Turtosity = 1.8
      
      REAL(8), PARAMETER :: Cf_Ward = 0.5           ! Ward (1964) set Cf_Ward=0.55; Studies suggest varies (Lu, 2020, Thesis)
!                                                   
      REAL(8), PARAMETER :: Pmib  = 1.0D-12         ! Mars soil permiability mdarcy 1e-12 (m2) in literature
      REAL(8), PARAMETER :: gacc  = 3.73d0          ! Mars gravity acceleration     [m/s2]
      REAL(8), PARAMETER :: S0    = 586.2           ! Mars solar constant           [W/m2]
      REAL(8), PARAMETER :: mAir  = 0.04323d0       ! Mars air mol weight      [kg/mol]
      REAL(8), PARAMETER :: mCO2  = 0.044009d0      !              water       [kg/mol]
      REAL(8), PARAMETER :: mH2O  = 0.018016d0      !              water       [kg/mol]

      REAL(8), PARAMETER :: pi    = 3.141592654d0

      REAL(8), PARAMETER :: rhoW  = 1000.d0        ! water density            [kg/m3]
      REAL(8), PARAMETER :: rhoQ  = 3200.d0        ! Mars soil particle density  [kg/m3] (doi.org/10.1073/pnas.0800202105)
      REAL(8), PARAMETER :: rden  = 160000.d0      ! Mars Particle-to-air density ratio
      REAL(8), PARAMETER :: rhoA  = 0.020d0        ! Mars default air
      REAL(8), PARAMETER :: rhoACO2   = 0.0193d0   ! CO2 density in Mars atmosphere
      REAL(8), PARAMETER :: rhoCO2ice = 1625.d0    ! Mars CO2 ice density     [kg/m3]
      REAL(8), PARAMETER :: Rgas  = 8.3143d0       ! universal gas constant   [J/mol/K]
      REAL(8), PARAMETER :: Rluft = 191.0d0        ! Mars specific air gas constant [J/mol/K]
      REAL(8), PARAMETER :: Rv    = 461.0d0        ! specific water vapor gas constant [J/K/kg]
      REAL(8), PARAMETER :: SBoltz = 5.67d-8        ! Stefan-Boltzmann const   [W/m2/K4]
      REAL(8), PARAMETER :: tday  = 88775.23104d0  ! Mars sol (solar days) = 24.6597 hrs    [s]
      REAL(8), PARAMETER :: tyear = 59301854.33472d0    ! 1 Mars year = 668 sols = 687 Earth days
      REAL(8), PARAMETER :: vonk  = 0.4d0          ! von Karman constant
      REAL(8), PARAMETER :: alf   = 3.33d+5        ! 
      REAL(8), PARAMETER :: Tc0   = 273.15d0       !

      CHARACTER*200  OutFile, hdr1, hdr2, hdr3

      INTEGER       nyear, nday, ntime
      INTEGER       iout, i, nrep, ifst, ilst

      REAL*8 :: vDiff                                        
      REAL*8 :: ss                                        
      REAL*8 :: sdelt, shdt, srdt                             
      INTEGER :: nsrep
     real(8), dimension(nmax, MXD) :: STsoil       
     real(8), dimension(nmax, MXD) :: SPsoil        
     real(8), dimension(nmax, MXD) :: SWsoil     
     real(8), dimension(nmax)      :: SolAbs      
     real(8), dimension(nmax)      :: Rnet          
     real(8), dimension(nmax)      :: SH             
     real(8), dimension(nmax)      :: FTsoil0          
     real(8), dimension(nmax)      :: Trad            
     real(8), dimension(nmax)      :: Rair            
     real(8) :: RsDnS, RlDnS, RlUpS, RnetS, albedS
     real(8) :: emisiS

     real(8) :: ustar, usb,ita
     real(8) :: gas             
     real(8) :: Eb, EbL, EbH  
     
     real(8), dimension(nmax, MXD) :: CVsoil          
     real(8), dimension(nmax, MXD) :: ECO2soil,Cp_soil,decrease_T  
     real(8), dimension(nmax, MXD) :: CWsoil        
     real(8), dimension(nmax, MXD) :: Vsoil,Psoil,dTsoil,dVsoil           
     real(8), dimension(nmax)      :: FCO2soil, QCO2soil        

     real(8), dimension(MXD)       :: ECO2             
     real(8), dimension(MXD)       :: ECO2_mst,Cp_sj1_mst         
     real(8), dimension(MXD)       :: vdelta,dP_inte,dT_inte,dV_inte          
     real(8), dimension(MXD)       :: vdelta_mst,P_mst,T_mst,V_mst       
     real(8), dimension(MXD)       :: T_sj1, T_sj
     real(8), dimension(MXD)       :: C_vj1, C_vj, C_vj_initial
     real(8), dimension(MXD)       :: C_wj1, C_wj
     real(8)                       :: fCO2, QCO2_pt           
     real(8)                       :: fCO2_mst, QCO2_mst       
     real(8)                       :: fCO2limit1, fCO2limit2

     real(8), dimension(MXD)       :: CVsoil_ave      
     real(8), dimension(MXD)       :: CWsoil_ave       
     real(8), dimension(MXD)       :: STsoil_ave     
   
     real(8), dimension(MXD)       :: aCv, bCv, cCv, rCv   
     real(8), dimension(MXD)       :: aT,  bT,  cT,  rT    
     real(8)                       :: bbj, bbj1, bbCvj
     real(8)                       :: aaa, bbb
          
     real(8), dimension(MXD)       :: Dfickian, d_sj, d_sj1 
     real(8), dimension(MXD)       :: deltj, deltj1        
     real(8), dimension(MXD)       :: Cp_sj1, Cp_sj       
     real(8), dimension(MXD)       :: k_sj1, k_sj          

     real(8) :: fave, fpsim, fpsis
     real(8) :: FickianD
     real(8) :: zRef, z0m, z0s, fhG

     real(8) :: TrH, TrL
     real(8) :: SHs
     
     real(8) :: ls, xwk1, xwk2, xwc, xwb
	 character(len=80) :: utc, dataID
     real(8) :: dummyhc, dummyhd
     
     INTEGER nn, nsmall, j, jj, ist, ifn
     INTEGER, PARAMETER :: ierror = 0     
     INTEGER idump, isign, iw

     character(200) :: DataFile, hiwi, hiwi1, hiwi2
     REAL    xx1,xx2,xx3,xx4,xx5,xx6

      irdr  = 10
      iprt  = 11
      idata = 12
      iout  = 13
	  debug = 58
      OPEN (irdr, file="mars_alsis.rdr")
      OPEN (iprt, file="mars_alsis.prt")

      READ (irdr, *) nrep          
      READ (irdr, *) wS           
      READ (irdr, *) m             
	  READ (irdr, *) ita             

      READ (irdr, *)                
      READ (irdr, *)                

      DO i = 1, m
         READ (irdr, *) dz(i), T_sj(i), C_vj(i), C_wj(i)
		 T_sj(i) = T_sj(i)*ita
		 C_vj_initial = C_vj(i)
      ENDDO
!
      READ (irdr, *) Tr                                   
      READ (irdr, *) delt       
      READ (irdr, *) ifst, ilst 
      READ (irdr, *) DataFile  
	  READ (irdr, *) dataID   
      READ (irdr, *) OutFile   
      READ (irdr, *) zRef     
      READ (irdr, *) z0m       
      READ (irdr, *) albedS    
      READ (irdr, *) emisiS     

      READ(irdr, *) nyear      
      READ(irdr, *) nday      
      READ(irdr, *) ntime      
      
      READ(irdr, *) hdr1       
      READ(irdr, *) hdr2    
      READ(irdr, *) hdr3        
      hiwi = TRIM(OutFile) // "_" // TRIM(dataID)
      OPEN (50, file=hiwi)
      write(50, *) ilst, nmax
      write(50, *) m, MXD
      write(50, *) nyear,nday, ntime, delt
      write(50, *) hdr1
      write(50, *) hdr2
      write(50, *) hdr3
      write(50, *) "Mars atmospheric and surface variable"
      write(50,'(7a12)') "Time", "SolAbs", "Rnet", "SH",  "FTsoil0", "Trad", "Rair"
      write(50,'(7a12)') "iw  ", "W/m2",   "W/m2", "W/m2","W/m2",    "K",    "m2/s"

      hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_STsoil" 
      OPEN (51, file=hiwi)
      write(51, *) ilst, nmax
      write(51, *) m, MXD
      write(51, *) nyear,nday, ntime, delt
      write(51, *) hdr1
      write(51, *) hdr2
      write(51, *) hdr3
      write(51, *) "Soil temperature of different layers"

      hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_CVsoil" 
      OPEN (52, file=hiwi)
      write(52, *) ilst, nmax
      write(52, *) m, MXD
      write(52, *) nyear,nday, ntime, delt
      write(52, *) hdr1
      write(52, *) hdr2
      write(52, *) hdr3
      write(52, *) "volumetric CO2 gas concentration in soil"

	  hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_CWsoil" 
      OPEN (53, file=hiwi)
      write(53, *) ilst, nmax
      write(53, *) m, MXD
      write(53, *) nyear,nday, ntime, delt
      write(53, *) hdr1
      write(53, *) hdr2
      write(53, *) hdr3

	  hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_ECO2soil" 
      OPEN (54, file=hiwi)
      write(54, *) ilst, nmax
      write(54, *) m, MXD
      write(54, *) nyear,nday, ntime, delt
      write(54, *) hdr1
      write(54, *) hdr2
      write(54, *) hdr3

      hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_Vsoil" 
      OPEN (55, file=hiwi)
      write(55, *) ilst, nmax
      write(55, *) m, MXD
      write(55, *) nyear,nday, ntime, delt
      write(55, *) hdr1
      write(55, *) hdr2
      write(55, *) hdr3

	  OPEN(120, file="FCO2_surface.prt")
	  write(120,'(4a13)') "Time ", "QCO2_soilsurf", "QCO2_air", "delta"
	  hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "Cp_soil.prt"
	  OPEN(121, file=hiwi)
      write(121,'(6a12)') "Time", "Cpsoil_1", "Cpsoil_2", "Cpsoil_3", "Cpsoil_4", "Cpsoil_5"
      write(121,'(2a12)') "iw  ", "J/kg"
	  hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "dec_T.prt"
	  OPEN(122, file=hiwi)
      write(122,'(6a12)') "Time", "dec_T1", "dec_T2", "dec_T3", "dec_T4", "dec_T5"
      write(122,'(6a12)') "iw  ", "K",    "K",      "K",      "K",      "K",      "K"
	  
	  hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_SPsoil"
	  OPEN(56, FILE=hiwi)

	  hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_dT"
	  OPEN(80, FILE=hiwi)
	  hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_dV"
	  OPEN(81, FILE=hiwi)

      IF ( (ifst.lt.1).or.(ifst.gt.ilst)) THEN
        WRITE(iprt, *) 'ifst = ', ifst, ' incorrect'
        STOP
      ENDIF
!
      IF ( (ilst.lt.1).or.(ilst.gt.nmax)) THEN
        WRITE(iprt, *) 'ilst = ', ilst, ' incorrect',nmax
        STOP
      ENDIF
!
      zC = 0.d0
      DO i = 1, m
         zC = zC + dz(i)
         zz(i) = zC
         zm(i) = zz(i) - 0.5*dz(i)
         rdz(i)  = 1.d0/dz(i)
      ENDDO

      DO i = 1, m-1
        rz(i) = dz(i)/( dz(i)+dz(i+1) )
        cz(i) = 2.d0/( dz(i)+dz(i+1) )
      ENDDO
      rz(m) = dz(m)/( dz(m)+dz(m-1) )     
      cz(m) = 0.d0
!
      r2dt = 0.5d0/delt            
      rdt = 0.1d0/delt             
      hdt = 0.001d0*delt

      sdelt = delt/120000.d0
      srdt  = 1.d0/sdelt
      shdt = 0.2d0*sdelt
      nsrep = int( delt/sdelt)

      hiwi = DataFile( 1:len(DataFile) )
      OPEN (idata, file=hiwi)

      READ(idata, *) 
	  temp = 0
	  DO i = 1, ilst
		READ(idata, *) utc, ls, xx1, xx2, xx3, xx4, xx5
		SolDn(i) = xx1
		LwDn(i)  = xx2
		Tair(i)  = xx3 * ita
		Uair(i)  = xx4
		Pres(i)  = xx5
	  ENDDO

      zeta = 0.d0


      ist = 1
      ifn = m
	  intemp = 0;

      DO i = 1, m
         xwk2 = C_wj(i)/rhoCO2ice                    
         CALL lookT(xwk2,wS,Cp_sj(i),k_sj(i),dummyhd)
      ENDDO

      DO 100 nn = 1, nrep

      print *, "nn = ", nn         
      IF (nn > 1) THEN
        T_sj(:) = STsoil_ave(:)
        C_vj(:) = CVsoil_ave(:)
        C_wj(:) = CWsoil_ave(:)
        write(201,*) (T_sj(i), i=1,m)
        write(202,*) (C_vj(i), i=1,m)
        write(203,*) (C_wj(i), i=1,m)
      ENDIF

      DO 200 iw = ifst, ilst

         z0s = z0m / 7.d0 
    
         aaa = (dlog(zRef / z0m) - fpsim(zeta))
         bbb = (dlog(zRef / z0s) - fpsis(zeta))
         
         ustar = vonk * Uair(iw) / aaa 
         gas   = (vonk * ustar)  / bbb 

         ustar = dmax1(ustar, 0.01d0)  
         gas   = dmax1(gas,   0.001d0)

         RsDnS = SolDn(iw) * (1.d0 - albedS)
         RlDnS = LwDn(iw)

         TrH = Tair(iw) + 15.d0
         TrL = Tair(iw) - 15.d0
!
         CALL BareSoil(EbH, TrH, T_sj(1), RsDnS, RlDnS,      &
             &     Tair(iw), rhoA, k_sj(1), dz(1), &
             &     emisiS, RlUpS, Rnets, SHs, fhG, gas)
         CALL BareSoil(EbL, TrL, T_sj(1), RsDnS, RlDnS,      &
             &     Tair(iw), rhoA, k_sj(1), dz(1), &
             &     emisiS, RlUpS, Rnets, SHs, fhG, gas)
         IF (  EbH * EbL .gt. 0.d0) THEN 
            IF (ierror == 1) THEN  
               Write(iprt,*)'EbL, gas, T_sj(1), RsDnS, RlDnS'
               Write(iprt,*) EbL, gas, T_sj(1), RsDnS, RlDnS
               Write(iprt,*)'Tair(iw),Pres(iw)'
               Write(iprt,*) Tair(iw),Pres(iw)
               Write(iprt,*)'rhoA, Lsblim, k_sj(1), dz(1), emisiS'
               Write(iprt,*) rhoA, Lsblim, k_sj(1), dz(1), emisiS
               Write(iprt,*)'RlUpS, Rnets, SHs, fhG'
               Write(iprt,*) RlUpS, Rnets, SHs, fhG
               WRITE(iprt,*)'Tr: F(Tr) HAVE SAME SIGN'
               idump = idump + 1
            ENDIF
         ENDIF

         IF ( (EbH - EbL) .eq. 0.d0) THEN
            WRITE(iprt,*) 'Tr: F(TrH)=F(TrL)'
            idump = idump + 1
         ENDIF
!
         ISIGN = 1
         IF (EbH .LT. EbL) ISIGN = -1
         EbH = ISIGN * EbH
         EbL = ISIGN * EbL
!c
         j = 0
123      j = j + 1
         Tr = ( TrH + TrL ) * 0.5d0
!
         CALL BareSoil(Eb, Tr, T_sj(1), RsDnS, RlDnS,   &
              &     Tair(iw), rhoA, k_sj(1), dz(1),     &
              &     emisiS, RlUpS, Rnets, SHs, fhG, gas)  
!
         Eb = ISIGN * Eb

         IF (Eb .gt. EbH .or. Eb .lt. EbL) THEN
            WRITE(iprt,*) 'Tr: F(Tr) NON-MONOTONIC'
            idump = idump + 1
         ENDIF
!
         IF ((Eb.EQ.0.d0).or.((TrH-TrL).le.1d-3)) GOTO 124
         IF (j .le. 200000) THEN
            IF (Eb .GT. 0.d0) TrH = Tr
            IF (Eb .LT. 0.d0) TrL = Tr
            GO TO 123
         ELSE
            WRITE(iprt, *) 'WARNING: Tr iteratn unconverge'
            WRITE(iprt, *) 'Tr = ', Tr
            WRITE(iprt, *) 'RnetS = ', RnetS
            WRITE(iprt, *) 'SHs = ', SHs
            WRITE(iprt, *) 'fhG = ', fhG
         ENDIF
124      CONTINUE

         zeta = - zRef * (vonk * gacc * SHs) / (rhoA * CpAir *&
              & Tair(iw) * ustar ** 3.d0)
         IF ( zeta .gt. 10.d0 ) zeta = 10.d0
         IF ( zeta .lt.-10.d0 ) zeta = -10.d0

         DO i = 1, m
            xwk2 = C_wj(i)/rhoCO2ice                  
            CALL lookT(xwk2,wS,Cp_sj1(i),k_sj1(i),dummyhd)
         ENDDO
     
         DO i = 2, m-1
            bbj1 = hdt / Cp_sj1(i) * rdz(i)
            aT(i) = - bbj1 * cz(i) * (k_sj1(i) + (k_sj1(i+1) -&
                 & k_sj1(i)) * rz(i))
            cT(i) = -bbj1 * cz(i-1) * (k_sj1(i-1) + (k_sj1(i) -&
                 & k_sj1(i-1)) * rz(i-1))
            bT(i) = 1.d0 - aT(i) - cT(i)

            bbj = hdt / Cp_sj(i) * rdz(i)
            rT(i) = T_sj(i) + bbj * ( cz(i) * (T_sj(i+1) - T_sj(i)) *&
                 & (k_sj(i) + (k_sj(i+1) - k_sj(i)) * rz(i)) - cz(i&
                 &-1) * (T_sj(i) - T_sj(i-1)) * (k_sj(i-1) + (k_sj(i)&
                 & - k_sj(i-1)) * rz(i-1)))
         ENDDO
     
         bbj1 = hdt / Cp_sj1(1) * rdz(1)
         aT(1) = -bbj1 * cz(1) * ( k_sj1(1) + (k_sj1(2) - k_sj1(1)) *&
              & rz(1) )
         cT(1) = 0.d0
         bT(1) = 1.d0 - aT(1)
     
         bbj = hdt / Cp_sj(1) * rdz(1)
         rT(1) = T_sj(1) + bbj1 * fhG + bbj * (cz(1) * (T_sj(2) -&
              & T_sj(1)) * (k_sj(1) + (k_sj(2) - k_sj(1)) * rz(1)) +&
              & fhG )

     
         bbj1  = hdt / Cp_sj1(m) * rdz(m)
         aT(m) = 0.d0
         cT(m) = -bbj1 * cz(m-1) * ( k_sj1(m-1) + (k_sj1(m) - k_sj1(m&
              &-1)) * rz(m-1) )
         bT(m) = 1.d0 - cT(m)

         bbj   = hdt / Cp_sj(m) * rdz(m)
         rT(m) = T_sj(m) - bbj * cz(m-1) * (T_sj(m) - T_sj(m-1)) *&
              & (k_sj(m-1) + (k_sj(m) - k_sj(m-1)) * rz(m-1))

         CALL tridiagonal(ist, ifn, cT, bT, aT, rT, T_sj1, ierror)

         IF (iw .eq. 1 ) THEN
			xwk1 = 0.966197                                          
		 ENDIF
		 IF (iw .gt. 1) THEN
		    xwk1 = QCO2(iw-1) + FCO2soil(iw-1)/43.54492
		 ENDIF
         ECO2_mst(:)   = 0.d0
         vdelta_mst(:) = 0.d0
		 P_mst(:) = 0.d0
		 Cp_sj1_mst(:) = 0.d0
         fCO2_mst      = 0.d0
		 QCO2_mst      = 0.d0
         DO nsmall = 1, nsrep

			intemp = intemp + 1
            CALL co2exchange(T_sj,C_vj,C_wj,wS,rhoA,xwk1,z0s,gas,Pmib,fCO2,ECO2,Dfickian,m,QCO2_pt,sdelt,dz,intemp)

            DO i = 1, m
               C_wj1(i) = C_wj(i)*exp( -ECO2(i)*sdelt )
            ENDDO

		hiwi = TRIM(OutFile) // "_" // TRIM(dataID) // "_SP_EACH"
            call advel(Dfickian,Pmib,Pres(iw),Tair(iw),T_sj,C_vj,C_wj,wS,rz,rdz,m,vdelta,2,iw,dz,hiwi,dP_inte,zm,dT_inte,dV_inte)       

         DO i = 1, m
            xwk2 = C_wj1(i)/rhoCO2ice                    
            CALL lookT(xwk2,wS,Cp_sj1(i),k_sj1(i),dummyhd)
            T_sj1(i) = T_sj1(i) - Lsblim*ECO2(i)*sdelt/Cp_sj1(i)
			
         ENDDO

         DO jj = 1, m
            deltj1(jj) = wS - C_wj1(jj)/rhoCO2ice
            deltj (jj) = wS - C_wj(jj) /rhoCO2ice
            xwk2 = C_wj1(i)/rhoCO2ice                    
            vDiff = FickianD(T_sj1(jj),xwk2,wS)
            d_sj1(jj)  = vDiff*deltj1(jj)               
            d_sj (jj)  = vDiff*deltj (jj)
         ENDDO

         DO i = 2, m-1
            bbCvj = shdt*rdz(i)
            aCv(i) = - bbCvj * cz(i) * (d_sj1(i) + (d_sj1(i+1) -  &
                 &        d_sj1(i)) * rz(i))
            cCv(i) = -bbCvj * cz(i-1) * (d_sj1(i-1) + (d_sj1(i) - &
                 &        d_sj1(i-1)) * rz(i-1))
            bCv(i) = 1.d0 - aCv(i) - cCv(i)
!
            rCv(i) = C_vj(i) + bbCvj * ( cz(i) * (C_vj(i+1) -     &
                 &        C_vj(i)) * (d_sj(i) + (d_sj(i+1) - d_sj(i))      &
                 &        * rz(i)) - cz(i-1) * (C_vj(i) - C_vj(i-1)) *     &
                 &        (d_sj(i-1) + (d_sj(i) - d_sj(i-1)) * rz(i-1))) 
         ENDDO

    
         bbCvj  = shdt * rdz(1)
         aCv(1) = -bbCvj * cz(1) * ( d_sj1(1) + (d_sj1(2) -      &
              &     d_sj1(1)) * rz(1) )
         cCv(1) = 0.d0
         bCv(1) = 1.d0 - aCv(1)    
         rCv(1) = C_vj(1) - bbCvj*fCO2 + bbCvj*( cz(1)*          &     
              &     ( C_vj(2)-C_vj(1) )*( d_sj(1)+( d_sj(2) -    &
              &     d_sj(1) )*rz(1) ) )                                 


         bbCvj  = shdt * rdz(m)
         aCv(m) = 0.d0
         cCv(m) = -bbCvj * cz(m-1) * ( d_sj1(m-1) + (d_sj1(m) -  &
              &     d_sj1(m-1)) * rz(m-1) )
         bCv(m) = 1.d0 - cCv(m)

         rCv(m) = C_vj(m) - bbCvj * cz(m-1) * (C_vj(m) - C_vj(m-1)) * &
              &     (d_sj(m-1) + (d_sj(m) - d_sj(m-1)) * rz(m-1))

         CALL tridiagonal(ist,ifn,cCv,bCv,aCv,rCv,C_vj1,ierror)

         DO i = 1, m-1
            C_vj1(i) = C_vj1(i) - sdelt*cz(i)*( vdelta(i+1)*C_vj1(i+1) - vdelta(i)*C_vj1(i) ) 
         ENDDO
         C_vj1(m) = C_vj1(m) - sdelt*cz(m)*( vdelta(m)*C_vj1(m) - vdelta(m-1)*C_vj1(m-1) )      
         
         DO i = 1, m
            IF ( C_vj1(i) < 0.d0 ) THEN
               C_vj1(i) = 0.d0
             ENDIF

         ENDDO

         DO i = 1, m
            C_vj1(i) = C_vj1(i) + ECO2(i)*sdelt
         ENDDO
         !
          DO i = 1, m
            IF ( C_vj1(i) < 0.d0 ) THEN
               C_vj1(i) = 0.d0
            ENDIF

         ENDDO

         bbCvj  = shdt * rdz(1)
         fCO2limit1 = 0.5*C_vj1(1)/bbCvj
         fCO2limit2 = -C_vj1(1)/bbCvj
         fCO2 = min(fCO2, fCO2limit1)
         fCO2 = max(fCO2, fCO2limit2)
         C_vj1(1) = C_vj1(1) - bbCvj*fCO2

         IF ( C_vj1(1) < 0.d0 ) THEN
            C_vj1(1) = 0.d0
         ENDIF

         C_vj(:) = C_vj1(:)
		 C_wj(:) = C_wj1(:)
         ECO2_mst(:) = ECO2_mst(:) + ECO2(:)
         vdelta_mst(:) = vdelta_mst(:) + vdelta(:)
		 P_mst(:) = P_mst(:) + dP_inte(:)
		 T_mst(:) = T_mst(:) + dT_inte(:)
		 V_mst(:) = V_mst(:) + dV_inte(:)
         fCO2_mst = fCO2_mst + fCO2
		 QCO2_mst = QCO2_mst + QCO2_pt
		 Cp_sj1_mst(:) = Cp_sj1_mst(:) +Cp_sj1(:)
      ENDDO                
      
      ECO2_mst(:) = ECO2_mst(:)/nsrep                   
      vdelta_mst(:) = vdelta_mst(:)/nsrep               
	  P_mst(:) = P_mst(:)/nsrep
	  T_mst(:) = T_mst(:)/nsrep
	  V_mst(:) = V_mst(:)/nsrep
      fCO2_mst = fCO2_mst/nsrep
	  QCO2_mst = QCO2_mst/nsrep
	  Cp_sj1_mst = Cp_sj1_mst/nsrep
   
      DO i = 1, m
         T_sj(i) = T_sj1(i)                   
         C_vj(i) = C_vj1(i)   
         C_wj(i) = C_wj1(i)   
      ENDDO
!
      SolAbs(iw) = RsDnS                    
      Rnet(iw)   = RnetS                     
      SH(iw)     = SHs                      
      FTsoil0(iw) = fhG                      
      DO i = 1, m
         STsoil(iw,i) = T_sj(i)             
         CVsoil(iw,i) = C_vj(i)             
         CWsoil(iw,i) = C_wj(i)              
         ECO2soil(iw,i) = ECO2_mst(i)      
         Vsoil(iw, i) = vdelta_mst(i)
		 Psoil(iw,i) = P_mst(i)
		 dTsoil(iw,i) = T_mst(i)
		 dVsoil(iw,i) = V_mst(i)
		 Cp_soil(iw, i) = Cp_sj1_mst(i)
		 decrease_T(iw,i) = Lsblim*ECO2_mst(i)/Cp_sj1_mst(i)
			IF (C_vj(i) .GT. 1d20) THEN
				PRINT *, 'C_vj exceeds 1e20, aborting calculation'
				STOP
			END IF
      ENDDO
      FCO2soil(iw) = fCO2_mst               
	  QCO2soil(iw) = QCO2_mst
      
      Trad(iw) = Tr
      Rair(iw) = 1/gas                      

      IF ( nn .eq.nyear) THEN 
         write(50,'(i8, 6f12.4)') iw, SolAbs(iw),Rnet(iw),SH(iw),FTsoil0(iw),Trad(iw),Rair(iw)
         write(51,'(i8, 21e12.4)') iw, (STsoil(iw,i), i=1,m)
         write(52,'(i8, 21e12.4)') iw, (CVsoil(iw,i), i=1,m)
         write(53,'(i8, 21e12.4)') iw, (CWsoil(iw,i), i=1,m)
         write(54,'(i8, 22e12.4)') iw, FCO2soil(iw), (ECO2soil(iw,i),i=1,m)
         write(55,'(i8, 21e12.4)') iw, (Vsoil(iw,i), i=1,m)
		 write(56,'(i8, 21e12.4)') iw, (Psoil(iw,i), i=1,m)
		 write(80,'(i8, 21e12.4)') iw, (dTsoil(iw,i), i=1,m)
		 write(81,'(i8, 21e12.4)') iw, (dVsoil(iw,i), i=1,m)
         write(120,'(i8, 5e12.4)') iw, QCO2soil(iw), QCO2(iw), QCO2soil(iw)-QCO2(iw)
      ENDIF
      
200 ENDDO

    DO i = 1, m
      STsoil_ave(i) = fave(STsoil(:,i),ilst)
      CVsoil_ave(i) = fave(CVsoil(:,i),ilst)
      CWsoil_ave(i) = fave(CWsoil(:,i),ilst)
    ENDDO

100 ENDDO

301 format (9(x,f12.6))

 STOP 'mars_alsis: normal execution and exit'

END Program mars_alsis



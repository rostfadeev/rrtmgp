SUBROUTINE SLM_RAD_DRIVER (KST, KEND, KLEV, KPOSI,&
 & PAPHI, PAPRS , PAPRSF , PCCO2 , PNEB,&
 & PQO3  , PDELP   , PEMIS  , PALBD,PALBP,PMU0,PITM,PIVEG, &
 & PDAER, PQ    , PQSAT   , PQICE , PQLI , PT    , PTS,&
 & PFRSO, FLX, PFRTH,PFRSODS,PFRTHDS,PFLSODS,P_STOCH_AEROS,P_STOCH_RADI,&
! workspace
 & ZTAUCLD, &
 & SWU_bot,TRT_top, TRU_top, TRD_top, TRU_bot, TRUC_bot, TRDC_bot, PALB,ZOVLP_test)  

!**** interface to radiation routines: CLIRAD(SW) and RRTM (LW)

!     PURPOSE.
!     --------
!           INTERFACE TO radiation routines

!**   INTERFACE.
!     ----------

!     EXPLICIT ARGUMENTS :
!        --------------------
! KST    : START INDEX OF DATA IN KLON-LONG VECTOR
! KEND   : END   INDEX OF DATA IN KLON-LONG VECTOR
! KLON : VECTOR LENGTH
! KLEV   : NUMBER OF LEVELS
! KPOSI - latitude
! PAER   : (KLON,6,KLEV )     ; OPTICAL THICKNESS OF THE AEROSOLS
! PALBD  : (KLON)           ; DIFFUSE ALBEDO 
! PALBP  : (KLON)           ; PARALLEL ALBEDO 
! PAPHIF (KLON,KLEV) : geopotential at full levels
! PAPRS  : (KLON,KLEV+1)      ; HALF LEVEL PRESSURE
! PAPRSF : (KLON,KLEV )       ; FULL LEVEL PRESSURE
! PCCO2  :                      ; CONCENTRATION IN CO2 (PA/PA)
! PNEB  : (KLON,KLEV )       ; CLOUD FRACTIONAL COVER
! PQO3   : (KLON,KLEV )       ; OZONE MIXING RATIO (MASS)
! PDELP    : (KLON,KLEV)        ; LAYER PRESSURE THICKNESS
! PEMIS  : (KLON)             ; SURFACE EMISSIVITY
! PMU0   : (KLON)             ; cos(SOLAR ANGLE)
! PQ     : (KLON,KLEV )       ; SPECIFIC HUMIDITY PA/PA
! PQSAT    : (KLON,KLEV )       ; SATURATION SPECIFIC HUMIDITY PA/PA
! PQICE  : (KLON,KLEV )       ; ICE    WATER KG/KG
! PQLI  : (KLON,KLEV )       ; LIQUID WATER KG/KG
! PITM   : (KLON)             ; LAND-SEA MASK
! PIVEG(KLON) : index of dominant surface type
! PT     : (KLON,KLEV)        ; FULL LEVEL TEMPERATURE
! PTS    : (KLON)             ; SURFACE TEMPERATURE

!     ==== OUTPUTS ===
! PFRSO(KLON,KLEV+1) : total sky SW flux
! PFRTHC (KLON,2)              ; CLEAR-SKY LONGWAVE FLUXES
! PFRTH (KLON,KLEV+1)         ; TOTAL-SKY LONGWAVE FLUXES


!  M.Tolstykh after RECMWF routine form ARPEGE/IFS code

!-----------------------------------------------------------------------



USE MODCTR, ONLY : LSTOCH_PAR

USE YOMPHY2  , ONLY : NTSHM,NTSML       
USE YOMCST   , ONLY : RG       ,RD       ,RTT      ,RPI, RMD      ,RMO3
USE YOMPHY3   , ONLY : RII0, RRAE   ,REPH2O, & 
 &  RLINLI,RFUETA, RCAEROS,RCSTBGA,RAER,RCARDI, RLWUH,&
 & REPLOG   ,REPSC    ,REPSCW   , CVDAES   ,CVDAEL   ,CVDAEU   ,CVDAED  
USE MODGEM, ONLY : PGEMU, PGSQM2
USE RRTMG_LW_RAD
USE MOD_AEROSOLS  , ONLY :       REPAER   ,  &
 & RTAEBC  ,RTAEOR   ,RTAESD   ,RTAESS   ,RTAESU, &
 & LAER_MACv2, RTAEZ, macnbands, MOD_AEROSOLS_get_vdfrac, RTAEMAC_aod, &
 & get_ssa, get_asy, RTAEMAC_ssa, RTAEMAC_asy, get_ssa2, get_asy2
USE MODPHY, ONLY : NSTOCHPAR, M_AEROS, R_STOCH_PAR
use clirad, only : sorad, clirad_wp => rb

Use moving_average, Only: movave
Use NETCDF_io_out,  Only: diag_eb
use modprmt, only: nlon


!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=4),INTENT(IN) :: KLEV, KST, KEND, KPOSI

REAL(KIND=8), INTENT(IN) :: PALB(KST:KEND)
REAL(KIND=8), INTENT(IN) :: PALBD(KST:KEND) 
REAL(KIND=8), INTENT(IN) :: PALBP(KST:KEND) 
REAL(KIND=8), INTENT(IN) :: PAPHI(KST:KEND,0:KLEV) 
REAL(KIND=8), INTENT(IN) :: PAPRS(KST:KEND,KLEV+1)
REAL(KIND=8), INTENT(IN) :: PAPRSF(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PCCO2
REAL(KIND=8), INTENT(IN) :: PNEB(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PQO3(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PDELP(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PEMIS(KST:KEND),PDAER(KLEV) 
REAL(KIND=8), INTENT(IN) :: PMU0(KST:KEND) 
REAL(KIND=8), INTENT(IN) :: PQ(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PQSAT(KST:KEND,KLEV)
REAL(KIND=8), INTENT(IN) :: PQICE(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PQLI(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PITM(KST:KEND)
REAL(KIND=8), INTENT(IN) :: PIVEG(KST:KEND) 
REAL(KIND=8), INTENT(IN) :: PT(KST:KEND,KLEV) 
REAL(KIND=8), INTENT(IN) :: PTS(KST:KEND)
REAL(KIND=8), INTENT(IN) :: P_STOCH_AEROS(KST:KEND)
REAL(KIND=8), INTENT(IN) :: P_STOCH_RADI(KST:KEND)
Real(kind=4), intent(inout) :: ZOVLP_test(KST:KEND,KLEV)

REAL(KIND=8)                :: PFRSO(KST:KEND,KLEV+1) 
REAL(KIND=4), INTENT(OUT)   :: PFRSODS(KST:KEND),PFLSODS(KST:KEND)
REAL(KIND=4), INTENT(OUT)   :: PFRTHDS(KST:KEND)
REAL(KIND=8), INTENT(INOUT) :: PFRTH(KST:KEND,KLEV+1)

REAL(KIND=4), INTENT(OUT) :: SWU_bot(KST:KEND)
REAL(KIND=4), INTENT(OUT) :: TRT_top(KST:KEND)
REAL(KIND=4), INTENT(OUT) :: TRU_top(KST:KEND)
REAL(KIND=4), INTENT(OUT) :: TRD_top(KST:KEND)
REAL(KIND=4), INTENT(OUT) :: TRU_bot(KST:KEND)
REAL(KIND=4), INTENT(OUT) :: TRUC_bot(KST:KEND)
REAL(KIND=4), INTENT(OUT) :: TRDC_bot(KST:KEND)
! 
!     ------------------------------------------------------------------
!*       LOCAL ARRAYS.
!              -------------
REAL(KIND=8) :: ZTH(KST:KEND,KLEV+1)
REAL(KIND=8) :: ZRMU0(kst:kend)

REAL(KIND=8) :: ZQS    
REAL(KIND=8) :: ZTAUCLD(KST:KEND,KLEV,16)  
REAL(KIND=8) :: ZCLDSW(KST:KEND,KLEV)

Real(Kind=8) :: ZTAUAER(KST:KEND,5,KLEV)
Real(Kind=8) :: Ztauaer_MacV2(KST:KEND,macnbands,KLEV)
Real(Kind=8) :: Zssaal_MacV2 (KST:KEND,macnbands,KLEV)
Real(Kind=8) :: Zasyal_MacV2 (KST:KEND,macnbands,KLEV)

REAL(KIND=8), DIMENSION(KST:KEND,KLEV+1) :: UFLX, DFLX
REAL(KIND=4), DIMENSION(KST:KEND,KLEV+1) :: UFLXC,DFLXC
REAL(KIND=8) :: ZFIWP(KST:KEND,KLEV), ZFLWP(KST:KEND,KLEV)
REAL(KIND=8) :: ZIWC(KST:KEND), ZLWC(KST:KEND)

INTEGER(KIND=4) :: IFDIA(2),IIDIA(2)
INTEGER(KIND=4) :: IAUCR, I, K, IB, JLON
LOGICAL :: LLFDIA, LLIDIA!
REAL(KIND=8) :: ZDESR(KST:KEND), ZRADIP(KST:KEND,KLEV), ZRADLP(KST:KEND,KLEV)       

INTEGER(KIND=4) :: IKL, JK, JKL, JKLP1, JL, JRTM, Isb, Islt, Ict, ICb, ICLD

REAL(KIND=8) ::   ZIWGKG, ZLWGKG,&
 & ZMSAID, ZMSAIU, ZMSALD, ZLWFUDG,  ZMULTL, Z1RADL, Z1RADI, &
 &  ZRWGKG,  ZTEMPC,  ZDPOG, ZPODT,ZAIWC,ZBIWC  
REAL(KIND=8) :: ZALND, ZASEA, ZD, ZDEN, ZNTOT, ZNUM
REAL (8), SAVE :: TAU_S(8),SSA_S(8),ASYM_S(8)
REAL (8), SAVE :: TAU_M(8),SSA_M(8),ASYM_M(8)       
REAL (8), SAVE :: TAU_C(8),SSA_C(8),ASYM_C(8)       
REAL (8) TAUA_V(KST:KEND,KLEV)
real(8) flx(KST:KEND,KLEV+1)
real(4) flx_d(kst:kend,KLEV+1)
logical lovercast,lcldwater,lcrel
real(8) ZCO2MWG,ZAIRMWG
REAL (8) ZAER(KLEV,6)

Integer(Kind=4) nrtae(5),n
Real(Kind=4) zrtae

REAL(KIND=8) :: ZCRAE, ZRII0, ZTCELS,ZFSR,  ZEMIW(KST:KEND)
REAL(kind=8) :: ZQ(KST:KEND,KLEV)

real(kind=clirad_wp) :: Ztaual(kst:kend,KLEV,8),Zssaal(kst:kend,KLEV,8), Zasyal(kst:kend,KLEV,8)
real(kind=clirad_wp) :: Zcwc(kst:kend,KLEV,3),Ztaucldsw(kst:kend,KLEV,3),Zreff(kst:kend,KLEV,3)
REAL(kind=clirad_wp) :: ZPL(kst:kend,KLEV+1)
real(kind=clirad_wp) :: Zco2
REAL(kind=clirad_wp) :: PQO3_clirad_wp(kst:kend,KLEV)
real(kind=clirad_wp) :: PALBP_clirad_wp(kst:kend)
real(kind=clirad_wp) :: PALBD_clirad_wp(kst:kend)
REAL(kind=clirad_wp) :: ZT(kst:kend,KLEV), ZNEB(kst:kend,KLEV)
REAL(kind=clirad_wp) :: ZQ_clirad_wp(kst:kend,KLEV), ZDEL(kst:kend,KLEV)

      Logical :: ifexist = .True.
      Real(Kind=4), Allocatable :: r4surf(:), r43d(:,:), r8zn(:), r4col(:)
      Real(Kind=8), Allocatable :: r82(:)
      Character(4) c4(6)
      Character(2) c2
      Integer(Kind=4) ddr(4,3), ldr(2,3)
      Integer(Kind=4) ddrL(4,3), ldrL(2,3)
      Integer(Kind=4) ddrL3d(4,3), ldrL3d(2,3)
      Integer(Kind=4) ddrL3dp(4,3), ldrL3dp(2,3)
      real(kind=8) CO2(KST:KEND,KLEV)
      real(kind=8) CH4(KST:KEND,KLEV)
      real(kind=8) N2O(KST:KEND,KLEV)
      real(kind=8) CFC11(KST:KEND,KLEV)
      real(kind=8) CFC12(KST:KEND,KLEV)
      real(kind=8) o2_conc
      real(kind=4) longitude(KST:KEND)
      real(kind=4) latitude(KST:KEND)

!                CONTINENTAL (C) AEROSOLS,   8 intervals:
     
    TAU_C=(/2.09880,1.76676,1.08667,0.87224,0.56770,0.25307, 0.28098,0.12890/)
     
     SSA_C=(/0.77936,0.89032,0.89711,0.86887,0.83504,0.75306,0.76509,0.69953/)
	
      ASYM_C=(/0.68684,0.65619,0.63822,0.63599,0.63209,0.68143,0.65978,0.77775/)
!
!
!              MARITIME (M)  AEROSOLS,   8 intervals
!
       TAU_M=(/1.23760,1.16595,1.01222,0.97019,0.90512,0.76353,0.80001,0.60126/)  
!
       SSA_M=(/0.92032,0.97951,0.99000,0.98782,0.98620,0.96099,0.98319,0.86221/)
!
       ASYM_M=(/0.75748,0.74395,0.74600,0.74975,0.75541,0.77743,0.77184,0.80229/)

!
!              STRATOSPHERIC (S)  AEROSOLS,  8 intervals
!
       TAU_S=(/1.55000,1.53946,1.09667,0.81293,0.45271,0.10020,0.11383,0.03955/) 
!
       SSA_S=(/0.9999,0.9999,0.9999,0.9999,0.9999,0.90023,0.99161,0.49372/)
!        
       ASYM_S=(/0.69380,0.71722,0.73244,0.68702,0.63582,0.39456,0.44335,0.17748/)

!*       1.1    LOCAL CONSTANTS
!                ---------------
 !    'option for cloud optical thickness scaling' 
!     'true if scaling is not required; false: if scaling is required'
     lovercast=.FALSE.
!'If lcldwater=.true., taucld is computed from cwc and reff.'  
!'If lcldwater=.false., taucld is an input parameter'
     lcldwater=.true.
!	'IF calculate without molecular scattering - false '
     lcrel=.TRUE.
     Isb=999
     Islt=999
     ict=NTSHM
     icb=NTSML
      ICLD=2
      ZAIRMWG = 28.970_8
      ZCO2MWG = 44.011_8
     Zco2=RCARDI/(1.E-6*ZCO2MWG/ZAIRMWG)
! 401. !323.86 !383.0 !co2 concentration ppmv
      ZAIRMWG = 28.970_8
      ZCO2MWG = 44.011_8

!      RCARDI  = Zco2*1.E-06_8*ZCO2MWG/ZAIRMWG

! s0=solar cons
      DO k= 1, KLEV
         DO i=KST, KEND
! if cldwater=.true. taucld is computed further
            ztaucld(i,k,1)=0.!  ZDEL0I(i,k) !ICE PARTICLES 
            ztaucld(i,k,2)=0.!ZDEL0L(i,k)   ! WATER PARTICLES 
            ztaucld(i,k,3)= 0.0  ! RAIN DROPS 
            Zreff(i,k,1)  =40.!REFI  ! ICE PARTICLES
            Zreff(i,k,2)  =5.25 !REFW  ! WATER PARTICLES
 !           fcld(i,k)    = PNEB(i,k)  ! CLOUD COVER 
            Zcwc(i,k,1)   = PQICE(i,k)!/PNEB(i,k)
            Zcwc(i,k,2)   =PQLI(i,k)!/PNEB(i,k)
            Zcwc(i,k,3)   = 0.0
         ENDDO
       ENDDO
   ZT(KST:KEND,:)=PT(KST:KEND,:)
   ZNEB(KST:KEND,:)=max(1.e-6,PNEB(KST:KEND,:))
   !height in km from geopotential height in m
!calculate layer depth, km
!     zl(KST:KEND,1:KLEV+1)=PAPHI(KST:KEND,0:KLEV)/(9.81*1000.0)

   DO k=1,KLEV
     do I=KST, KEND
       zdel(i,k)=(PAPHI(I,k-1)-PAPHI(I,k))/9.80665E3_8
      end do
     end do
     Zpl(KST:KEND,1:KLEV+1)=PAPRS(KST:KEND,1:KLEV+1)*0.01
!           SOLAR ZENITH ANGLE IS EARTH'S CURVATURE CORRECTED
      ZRMU0(KST:KEND)=PMU0(KST:KEND)
! places to calculate solar radiation (from APLPAR)
      IAUCR=0
      LLIDIA=ZRMU0(KST) > 0._8
      LLFDIA=.FALSE.

      IF (LLIDIA) THEN
      IAUCR=1
      IIDIA(1)=KST
      ENDIF
  
  !     NON VECTORIZED LOOP.
  
      DO JLON=KST,KEND
      LLFDIA=ZRMU0(JLON) > 0._8
  
      IF (LLFDIA.NEQV.LLIDIA) THEN
  
      IF (LLFDIA) THEN
        IAUCR=IAUCR+1
        IIDIA(IAUCR)=JLON
      ELSE
        IFDIA(IAUCR)=JLON-1
      ENDIF
  
      LLIDIA=LLFDIA
      ENDIF
  
      ENDDO
  
      IF (LLFDIA) THEN
      IFDIA(IAUCR)=KEND
      ENDIF

! inhomogenity correction
!  IF (KLEV == 51) THEN  
    ZLWFUDG=RLWUH !0.67_8
! preparation of arrays
!   IF (KPOSI == 111 .AND. KST ==2) write(*,*)  RLWUH
DO JK=1,KLEV
  DO JL=KST, KEND
    ZQS =MAX(2.0_8*REPH2O,PQSAT(JL,JK))
    ZQ  (JL,JK)=MAX(REPH2O,MIN(PQ(JL,JK),ZQS*(1.0_8-REPH2O)))
    ZQ_clirad_wp(JL,JK) = ZQ(JL,JK)
    ZEMIW(JL)= 0.99_8
  ENDDO
ENDDO

!  interpolation to half-levels

DO JK=2,KLEV
  DO JL=KST, KEND
    ZTH(JL,JK)= &!MAX(168.,MIN(355., &
     & (PT(JL,JK-1)*PAPRSF(JL,JK-1)*(PAPRSF(JL,JK)-PAPRS(JL,JK))&
     & +PT(JL,JK)*PAPRSF(JL,JK)*(PAPRS(JL,JK)-PAPRSF(JL,JK-1)))&
     & *(1.0_8/(PAPRS(JL,JK)*(PAPRSF(JL,JK)-PAPRSF(JL,JK-1))))!))  
  ENDDO
ENDDO

! 
UFLX=0.
DFLX=0.
UFLXC=0.
DFLXC=0.
! vertical boundaries
DO JL= KST, KEND
  ZTH(JL,KLEV+1)=PTS(JL)
  ZTH(JL,1)=PT(JL,1)-PAPRSF(JL,1)*(PT(JL,1)-ZTH(JL,2))&
   & /(PAPRSF(JL,1)-PAPRS(JL,2))  
ENDDO
!------------ FORMATION OF AEROSOL VERTICAL PROFILE IN STRATOSPHERE (T.Tarasova)
if (LAER_MACv2) then
        
    DO I = KST, KEND
        nrtae(:)=0
        do JK = 1, KLEV
            if (0._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) <  1._4) nrtae(1) = nrtae(1)+1
            if (1._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) <  3._4) nrtae(2) = nrtae(2)+1
            if (3._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) <  6._4) nrtae(3) = nrtae(3)+1
            if (6._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) < 12._4) nrtae(4) = nrtae(4)+1
        end do
          
        DO JK = 1, KLEV
            zrtae = 0._4
            if (0._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) <  1._4) zrtae = 1._4/Real(nrtae(1),4)
            if (1._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) <  3._4) zrtae = 1._4/Real(nrtae(2),4)
            if (3._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) <  6._4) zrtae = 1._4/Real(nrtae(3),4)
            if (6._4 <= RTAEZ(I,JK,KPOSI) .and. RTAEZ(I,JK,KPOSI) < 12._4) zrtae = 1._4/Real(nrtae(4),4)
            
            Zssaal_MacV2 (I,1:macnbands,JK) = SSA_S(1)
            Zasyal_MacV2 (I,1:macnbands,JK) = ASYM_S(1)
            
            Ztauaer_MacV2(I,1:macnbands,JK) = zrtae*MOD_AEROSOLS_get_vdfrac(RTAEZ(I,JK,KPOSI))*RTAEMAC_aod(I,KPOSI,1:macnbands)
            !Zssaal_MacV2 (I,1:macnbands,JK) = zrtae*MOD_AEROSOLS_get_vdfrac(RTAEZ(I,JK,KPOSI))*RTAEMAC_ssa(I,KPOSI,1:macnbands)
            !Zasyal_MacV2 (I,1:macnbands,JK) = zrtae*MOD_AEROSOLS_get_vdfrac(RTAEZ(I,JK,KPOSI))*RTAEMAC_asy(I,KPOSI,1:macnbands)
            
            do n = 1,macnbands
                if (Ztauaer_MacV2(I,n,JK) < RCAEROS) Ztauaer_MacV2(I,n,JK) = RCAEROS
                !if (Zssaal_MacV2 (I,n,JK) < RCAEROS) Zssaal_MacV2 (I,n,JK) = 1._8
                !if (Zasyal_MacV2 (I,n,JK) < RCAEROS) Zasyal_MacV2 (I,n,JK) = 1._8
            end do
        end do
        !Ztaual(KLON,KLEV,8)
        ! soluv nband=3:  0.2-0.303, 0.303-0.323, 0.323-0.7  mcm
        ! solir
        !1                       0.323-1.22
        !2     14280-8200        0.70-1.22
        !3     8200-1000         1.22-10.0
        !4     8200-4400         1.22-2.27
        !5     4400-1000         2.27-10.0
        !                               1    2    3    4    5     6     7     8     9    10    11    12    13    14
        !RTAEMAC_nm(1:macnbands) = (/ 230, 300, 400, 550, 700, 1000, 1270, 1460, 1780, 2050, 2320, 2790, 3470, 8000 /)
          
        DO JK = 1, KLEV
            Ztaual(I,JK,1) = & ! 33000-50000 cm-1 200-  303 nm soluv
                             & Ztauaer_MacV2(I, 1,JK)
             Ztaual(I,JK,2) = & ! 31000-33000      303-  323    soluv
                              & Ztauaer_MacV2(I, 2,JK)
             Ztaual(I,JK,3) = & ! 14280-31000      323-  700    soluv
                              & Ztauaer_MacV2(I, 3,JK)+Ztauaer_MacV2(I, 4,JK)+Ztauaer_MacV2(I, 5,JK)
             Ztaual(I,JK,4) = & !  8200-31000      323- 1220    solir
                              & Ztauaer_MacV2(I, 3,JK)+Ztauaer_MacV2(I, 4,JK)+Ztauaer_MacV2(I, 5,JK)+Ztauaer_MacV2(I, 6,JK)
             Ztaual(I,JK,5) = & !  8200-14280      700- 1220    solir
                              & Ztauaer_MacV2(I, 5,JK)+Ztauaer_MacV2(I, 6,JK)
             Ztaual(I,JK,6) = & !  1000- 8200     1220-10000    solir
                              & Ztauaer_MacV2(I, 7,JK)+Ztauaer_MacV2(I, 8,JK)+Ztauaer_MacV2(I, 9,JK)+Ztauaer_MacV2(I,10,JK) &
                              &+Ztauaer_MacV2(I,11,JK)+Ztauaer_MacV2(I,12,JK)+Ztauaer_MacV2(I,13,JK)+Ztauaer_MacV2(I,14,JK)
             Ztaual(I,JK,7) = &       !  4400- 8200     1220- 2270    solir
                              & Ztauaer_MacV2(I, 7,JK)+Ztauaer_MacV2(I, 8,JK)+Ztauaer_MacV2(I, 9,JK)+Ztauaer_MacV2(I,10,JK)
             Ztaual(I,JK,8) = &       !  1000- 4400     2270-10000    solir
                              & Ztauaer_MacV2(I,11,JK)+Ztauaer_MacV2(I,12,JK)+Ztauaer_MacV2(I,13,JK)+Ztauaer_MacV2(I,14,JK)
             !Ztaual(I,JK,1:8) = 0.9_8*Ztaual(I,JK,1:8)
        end do
    end do ! i

    DO JK = 1, KLEV        
        DO I = KST, KEND
            Zssaal(i,JK,1:8) = SSA_S (1:8)
            Zasyal(i,JK,1:8) = ASYM_S(1:8)
        end do
    end do
        
else ! LAER_MACv2 == False
    DO I=KST, KEND
        ZAER=0.
          
        Isb=KLEV
        Islt=KLEV
        DO k=1, KLEV
            ZAER(K,6) = RCSTBGA
            IF (PAPHI(i,k-1) > 12.0*9.81E3.AND.PAPHI(i,k-1) <=20.0*9.81E3) THEN
                IF (LSTOCH_PAR) THEN
                    TAUA_V(i,k)=P_STOCH_AEROS(I)*zdel(i,k)
                ELSE
                    TAUA_V(i,k)=0.000218*zdel(i,k)
                ENDIF
                Isb=k
            ELSE IF (PAPHI(i,k-1) > 20.0*9.81E3) THEN
                TAUA_V(i,k)=0.
            END IF
        END DO

        IF (LSTOCH_PAR) THEN
            do ib=1,8
                do k= 1, Isb
!                    Ztaual(i,k,ib) = TAU_S(ib)*TAUA_V(i,k)*(0.05+PGSQM2(KPOSI))
                    Ztaual(i,k,ib) = TAU_S(ib)*PDAER(k)* P_STOCH_AEROS(I)!*(0.05+PGSQM2(KPOSI))        
                    Zssaal(i,k,ib) = SSA_S(ib)
                    Zasyal(i,k,ib) = ASYM_S(ib)
                end do
            end do
        ELSE
            do ib = 1,8
                do k = 1, Isb
!                   Ztaual(i,k,ib) = TAU_S(ib)*TAUA_V(i,k)*(0.05+PGSQM2(KPOSI))
                    Ztaual(i,k,ib) = TAU_S(ib)*PDAER(k)!*(0.05+PGSQM2(KPOSI))
                    Zssaal(i,k,ib) = SSA_S(ib)
                    Zasyal(i,k,ib) = ASYM_S(ib)
                end do
            end do
        END IF
        do k = 1, Isb
!           PAER(I,6,K) = TAU_S(8)*TAUA_V(i,k)*(0.05+PGSQM2(KPOSI))
            ZAER(K,6) = TAU_S(8)*PDAER(k)!*(0.05+PGSQM2(KPOSI))
        end do
!---------IF MAR-I aerosol profile in troposphere
!
        IF (PITM(I) < 0.5_8 .OR. NINT(PIVEG(I)) == 2) THEN
            DO k = Isb+1, KLEV        
                IF (PAPHI(I,k-1) > 0._8.AND.PAPHI(I,k-1) <=2.0*9.81E3) THEN
                    TAUA_V(I,k)=0.025_8*zdel(I,k)
                END IF
                IF (PAPHI(I,k-1) >2.*9.81E3.AND.PAPHI(I,k-1) <= 12.0*9.81E3) THEN
                    IF (PAPHI(I,k-1) > 2.*9.81E3.AND.PAPHI(I,k-1) <= 4.*9.81E3) THEN
                        TAUA_V(I,k)=0.0025_8*zdel(I,k)
                    ELSE IF (PAPHI(I,k-1) > 4.*9.81E3.AND.PAPHI(I,k-1) <= 8.*9.81E3) THEN
                        TAUA_V(I,k)=0.0001_8*zdel(I,k)
                    ELSE
                        TAUA_V(I,k)=0.0!*zdel(I,k)
                    END IF 
                    Islt=k
                END IF
            END DO
            do ib=1, 8
                do k= Isb+1, Islt
!                   Ztaual(i,k,ib) = TAU_C(ib)*TAUA_V(i,k)
                    Ztaual(i,k,ib) = TAU_C(ib)*PDAER(k)*(RTAESD(i,KPOSI))/TAU_C(8)
                    Zssaal(i,k,ib) = SSA_C(ib)
                    Zasyal(i,k,ib) = ASYM_C(ib)
                end do
            end do
            do k= Isb+1, Islt
!               ZAER(K,3) = TAU_C(8)*TAUA_V(i,k)
                ZAER(K,3) = (RTAESD(i,KPOSI))*PDAER(K)
            end do
            do ib=1,8
                do k= Islt+1, KLEV
!                   Ztaual(i,k,ib) = TAUA_V(i,k)*TAU_M(ib)
                    Ztaual(i,k,ib) = TAUA_V(i,k)*(TAU_M(ib)*(RTAESS(i,KPOSI))/TAU_M(8))
                    Zssaal(i,k,ib) = SSA_M(ib)
                    Zasyal(i,k,ib) = ASYM_M(ib)
                end do
            end do
            do k= Islt+1, KLEV
!               ZAER(K,2) = TAU_M(8)*TAUA_V(i,k)! PDAER(K)
                ZAER(K,2) = (RTAESS(i,KPOSI))*TAUA_V(i,k)
            end do
        ELSE
!---------IF CONT-I aerosol profile in troposphere        
!         
            DO k=Isb+1, KLEV
                IF (PAPHI(I,k-1) > 0.0.AND.PAPHI(I,k-1) <= 2.0*9.81E3) THEN
                    TAUA_V(I,k)=0.1_8*zdel(I,k)
                END IF
                IF (PAPHI(i,K-1) > 2.*9.81E3.AND.PAPHI(I,k-1) <= 12.0*9.81E3) THEN
                    TAUA_V(I,k)=0.0025_8*zdel(I,k)
                END IF
            END DO

            do ib=1,8
                do K= Isb+1, KLEV
!                   Ztaual(i,k,ib) = TAU_C(ib)*TAUA_V(i,k)
                    Ztaual(i,k,ib) = TAU_C(ib)*TAUA_V(i,k)*(RTAESS(i,KPOSI)+RTAESD(i,KPOSI))/TAU_C(8)
                    Zssaal(i,k,ib) = SSA_C(ib)
                    Zasyal(i,k,ib) = ASYM_C(ib)
                end do
            end do           
            do K= Isb+1, KLEV
!               ZAER(K,3) = TAU_C(8)*TAUA_V(i,k)
                ZAER(K,3) = (RTAESS(i,KPOSI)+RTAESD(i,KPOSI))*TAUA_V(i,k)
            end do
!      
        END IF ! (PITM(I) < 0.5_8 .OR. NINT(PIVEG(I)) == 2)
! aerosol optical depth is converted for RRTM bands
        DO JK = 1, KLEV
            ZAER(JK,1)=RCAEROS
            ZAER(JK,4) =RCAEROS
            ZAER(JK,5)=RCAEROS
            ZTAUAER(I,1,JK) = RAER(1,1)*ZAER(JK,1)+RAER(1,2)*ZAER(JK,2)&
                            &+RAER(1,3)*ZAER(JK,3)+RAER(1,4)*ZAER(JK,4)&
                            &+RAER(1,5)*ZAER(JK,5)+RAER(1,6)*ZAER(JK,6)  
            ZTAUAER(I,2,JK) = RAER(2,1)*ZAER(JK,1)+RAER(2,2)*ZAER(JK,2)&
                            &+RAER(2,3)*ZAER(JK,3)+RAER(2,4)*ZAER(JK,4)&
                            &+RAER(2,5)*ZAER(JK,5)+RAER(2,6)*ZAER(JK,6)  
            ZTAUAER(I,3,JK) = RAER(3,1)*ZAER(JK,1)+RAER(3,2)*ZAER(JK,2)&
                            &+RAER(3,3)*ZAER(JK,3)+RAER(3,4)*ZAER(JK,4)&
                            &+RAER(3,5)*ZAER(JK,5)+RAER(3,6)*ZAER(JK,6)  
            ZTAUAER(I,4,JK) = RAER(4,1)*ZAER(JK,1)+RAER(4,2)*ZAER(JK,2)&
                            &+RAER(4,3)*ZAER(JK,3)+RAER(4,4)*ZAER(JK,4)&
                            &+RAER(4,5)*ZAER(JK,5)+RAER(4,6)*ZAER(JK,6)  
            ZTAUAER(I,5,JK) = RAER(5,1)*ZAER(JK,1)+RAER(5,2)*ZAER(JK,2)&
                            &+RAER(5,3)*ZAER(JK,3)+RAER(5,4)*ZAER(JK,4)&
                            &+RAER(5,5)*ZAER(JK,5)+RAER(5,6)*ZAER(JK,6)
        END DO
    END DO !I
    
end if ! ! LAER_MACv2
    
! Ozone at half levels (index will be inverted in RRTMG-LW, s/r inatm)
!  DO JK = 1, KLEV
!    DO JL = KST, KEND
!      ZZOZ(JL,JK)=0.5_8*(PQO3(JL,JK-1)+PQO3(JL,JK))
!    ENDDO
!  ENDDO

!Invert cloudiness for RRTM
DO JK = 1 , KLEV
    JKL = KLEV+ 1 - JK
    DO JL = KST,KEND
        ZCLDSW(JL,JK)  = MAX( 0.0_8 ,MIN( 1.0_8 ,PNEB(JL,JKL)))
    END DO
END DO
!   The same for LWP, IWP, eff. ice and droplet radius  and optical depths
DO JK = 1 , KLEV
    IKL = KLEV + 1 - JK
    DO JL = KST,KEND
! --- LiqWaterContent (g.m-3) AND LiqWaterPath (g.m-2)
        IF (PNEB(JL,IKL) > REPSC ) THEN
            ZLWGKG=MAX(PQLI(JL,IKL)*1000._8,0.0_8)
            ZIWGKG=MAX(PQICE(JL,IKL)*916._8,0.0_8)
        ELSE
            ZLWGKG=0.0_8
            ZIWGKG=0.0_8
        END IF
        ZDPOG=PDELP(JL,IKL)/RG
        ZFLWP(JL,JK)= ZLWGKG*ZDPOG
        ZFIWP(JL,JK)= ZIWGKG*ZDPOG
        ZPODT=PAPRSF(JL,IKL)/(RD*PT(JL,IKL))
        ZLWC(JL)=ZLWGKG*ZPODT
        ZIWC(JL)=ZIWGKG*ZPODT
    END DO
    DO JL = KST,KEND
        IF (PITM(JL) < 0.5_8) THEN
            ZASEA=50._8
            ZD=0.33_8
            ZNTOT=-1.15E-03_8*ZASEA*ZASEA+0.963_8*ZASEA+5.30_8
        ELSE
            ZALND=900._8 !900._8
            ZD=0.43_8
            ZNTOT=-2.10E-04_8*ZALND*ZALND+0.568_8*ZALND-27.9_8
        END IF
        ZNUM=3._8*ZLWC(JL)*(1._8+3._8*ZD*ZD)**2
        ZDEN=4._8*RPI*ZNTOT*(1._8+ZD*ZD)**3
        IF ((ZNUM/ZDEN) > REPLOG) THEN
            IF (LSTOCH_PAR) THEN
                ZRADLP(JL,JK)=P_STOCH_RADI(JL)*100.*EXP(0.333*LOG(ZNUM/ZDEN))
            ELSE
                ZRADLP(JL,JK)=100.*EXP(0.333*LOG(ZNUM/ZDEN))
            END IF
            ZRADLP(JL,JK)=MAX(ZRADLP(JL,JK), 4._8) ! orig - 4
            !ZRADLP(JL,JK)=MAX(ZRADLP(JL,JK), 5._8) ! orig - 4
            !ZRADLP(JL,JK)=MIN(ZRADLP(JL,JK),16._8)! bestsofar25._8) !orig - 16
            ZRADLP(JL,JK)=MIN(ZRADLP(JL,JK),25._8)! bestsofar25._8) !orig - 16
            Zreff(JL,IKL,2)=ZRADLP(JL,JK)
            Zreff(JL,JK,3)=0.
        ELSE
            ZRADLP(JL,JK)=4._8
        END IF
        ! diagnosing the ice particle effective radius/diameter
        IF (ZIWC(JL) > 0.0001_8 ) THEN
            ZTEMPC = PT(JL,IKL)-83.15_8
            ZTCELS = PT(JL,IKL)-RTT
            ZFSR = 1.2351_8 +0.0105_8 * ZTCELS
            !Sun, 2001 (corrected from Sun & Rikus, 1999)
            ZAIWC = 45.8966_8 * ZIWC(JL)**0.2214_8
            ZBIWC = 0.7957_8 * ZIWC(JL)**0.2535_8
            IF (LSTOCH_PAR) THEN
                ZDESR(JL) = (ZFSR * (ZAIWC + ZBIWC*ZTEMPC))*P_STOCH_RADI(JL)
            ELSE
                ZDESR(JL) = ZFSR * (ZAIWC + ZBIWC*ZTEMPC)
            ENDIF
            !ZDESR(JL) = MIN ( MAX( ZDESR(JL), 85._8), 350._8) !85,350
            ZRADIP(JL,JK)= 3._8*sqrt(3._8)*0.125_8*MIN ( MAX( ZDESR(JL), 85._8), 199._8) !ZDESR(JL)
        ELSE
            ZDESR(JL) = 45._8
            ZRADIP(JL,JK)= 25._8
        END IF  
        Zreff(JL,IKL,1)=ZRADIP(JL,JK)     
    END DO !JL
!   optical cloud depth for RRTMG LW
    DO JRTM=1,16
        DO JL = KST,KEND
            ZTAUCLD(JL,JK,JRTM) = 0.0_8
            ZMSALD = 0.0_8
            ZMSAID = 0.0_8
            IF (ZFLWP(JL,JK)+ZFIWP(JL,JK) > REPSCW) THEN
                Z1RADL = 1.0_8 / ZRADLP(JL,JK)
                ZMSALD = RLINLI(JRTM,1)+ZRADLP(JL,JK)*RLINLI(JRTM,2)+ Z1RADL*&
                & (RLINLI(JRTM,3) + Z1RADL*(RLINLI(JRTM,4) + Z1RADL*&
                & RLINLI(JRTM,5) ))  
          
                Z1RADI = 1._8 / ZDESR(JL)
                ZMSAID = RFUETA(JRTM,1) + Z1RADI*(RFUETA(JRTM,2) + Z1RADI*RFUETA(JRTM,3))
                ZTAUCLD(JL,JK,JRTM) = ZLWFUDG*(ZMSALD*ZFLWP(JL,JK)+ZMSAID*ZFIWP(JL,JK))
            END IF
        END DO
    END DO
END DO ! JK
! if (maxval(ZRADIP)> 140._8) write(*,*) KPOSI,KST,maxval(ZRADIP),maxloc(ZRADIP)
! SW radiation - CALL CLIRAD for every piece
flx(KST:KEND,:) = 0.
flx_d(KST:KEND,:)=0.
SWU_bot(KST:KEND)=0.
       
do k = 1, KLEV 
    do i = KST, KEND
        PQO3_clirad_wp(i, k) = PQO3(i, k)
    end do
end do

do i = KST, KEND
    PALBP_clirad_wp(i) = PALBP(i)
end do

do i = KST, KEND
    PALBD_clirad_wp(i) = PALBD(i)
end do
       
IF (IAUCR >= 1) THEN
    DO i=1,IAUCR
        CALL SORAD (IIDIA(i),IFDIA(i),kst, kend, KLEV,lcrel, &
     &              Zpl,ZT,ZQ_clirad_wp,PQO3_clirad_wp,Zco2, &
     &              lovercast,lcldwater,Zcwc,Ztaucldsw,Zreff,ZNEB,ict,icb, &
     &              Ztaual,Zssaal,Zasyal,PITM(KST), &
     &              ZRMU0,PALBP_clirad_wp,PALBD_clirad_wp,&
     &              PALBP_clirad_wp,PALBD_clirad_wp, &
!     &             fdiruv,fdifuv,fdirpar,fdifpar,fdirir,fdifir,flx, flc,flx_d,SWU_bot)
                    flx, flx_d,SWU_bot)
    END DO
END IF 
! downward SW flux at the surface - full and intermediate values to store
DO I = KST, KEND
    PFRSODS(I) = flx_d(I,KLEV+1)*RII0*ZRMU0(i)
    PFLSODS(I) = flx_d(I,KLEV+1)
!   SWU_bot(i)=SWU_bot(i)*RII0*ZRMU0(i)
!   IF (ABS(SWU_bot(i)) > 1000.) write(*,*) KST,I,KPOSI,SWU_bot(i)           
 !! PFRSOPS(i)=real((fdiruv(i)+fdirpar(i)+fdirir(i)),8)*RII0*ZRMU0(i)
!   ZFRSOD(i)=(fdifuv(i)+fdifpar(i)+fdifir(i))*RII0*ZRMU0(i)
END DO
! Total SW flux
DO K = 1, KLEV+1
    DO I = KST, KEND
        PFRSO(i,k)=flx(i,k)*RII0*ZRMU0(i)
    END DO
END DO

Ztauaer = 0.1e-8

! LW radiation   -  CALL  RRTMG LW
if (LAER_MACv2) then

    Ztauaer_MacV2 = 0.1e-8

    CALL rrtmg_lw &
        (KST,KEND,KLEV, ICLD    ,0    , &
        PAPRSF    , PAPRS    ,PT    , ZTH    ,PTS    , &
        ZQ  ,PQO3   , PCCO2,& !co2vmr  ,  ch4vmr  ,n2ovmr  ,o2vmr, &
!       cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
!       PEMIS,Zemiw, 0 ,2,1, ZCLDSW   , &
!       PEMIS,Zemiw, 2 ,3,1, ZCLDSW   , &
        PEMIS,Zemiw, 2 ,2,1, ZCLDSW   , &
        Ztaucld  ,ZFIWP  ,ZFLWP  , ZRADIP   ,ZRADLP   , &
        tauaer_MacV2 = Ztauaer_MacV2  ,  &
        uflx = uflx, dflx = dflx, uflxc = uflxc, dflxc = dflxc)


      subroutine rrtmg_lw &
            (KST,KEND, nlay    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   , co2vmr  , & !  ch4vmr  ,n2ovmr  ,o2vmr, &
!             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
             emis,emiw, inflglw ,iceflglw,liqflglw,cldfr   , &
             taucld  ,cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer, tauaer_MacV2  , &
            uflx    ,dflx    ,uflxc   ,dflxc,  &
!             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
             duflx_dt,duflxc_dt )



        
else
    Ztauaer = 0.1e-8

    CALL rrtmg_lw &
        (KST,KEND,KLEV, ICLD    ,0    , &
        PAPRSF    , PAPRS    ,PT    , ZTH    ,PTS    , &
        ZQ  ,PQO3   , PCCO2,& !co2vmr  ,  ch4vmr  ,n2ovmr  ,o2vmr, &
!       cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
!       PEMIS,Zemiw, 0 ,2,1, ZCLDSW   , &
!       PEMIS,Zemiw, 2 ,3,1, ZCLDSW   , &
        PEMIS,Zemiw, 2 ,2,1, ZCLDSW   , &
        Ztaucld  ,ZFIWP  ,ZFLWP  , ZRADIP   ,ZRADLP   , &
        tauaer = Ztauaer, &
        uflx = uflx, dflx = dflx, uflxc = uflxc, dflxc = dflxc)
end if

! total LW flux; also invert index in vertical
DO JKL = 1 , KLEV+1
    JK = KLEV+1 + 1 - JKL
    DO JL = KST,KEND
        PFRTH(JL,JKL) = -(uflx(JL,JK)-dflx(JL,JK))
    END DO 
END DO

DO JL = KST,KEND
    TRT_top(JL) = -(uflx(JL,KLEV+1)-dflx(JL,KLEV+1))
    TRU_top(JL) = -uflx(JL,KLEV+1)
    TRD_top(JL) =  dflx(JL,KLEV+1)
    TRU_bot(JL) = -uflx(JL,1)
!   TRD_bot(JL) = -dflx(JL,1)
    TRUC_bot(JL) = -uflxc(JL,1)
    TRDC_bot(JL) = -dflxc(JL,1)
END DO 

! downward LW flux at the surface
DO JL=KST, KEND
    PFRTHDS(JL) =  dflx(JL,1) ! - before
END DO


        if (.True.) then
        if (movave%ready) then
          
          ifexist = .True.
          c4(1:6) = (/ "tst_", "sim_", "hr6_", "day_", "mon_", "yea_" /)
          
          ddrL = Reshape ((/ KST,KST,KEND,KEND, 1,1,1,1,  1,    1,    1,1 /), (/ 4,3 /))
          ldrL = Reshape ((/     KST,KEND,        1,1,      KPOSI,KPOSI   /), (/ 2,3 /))
          
          ddrL3d = Reshape ((/ KST,KST,KEND,KEND, 1,1,KLEV,KLEV,  1,    1,    1,1 /), (/ 4,3 /))
          ldrL3d = Reshape ((/     KST,KEND,        1,KLEV,           KPOSI,KPOSI   /), (/ 2,3 /))
          
          ddrL3dp = Reshape ((/ KST,KST,KEND,KEND, 1,1,KLEV+1,KLEV+1,  1,    1,    1,1 /), (/ 4,3 /))
          ldrL3dp = Reshape ((/     KST,KEND,        1,KLEV+1,           KPOSI,KPOSI   /), (/ 2,3 /))

          co2 = 0.0004_8
          n2o = 2.5e-7_8
          ch4 = 1.e-6_8
          cfc11 = 7.e-10_8
          cfc12 = cfc11
          o2_conc = 0.209488
          !longitude = 5._8
          !do i = KST, KEND
          !  latitude(i) = -180._4*(real(i,4)-1)/(real(nlon,4)-1) + 90
          !end do
          
          ZOVLP_test(KST:KEND,1:KLEV) = 0._4
          
          !REAL(KIND=8), DIMENSION(KST:KEND,KLEV+1) :: UFLX, DFLX
          !REAL(KIND=4), DIMENSION(KST:KEND,KLEV+1) :: UFLXC,DFLXC
          !do i = 2,2; call movave % add_data(id = c4(i)//"uflx", data = uflx(KST:KEND,2:KLEV+1), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist); enddo;
          !do i = 2,2; call movave % add_data(id = c4(i)//"dflx", data = dflx(KST:KEND,2:KLEV+1), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist); enddo;
          !call movave % add_data(id = "longitude", data = longitude(KST:KEND), ldr = ldrL, ddr = ddrL, ifexist = ifexist)
          !call movave % add_data(id = "latitude", data = latitude(KST:KEND), ldr = ldrL, ddr = ddrL, ifexist = ifexist)
          call movave % add_data(id = "lw_emissivity", data = PEMIS(KST:KEND), ldr = ldrL, ddr = ddrL, ifexist = ifexist)
          call movave % add_data(id = "sw_albedo", data = PALB(KST:KEND), ldr = ldrL, ddr = ddrL, ifexist = ifexist)
          call movave % add_data(id = "skin_temperature", data = PTS(KST:KEND), ldr = ldrL, ddr = ddrL, ifexist = ifexist)
          call movave % add_data(id = "cos_solar_zenith_angle", data = PMU0(KST:KEND), ldr = ldrL, ddr = ddrL, ifexist = ifexist)
          call movave % add_data(id = "pressure_hl", data = PAPRS(KST:KEND,2:KLEV+1), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "temperature_hl", data = ZTH(KST:KEND,2:KLEV+1), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "height_hl", data = RTAEZ(KST:KEND,1:KLEV,KPOSI), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)          
          call movave % add_data(id = "q", data = PQ(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "o3_mmr", data = PQO3(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "cloud_fraction", data = PNEB(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "q_ice", data = PQICE(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "q_liquid", data = PQLI(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "re_ice", data = ZRADIP(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "re_liquid", data = ZRADLP(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "overlap_param", data = ZOVLP_test(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "co2_vmr", data = CO2(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "ch4_vmr", data = CH4(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "n2o_vmr", data = N2O(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "cfc11_vmr", data = CFC11(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "cfc12_vmr", data = CFC12(KST:KEND,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "Ztauaer1", data = Ztauaer(KST:KEND,1,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "Ztauaer2", data = Ztauaer(KST:KEND,2,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "Ztauaer3", data = Ztauaer(KST:KEND,3,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "Ztauaer4", data = Ztauaer(KST:KEND,4,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)
          call movave % add_data(id = "Ztauaer5", data = Ztauaer(KST:KEND,5,1:KLEV), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist)

          
          !REAL(KIND=8) :: ZTAUCLD(KST:KEND,KLEV,16)  
          do k = 1,16
            if (k < 10) then
              write(c2,'(a,i1)') '0',k
            else
              write(c2,'(i2)') k
            end if
            do i = 2,2; call movave % add_data(id = c4(i)//"ZTAUCLD_"//c2, data = ZTAUCLD(KST:KEND,1:KLEV,k), ldr = ldrL3d, ddr = ddrL3d, ifexist = ifexist); enddo;
          end do
          
        end if
        end if

END SUBROUTINE SLM_RAD_DRIVER


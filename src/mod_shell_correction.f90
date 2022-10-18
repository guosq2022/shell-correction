module mod_shell_correction
   !use mod_global
   implicit none
   integer, parameter :: r8 = selected_real_kind(15, 307)
   integer, parameter :: r4 = selected_real_kind(6, 37)

   integer, parameter :: NDPLOT=90, NDBLOK=300, NDMAIN=30, NDXEKT=NDMAIN+1,     &
   &                     NDOMEG=2*NDMAIN+1
   integer, parameter :: NDENER=1771, NDLAGR=19, NDHERM=36, NDPARI=2,           &
   &                     NDFACT=40, NDMAUX=NDMAIN+1, NDINTG=NDMAIN+2,           &
   &                     NDNZMX=NDMAIN+1, NDROMX=int((real(NDMAIN)+1)/2),       &
   &                     NDFLPL=(NDMAIN+2)/2, NDLIM=NDMAIN+1

   !# single particle levels
   real(r8) :: SPL_N(NDENER), SPL_P(NDENER)

   !# wood-saxon
   real :: AOC, AOSO, ROC, ROSO, VO, VOSO, ROSOQ

   !# deformation parameters
   real :: beta1, beta2, beta3, beta4, beta5, beta6

   ! user defined configurations of wsbeta
   !======================================================================
   !         P R E S E T I N G   T H E   D E F A U L T   V A L U E S
   !======================================================================
   !         N A M E L I S T                             I N T D A T
   !======================================================================
   integer :: lcount    = 0
   integer :: KOMple    = 1
   integer :: iarray    = 0
   integer :: IDEcom    = 0
   integer :: INCrea    = 0
   integer :: ISPh      = 0
   integer :: iquadr    = 0
   integer :: ILPunl    = 19
   integer :: ILPunh    = 35
   integer :: iepsl     = 10
   integer :: iepsh     = 10
   integer :: NCUtl     = 0
   integer :: MIZol     = 0
   integer :: NPRint    = 0
   integer :: IPRint    = 25
   integer :: iplot     = 0
   integer :: ishape    = 0
   integer :: NPArit    = 2
   integer :: ICOuso    = 0
   integer :: iwrite    = 0

   !======================================================================
   !         N A M E L I S T                             D Y N A M C
   !======================================================================
   ! CONSTANTS BELOW BASED ON NUCL. PHYS. A412(1984)61
   real :: ECUtb        = 500.00
   real :: EDEcdo       = -10.0
   real :: EDEcgo       = -5.0
   real :: FINcla       = 0.144018
   real :: FINcro       = 0.276515
   real :: FACc         = 1.200000

   integer :: NMAxp     = 15
   integer :: NMAxn     = 15
   integer :: iendd     = 1
   integer :: ichoic    = 6

   integer :: iz, in, ia

   !# global
   integer :: IProt

   real :: EPSl, EPSH
   real :: CSTr, CMAle, UAVg

   real :: RSP(NDBLOK), COTsp(NDBLOK),  SITsp(NDBLOK), DRDtsp(NDBLOK)
   integer :: NRAd, NRAp, NSTep
   real :: R01, COT0 , SIT0 , DFF
   real :: XEKat(NDXEKT), UMAxx, VMAxx
   real :: XLAm

   real :: STU , CTU , RU , RU2
   integer :: IMIn(NDROMX,NDNZMX) , IPLU(NDROMX,NDNZMX) , NMAXI
   integer :: NDImom(NDOMEG,NDPARI) ,  NOSMAX(NDOMEG,NDPARI)
   real :: ELM(NDBLOK,NDBLOK)

contains



   !-------------------------------------------------------------------
   ! SINGLE PARTICLE ENERGIES, WAVE FUNCTIONS, QUADRUPOLE
   ! MOMENTS, AND G-FACTORS IN AXIALLY DEFORMED WOODS-SAXON POTENTIAL
   ! WITH APPLICATIONS IN THE TWO-CENTRE-TYPE NUCLEAR PROBLEMS.
   ! S. CWIOK, J. DUDEK, W. NAZAREWICZ, J. SKALSKI, T.R. WERNER.
   ! CREF. IN COMP. PHYS. COMMUN. 46 (1987) 379
   !-------------------------------------------------------------------
   SUBROUTINE WSBETA
      IMPLICIT NONE
      real :: AA, AKAppa, AN, axn, axp, AZ, baeta2
      real :: delt, delta, dr, enc, fae, gae
      real :: hompr, homz, pi, PIStar, pppp, r, sx, teta, u
      real :: uaoc, uaoso, ubar, uroc, urosoq, uvo, uxappa, uxlam
      real :: V0, vvvd, vvvg, wmaxs, wmins, x
      real :: zd, zg, zkrok
      INTEGER :: i, ib, ico, iczy, ilo, ilpuns, ilpunz, ILS
      INTEGER :: iw, kli, kq, lcou, lcpom, li, licz
      INTEGER :: model, ndnega, ndposi, nmax, nnnmm, nrapis, nrxp
      INTEGER :: nshap, nw
      DIMENSION pppp(NDPLOT), vvvd(NDPLOT), vvvg(NDPLOT), enc(NDPLOT), iczy(NDXEKT)
      DIMENSION ndposi(NDOMEG), ndnega(NDOMEG)
      NAMELIST /INTDAT/ lcount, KOMple, iarray, IDEcom, INCrea,     &
      & ISPh, iquadr, iepsl, iepsh, NCUtl, MIZol, NPRint, IPRint,   &
      & iplot, ishape, NPArit, ICOuso
      NAMELIST /DYNAMC/ ECUtb, EDEcdo, EDEcgo, FINcla, FINcro,      &
      & FACc, IZ, IN, NMAxp, NMAxn, iendd, ichoic
      NAMELIST /USERDF/ uvo, uxappa, uxlam, uroc, urosoq, uaoc, uaoso
      DATA axn, axp/939.5527, 938.2592/ , pi/3.14159265/
      DATA nrxp, nrapis/50, 250/ , nw/10/

      NSTep = 2

      NRAp = nrxp

      DO kli = 1 , NDXEKT
         XEKat(kli) = 0.0000
      ENDDO

      IF ( iarray==0 .and. iwrite/=0) WRITE (12,99001)

99001 FORMAT (///,1X,'THE WAVE FUNCTIONS ARE NOT RECORDED',//)
      IF ( iarray/=0 .and. iwrite/=0) WRITE (12,99002)
99002 FORMAT (///,1X,'THE WAVE FUNCTIONS ARE RECORDED',//)

      IF ( ISPh/=0 .and. iwrite/=0) THEN

         WRITE (12,99013)

         DO i = 1 , 5
            WRITE (12,99003)
99003       FORMAT ('+','THE FIXED SPHERICAL BASIS USED   ',            &
            &'CHECK WHETHER YOU REALLY NEED IT')
         ENDDO

         WRITE (12,99013)
      ENDIF

      model = 3
      licz = 0

      if(iwrite/=0) WRITE (12,99004) KOMple
99004 FORMAT (1X,'KOMPLET = ',I3)

      lcpom = lcount/2

      if(iwrite/=0) WRITE (12,99005) lcpom , licz
99005 FORMAT (///////,1X,'DATE  ',10('.'),'  TIME  ',10('.'),29X,       &
      &'THE LAST DATA SET RUN ...',I3,10X,                       &
      &'         RUN NUMBER ....',I3,///)

      !         SKIPPING THE LCOUNT/2 SETS OF DATA FROM THE INPUT FILE

      !IF ( lcount/=0 ) THEN
      !
      !         S K I P P I N G      D A T A
      !
      !DO iiii = 1 , lcpom
      !
      !READ (5,DYNAMC)
      !IF ( ichoic>4 ) READ (5,USERDF)
      !IF ( ichoic>4 ) READ (5,USERDF)
      !
      !ENDDO
      !ENDIF
      !======================================================================
      !       A C T U A L   D A T A   S E Q U E N C E   F O R   T H E   R U N
      !======================================================================
      !
100   licz = licz + 1
      !
      !======================================================================
      !         I N P U T   I N P U T   I N P U T   I N P U T   I N P U T
      !======================================================================
      !
      !READ (5,DYNAMC)
      !READ (11,DYNAMC)
      !
      IF ( licz==1 .and. iwrite/=0 ) WRITE (12,99006) ECUtb
99006 FORMAT (///,1X,'THE ENERGY CUT OFF PARAMETER ',F7.1)
      !
      IF ( INCrea==1 ) THEN
         !
         IF ( FINcla<0.0001 ) FINcla = 0.144018
         IF ( FINcro<0.0001 ) FINcro = 0.276515
         !
         DO i = 1 , 15
            !
            baeta2 = (i-1)*0.075
            fae = FFFF(baeta2)
            gae = GGGG(baeta2)
            !
            if(iwrite/=0) WRITE (12,99007) FINcla , FINcro , baeta2 , fae , gae
99007       FORMAT (1X,6('*'),'  FINCLA=',F6.3,'   FINCOR=',F6.3,       &
            &'   AT BETA2 =',F8.4,'   F =',F8.4,'   G =',F8.4)
            !
         ENDDO
      ENDIF
      !
      !
      ILS = NMAxp
      IF ( ILS<NMAxn ) ILS = NMAxn
      !
      !         EPSL=10.0**(-IEPSL),   EPSH=10.0**(-IEPSH)
      !
      EPSl = 10.0**(-iepsl)
      EPSh = 10.0**(-iepsh)
      !
      IF ( ABS(BETa3)>1.0E-9 .OR. ABS(BETa5)>1.0E-9 ) NSTep = 1
      IF ( ABS(BETa1)>1.0E-9  ) NSTep = 1
      !
      IF ( NPArit==1 ) NSTep = NPArit
      !
      !
      IA = IZ + IN
      AA = FLOAT(IA)
      AZ = FLOAT(IZ)
      AN = FLOAT(IN)

      !
      !          NEUTRON-PROTON LOOP (IW=1 NEUTRONS, IW=2 PROTONS)
      !
      DO iw = 1 , 2
         !
         IF ( (MIZol/=1) .OR. (iw/=2) ) THEN
            !
            lcount = lcount + 1
            IPRot = iw - 1
            !
            nmax = NMAxn
            IF ( IPRot==1 ) nmax = NMAxp
            !
            IF ( MOD(NRAp,2)/=0 ) NRAp = NRAp + 1
            !
            PIStar = FLOAT(NRAp)
            !
            CALL PAR(AKAppa, BETa2)
            !
            V0 = -ABS(VO)
            !
            IF ( (MIZol/=2) .OR. (iw/=1) ) THEN
               !
               if(iwrite/=0) WRITE (12,99008)
99008          FORMAT ('1')
            ENDIF
            !
            !
            !=======================================================================AAXX0525
            !          *********        WRITE   ON   DISC        **********
            !=======================================================================AAXX0527
            !
            if(iwrite/=0) WRITE (nw,*) model , KOMple

            !
            !=======================================================================AAXX0531
            !          *********     WRITE ON DISC (IN INTRO)    **********
            !=======================================================================AAXX0533
            !
            CALL INTRO(hompr,homz,nw,ILS,AKAppa)
            !
            IF ( (MIZol/=2) .OR. (iw/=1) ) THEN
               !
               nnnmm = (nmax+1)*(nmax+2)*(nmax+3)/6
               !
               CALL DIMDEF(nmax,homz,hompr,ndposi,ndnega)
               !
               !=======================================================================AAXX0543
               !          *********        WRITE   ON   DISC        **********
               !=======================================================================AAXX0545
               !

               if(iwrite/=0) WRITE (nw,*) nnnmm , NDOMEG , (ndnega(kq),kq=1,NDOMEG) , &
               & (ndposi(kq),kq=1,NDOMEG)
               !
               NRAd = NRAp + 1
               delt = -pi/NRAp
               teta = pi + delt
               x = COS(0.001)
               COTsp(1) = -x
               COTsp(NRAd) = x
               SITsp(1) = -SQRT(1.0-x*x)
               SITsp(NRAd) = SITsp(1)
               RSP(1) = RSUR(-x)
               RSP(NRAd) = RSUR(x)
               DRDtsp(1) = -DRSUR(-x)*SITsp(1)
               DRDtsp(NRAd) = -DRSUR(x)*SITsp(NRAd)
               !
               DO i = 2 , NRAp
                  x = COS(teta)
                  COTsp(i) = x
                  SITsp(i) = SQRT(1.-x*x)
                  RSP(i) = RSUR(x)
                  dr = DRSUR(x)
                  DRDtsp(i) = -dr*SITsp(i)
                  teta = teta + delt
               ENDDO
               !
               IF ( NSTep==2 ) NRAd = NRAp/2 + 2
               !
               CALL MAXRO(UMAxx,VMAxx,nrapis)
               !
               CALL EXTREM(XEKat,iczy,nmax,pppp,vvvg)
               !
               ubar = UAVg
               !
               CALL ELMAT(hompr,homz,nmax,CSTr,nw,ILS,AKAppa)
            ENDIF
         ENDIF
         !
         !
         !          NUCLEAR SHAPE PLOTTING
         !
      ENDDO
      !
      IF ( ishape/=0 ) THEN
         !
         delta = pi/64
         !
         DO i = 1 , 65
            teta = (i-1)*delta
            x = COS(teta)
            sx = SIN(teta)
            r = RSUR(x)
            pppp(i) = CSTr*x*r
            vvvg(i) = CSTr*sx*r
            vvvd(i) = -vvvg(i)
         ENDDO
         !
         if(iwrite/=0) WRITE (12,99009) BETa2 , BETa3 , BETa4 , BETa5 , BETa6
99009    FORMAT ('1',/,7X,'Z (IN FERMI)',6X,' THE NUCLEAR SHAPE   ',1X, &
         &'BETA2 = ',F6.3,'    BETA3 = ',F5.2,'    BETA4 = ',    &
         & F6.3,'    BETA5 = ',F5.2,'    BETA6 = ',F5.2,/)
         !
         ib = 1
         ico = 2
         ilo = 2
         !
         nshap = 50
         !
         CALL SHAPE(pppp,vvvd,vvvg,65,nshap)
      ENDIF
      !
      !
      IF ( iplot/=0 ) THEN
         !
         wmins = V0
         wmaxs = 0.0
         !
         zg = 3.0
         zd = -zg
         !
         ilpunz = 47
         zkrok = (zg-zd)/ilpunz
         !
         DO li = 1 , ilpunz
            enc(li) = zd + (li-1)*zkrok
         ENDDO
         !
         ico = 3
         ib = 1
         !
         wmins = -70.0
         wmaxs = +10.0
         !
         !         SECTION OF THE POTENTIAL FOR  V = CONST
         !
         DO li = 1 , ilpunz
            !
            u = enc(li)
            !
            vvvd(li) = POTF(u,0.1)
            vvvg(li) = POTF(u,0.5)
            pppp(li) = POTF(u,1.0)
            !
         ENDDO
         !
         if(iwrite/=0) WRITE (12,99010)
99010    FORMAT ('1',2X,'SECTION OF THE POTENTIAL FOR V = CONST       ',&
         & 20X,                                                   &
         &'    V=1.0    (A)          V=0.5    (B)          V=0.1  '&
         & ,'     (C)')
         !
         CALL WYKRES(pppp,vvvg,vvvd,ilpunz,ico,ilo,wmins,wmaxs,ib)
         !
         !         SECTION OF THE POTENTIAL FOR  U = CONST (U.GT.0.)
         !
         ilpuns = 0
         zkrok = zg/(ilpunz-1)
         !
         DO li = 1 , ilpunz
            !
            u = (li-1)*zkrok
            !
            enc(li) = u
            !
            IF ( u>=0. ) THEN
               !
               ilpuns = ilpuns + 1
               !
               vvvd(ilpuns) = POTF(0.1,u)
               vvvg(ilpuns) = POTF(0.5,u)
               pppp(ilpuns) = POTF(1.0,u)
            ENDIF
            !
         ENDDO
         !
         if(iwrite/=0) WRITE (12,99011)
99011    FORMAT ('1',//,2X,'SECTION OF THE POTENTIAL FOR U = CONST  ',  &
         & 26X,                                                   &
         &'   U=1.0   CASE (A),      U=0.5   CASE (B),     U=0.1  CASE (C)'&
         & )
         !
         CALL WYKRES(pppp,vvvg,vvvd,ilpuns,ico,ilo,wmins,wmaxs,ib)
         !
         !         SECTION OF THE POTENTIAL FOR  U = CONST (U.LT.0.)
         !
         ilpuns = 0
         !
         DO li = 1 , ilpunz
            !
            u = enc(li)
            !
            IF ( u>=0.0 ) THEN
               !
               ilpuns = ilpuns + 1
               !
               vvvd(ilpuns) = POTF(-0.10,u)
               vvvg(ilpuns) = POTF(-0.70,u)
               pppp(ilpuns) = POTF(-1.20,u)
            ENDIF
            !
         ENDDO
         !
         if(iwrite/=0) WRITE (12,99012)
99012    FORMAT ('1',//,2X,'SECTION OF THE POTENTIAL FOR U = CONST  ',  &
         & 26X,                                                   &
         &'   U=-1.2   CASE (A),    U=-0.7   CASE (B),   U=-0.1   CASE (C)'&
         & )
         !
         CALL WYKRES(pppp,vvvg,vvvd,ilpuns,ico,ilo,wmins,wmaxs,ib)
      ENDIF
      !
      !
      lcou = lcount/2
      !
      IF ( iendd==0 ) GOTO 100
      !
      !
      !
      CLOSE (nw)
99013 FORMAT (/////)
   END SUBROUTINE WSBETA

   SUBROUTINE INTRO(Hompr,Homz,Nw,ILS,AKAppa)
      IMPLICIT NONE
      REAL AA, agrad, AKAppa, AMN, AMP, BC, BS, c, R0
      REAL EMC2, ende, hom0, Hompr, Homz, r
      REAL U1 , U2 , vef
      INTEGER  ILS , JSIgn
      INTEGER Nw

      !         INTRO DEFINES ALL THE AUXILIARY PARAMETERS FOR W-S POTENTIAL

      !         U1, U2, CMALE, BS, BC, PISTAR, UAVG, HOM0, HOMZ, HOMPR  ETC.

      !         CALLED FROM WSBGIG (THE MAIN PROGRAM)

      real, parameter :: axn = 939.5527
      real, parameter :: axp = 938.2592

      JSIgn = ichoic

      AA = float(IA)
      R0 = FLOAT(iendd)

      AMN = axn
      AMP = axp

      EMC2 = AMP

      IF ( IPRot==0 ) EMC2 = AMN

      IF ( IPRot==1 ) GOTO 200

      !CALL STALE

      CALL ZER(U1,U2)

      CALL CVOLUM(CMAle)

      IF ( MOD(ILPunh,2)/=0 ) ILPunh = ILPunh + 1

      BC = BC1(CMAle)
      BS = BS1(CMAle)

      agrad = FACc
      UAVg = 0.000

      IF ( ABS(BETa1)<1.0E-9 ) THEN
         IF ( ABS(BETa3)<1.0E-9 .AND. ABS(BETa5)<1.0E-9 ) GOTO 100
      ENDIF

      CALL CENTR(UAVg,CMAle)

100   IF ( FACc<0.1 ) FACc = 1.15

      hom0 = FACc*41./AA**(1./3.)

      CALL HMZ(hom0,Homz,Hompr)

200   r = ROC
      CSTr = CMAle*r
      c = CSTr
      vef = AKAppa*VOSo
      ende = R0

      IF ( (MIZol/=2) .OR. (IPRot/=0) ) THEN

         if(iwrite/=0) WRITE (12,99001) IZ , IN , IA ,   &
         & VO , AOC , ROC , VOSo , AOSo , &
         & ROSo , BETa2 , BETa3 , BETa4 , BETa5 , BETa6 ,&
         & U1 , U2 , ICOuso , BC , Homz , Hompr ,        &
         & ILPunh , ILPunl , BS , CMAle , agrad , NMAxp ,&
         & NMAxn , ILS , IPRot , JSIgn , UAVg , ROSoq ,  &
         & hom0 , XLAm

99001    FORMAT (1X,129('*'),/,1X,'*',127(' '),'*',/,1X,'*',     &
         &' SHELL MODEL - DEFORMED WOODS-SAXON POTENTIAL',54X,   &
         &'  NUCLEUS Z=',I3,' N=',I3,' A=',I3,' *',/,1X,'*',     &
         &127(' '),'*',/,1X,'*',                                 &
         &' CENTRAL    POTENTIAL  -  DEPTH  VO = ',F5.1,21X,     &
         &'DIFFUSENESS AVO = ',F6.3,20X,'RADIUS RO = ',F6.2,' *',&
         & /,1X,'*',127(' '),'*',/,1X,'*',                        &
         &' SPIN-ORBIT POTENTIAL  -  DEPTH VSO = ',F5.1,21X,     &
         &'DIFFUSENESS ASO = ',F6.3,20X,'RADIUS RV = ',F6.2,' *',&
         & /,1X,'*',127(' '),'*',/,1X,'*',                        &
         &' D E F O R M A T I O N    BETA2 = ',F6.3,5X,          &
         & 'BETA3 = ',F6.3,5X,'BETA4 = ',F6.3,5X,'BETA5 = ',F7.4, &
         & 4X,'BETA6 = ',F6.3,10X,' *',/,1X,'*',127X,'*',/,1X,'*',&
         &' AUXILIARY PARAMETERS  ','   U1 =',F8.4,'   U2 =',    &
         & F7.4,'    COUSO =',I2,6X,'     BC =',F8.5,2X,'  HMZ=', &
         & F7.3,'  HMP=',F7.3,10X,' *',/,1X,'*',127X,'*',/,1X,'*',&
         & 26X,'IH =',I3,'   IL =',I3,'  BS =',F8.6,'  CVOLUM =', &
         & F8.5,' FACC =',F6.4,'    NMP=',I2,'   NMN=',I2,        &
         &'   ILS=',I2,11X,'*',/,1X,'*',127(' '),'*',/,1X,'*',   &
         & 26X,'IPROT = ',I1,'   CHOICE = ',I2,'   UBAR = ',F9.6, &
         & 3X,'     ROSO =',F6.4,'    HOMO = ',F7.3,'  XL=',F5.2, &
         & 10X,' *',/,1X,'*',127(' '),'*',/,1X,129('*'))

         IF ( JSIgn/=(-9) .and. (iwrite/=0)) WRITE (Nw,*) IZ , &
         & IN , IA , VO , AOC , ROC ,&
         & VOSo , AOSo , ROSo , BETa2 , BETa3 , &
         & BETa4 , BETa5 , U1 , U2 , R0 , BC ,  &
         & BS , Homz , Hompr , ILPunh , ILPunl ,&
         & c , CMAle , ende , NMAxp , NMAxn ,   &
         & ILS , IPRot , JSIgn , UAVg , vef ,   &
         & AKAppa , ECUtb , IDEcom , ICOuso ,   &
         & NPArit , Iarray
      ENDIF
   END SUBROUTINE INTRO

   SUBROUTINE STALE(A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,C0,C1,C2,C3,C4)
      IMPLICIT NONE
      REAL A0, A1, A2, A3, A4, A5, A6
      real B0, B1, B2, B3, B4, B5
      real C0, C1, C2, C3, C4
      REAL g1, g2, g3, g4, g5, g6, pi4

      !         CALLED FROM INTRO FOR DEFINING GEOMETRICAL CONSTANTS

      pi4 = real(16.*ATAN(1.0D0))

      g1 = SQRT(3./pi4)*BETa1
      g2 = SQRT(5./pi4)/2.*BETa2
      g3 = SQRT(7./pi4)/2.*BETa3
      g4 = SQRT(9./pi4)/8.*BETa4
      g5 = SQRT(11./pi4)/8.*BETa5
      g6 = SQRT(13./pi4)/16.*BETa6

      A0 = 1.0 - g2 + 3.*g4 - 5.0*g6
      A1 = 15.*g5 - 3.*g3 + g1
      A2 = 3.*g2 - 30.*g4 + 105.*g6
      A3 = 5.0*g3 - 70.0*g5
      A4 = 35.*g4 - 315.*g6
      A5 = 63.0*g5
      A6 = 231.*g6

      B0 = A1
      B1 = 2.*A2
      B2 = 3.*A3
      B3 = 4.*A4
      B4 = 5.*A5
      B5 = 6.*A6

      C0 = B1
      C1 = 2.*B2
      C2 = 3.*B3
      C3 = 4.*B4
      C4 = 5.*B5
   END SUBROUTINE STALE

   SUBROUTINE ZER(U1,U2)
      IMPLICIT NONE
      REAL s0, s1, s2, s3, s4, s5, s6, U1, U2

      !         CALLED FROM INTRO TO DEFINE EXTREME DIMENSIONS ALONG Z-AXIS

      s0 = SQRT(4.*3.1415926)
      s1 = SQRT(3.)/s0*BETa1
      s2 = SQRT(5.)/s0*BETa2
      s3 = SQRT(7.)/s0*BETa3
      s4 = SQRT(9.)/s0*BETa4
      s5 = SQRT(11.)/s0*BETa5
      s6 = SQRT(13.)/s0*BETa6

      U1 = -(1.-s1+s2-s3+s4-s5+s6)
      U2 = +(1.+s1+s2+s3+s4+s5+s6)

      IF ( U1>=0.0 ) THEN

         if(iwrite/=0) WRITE (12,99001) U1
99001    FORMAT (///,1X,'INCORRECT U1=',E12.5,///)

         STOP 'ZERU1U2'

      ELSEIF ( U2>0.0 ) THEN
         GOTO 99999
      ENDIF
      if(iwrite/=0) WRITE (12,99002) U2
99002 FORMAT (///,1X,'INCORRECT U2=',E12.5,///)
      STOP 'ZERU1U2'
99999 END SUBROUTINE ZER

   SUBROUTINE CVOLUM(Cmale_)
      IMPLICIT NONE
      REAL Cmale_
      DOUBLE PRECISION  cmall

      !         EXTERNAL CMAL
      !         CVOLUM COMPUTES THE CONSTANT 'CMALE'  -  FROM THE CONSTANT
      !         VOLUME CONDITION (CALLED FROM INTRO)

      CALL QG16(-1.0D0,1.0D0,CMAL,cmall)

      Cmale_ = real((2./cmall)**(1./3.))
   END SUBROUTINE CVOLUM

   FUNCTION CMAL(X)
      IMPLICIT NONE
      DOUBLE PRECISION CMAL , X , r

      !         CALLED FROM QG16 (IN CVOLUM)

      r = RSUR(SNGL(X))
      CMAL = r*r*r
   END FUNCTION CMAL

   SUBROUTINE HMZ(Hom0,Homz,Hompr)
      IMPLICIT NONE
      REAL Hom0 , Hompr , Homz , stos
      DOUBLE PRECISION z2 , ro2 , xd , xg

      !         EXTERNAL Z2AV , RO2AV
      !         CALLED FROM INTRO TO DEFINE THE OPTIMAL CHOICE OF THE H.O.
      !         POTENTIAL (THE  BASIS  GENERATOR  FOR THE DIAGONALISATION)

      xd = -1.0D0
      xg = +1.0D0

      Homz = Hom0
      Hompr = Hom0

      IF ( ISPh/=0 ) RETURN

      CALL QG16(xd,xg,Z2AV,z2)
      CALL QG16(xd,xg,RO2AV,ro2)

      ro2 = ro2 - z2
      stos = real(SQRT(ro2*0.5/z2))

      Hompr = Hom0*stos**(-1./3.)
      Homz = Hompr*stos
   END SUBROUTINE HMZ

   FUNCTION Z2AV(X)
      IMPLICIT NONE
      DOUBLE PRECISION Z2AV , r , X

      !         CALLED FROM HMZ (AUXILIARY ROUTINE)

      r = RSUR(SNGL(X))
      Z2AV = X*X*r**5
   END FUNCTION Z2AV

   FUNCTION RO2AV(X)
      IMPLICIT NONE
      DOUBLE PRECISION RO2AV , r , X

      !         CALLED FROM HMZ (AUXILIARY ROUTINE)

      r = RSUR(SNGL(X))
      RO2AV = r**5
   END FUNCTION RO2AV

   SUBROUTINE MAXRO(Ux,Vx,Nx)
      IMPLICIT NONE
      REAL delt, r, teta, u, Ux, uy, v, Vx, vy, x
      INTEGER i, nrad0, Nx
      DOUBLE PRECISION pi
      DATA pi/3.1415926531D0/

      !          MAXRO  COMPUTES  MAXIMUM VALUE OF RO VARIABLE ON NUCLEAR
      !          SURFACE. CALLED FROM SWBETA FOR THE LINE PRINTER PLOTING

      nrad0 = Nx + 1
      delt = real(-pi/Nx)
      teta = real(pi)
      v = 0.
      u = 0.

      DO i = 1 , nrad0
         x = COS(teta)
         r = RSUR(x)
         uy = r*x
         vy = r*SQRT(1.-x*x)
         IF ( vy>=v ) THEN
            v = vy
            u = uy
         ENDIF
         teta = teta + delt
      ENDDO

      Ux = u
      Vx = v
   END SUBROUTINE MAXRO

   SUBROUTINE PAR(Akappa, BBeta2)
      IMPLICIT NONE
      REAL aa, Akappa, akonc, alp, amc, arp, asrod, &
      & BBeta2, facl, facr, hc, r, rc
      REAL uaoc, uaoso, uroc, urosoq, uvo, uxappa, uxlam, &
      & voo,  x, xappa, xmred
      INTEGER i, ip, isig, kol, kor, nol, nor
      REAL ASYM, CSYM, XLAM0
      !*-- changed by guo 2018-12-28
      DOUBLE PRECISION xl , rs
      DIMENSION xl(6,6) , rs(6,6) , alp(2) , arp(2)

      NAMELIST /USERDF/ uvo , uxappa , uxlam , uroc , urosoq , uaoc ,   &
      & uaoso

      DATA hc , amc/197.32891 , 938.9059/
      DATA asrod , akonc , alp , arp/207.5 , 277.5 , 114.5 , 110.5 ,    &
      & 126.5 , 120.5/

      DATA xl/ + .108043086D-09 , -.229525357D-07 , -.420955282D-05 ,   &
      & +.316904416D-03 , +.225774528D+00 , +.961727722D+01 ,        &
      & +.349062307D-09 , -.216539156D-06 , +.518987481D-04 ,        &
      & -.640940855D-02 , +.460069902D+00 , +.179017995D+01 ,        &
      & 2*0.0D0 , -.177834815D-5 , .649404259D-3 , .517376308D-3 ,   &
      & .193547876D2 , 2*0.0D0 , -.365074656D-6 , -.991050669D-4 ,   &
      & .112741878D0 , .195050080D1 , 5*0.0D0 , 31.5D0 , 5*0.0D0 ,   &
      & 17.8D0/

      DATA rs/ - .111325573D-10 , +.309669940D-08 , +.273809305D-06 ,   &
      & -.961087376D-04 , -.393589891D-02 , +.232744244D+01 ,        &
      & -.147873623D-11 , +.393871857D-09 , +.566738200D-07 ,        &
      & -.810119821D-05 , -.337081241D-02 , +.131134118D+01 ,        &
      & 2*0.0D0 , .571084434D-7 , -.243003373D-4 , +.165347512D-2 ,  &
      & .147195693D+1 , 2*0.0D0 , .628993602D-8 , +.820360178D-5 ,   &
      & -.508831022D-2 , .157829502D+1 , 5*0.0D0 , 1.28D0 , 5*0.0D0 ,&
      & 0.932D0/

      aa = FLOAT(ia)
      xmred = (aa-1)/aa

      isig = -1
      IF ( Iprot==1 ) isig = +1

      Akappa = 0.5

      facl = FFFF(BBeta2)
      facr = GGGG(BBeta2)

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !+         CHOICE OF THE WOODS-SAXON POTENTIAL PARAMETERS +
      !+                                                        +
      !+         ICHOIC.EQ.0  -  WAHLBORN  PARAMETERS           +
      !+         ICHOIC.EQ.1  -      ROST  PARAMETERS           +
      !+         ICHOIC.EQ.2  -  CEPURNOV  PARAMETERS           +
      !+         ICHOIC.EQ.3  -  " NEW "   PARAMETERS           +
      !+         ICHOIC.EQ.4  -  UNIVERSAL PARAMETERS           +
      !+         ICHOIC.EQ.5  -  PARAMETERS   DEFINED           +
      !+                         VIA NAMELIST  USERDF           +
      !+         ICHOIC.EQ.6  -    WS3.2   PARAMETERS           +
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF ( Ichoic==0 ) THEN

         !        WAHLBORN PARAMETERS

         Vo = 51.
         voo = Vo

         Roc = 1.27*aa**(1./3.)
         Roso = 1.27*aa**(1./3.)*facr

         Aoc = 0.67
         Aoso = 0.67

         !        A C C I D E N T A L |||

         xappa = 0.67
         Vo = -Vo*(1.+isig*xappa*(In-Iz)/aa)

         XLAm = 32.0*facl
         Voso = -0.25*XLAm*Vo*(hc/amc)**2/xmred**2

         ROSoq = Roso*aa**(-1./3.)

         RETURN
      ENDIF

      !--------------------------------------------------------------
      IF ( Ichoic==1 ) THEN

         !        ROST PARAMETERS

         IF ( Iprot==1 ) THEN

            Roc = 1.275*aa**(1./3.)
            Roso = 0.932*aa**(1./3.)*facr
            XLAm = 17.8*facl
         ELSE

            Roc = 1.347*aa**(1./3.)
            Roso = 1.280*aa**(1./3.)*facr

            XLAm = 31.5*facl
         ENDIF

         Aoc = 0.70
         Aoso = 0.70

         Vo = 49.6
         voo = Vo

         xappa = 0.86
         Vo = -Vo*(1.+isig*xappa*(In-Iz)/aa)
         Voso = -0.25*XLAm*Vo*(hc/amc)**2/xmred**2
         ROSoq = Roso*aa**(-1./3.)

         RETURN
      ENDIF

      !----------------------------------------------------------------

      IF(Ichoic==2) THEN

         !        CHEPURNOV PARAMETERS

         Roc = 1.24*aa**(1./3.)
         Roso = 1.24*aa**(1./3.)*facr

         Aoc = 0.63
         Aoso = 0.63

         Vo = 53.3
         voo = Vo

         xappa = 0.63
         Vo = -Vo*(1.+isig*xappa*(In-Iz)/aa)
         XLAm = 23.8*(1.+2.*(In-Iz)/aa)*facl
         Voso = -.25*XLAm*Vo*(hc/amc)**2/xmred**2
         ROSoq = Roso*aa**(-1./3.)

         RETURN
      ENDIF

      !------------------------------------------------------------------

      IF ( Ichoic==3 ) THEN

         !        NEW PARAMETERS

         ip = Iprot + 1

         nol = 0
         nor = 0

         IF ( aa<=alp(ip) .OR. aa>=akonc ) nol = 2
         IF ( aa<=arp(ip) .OR. aa>=akonc ) nor = 2
         IF ( aa>=asrod .AND. aa<=akonc ) nol = 1
         IF ( aa>=asrod .AND. aa<=akonc ) nor = 1

         kol = ip + 2*nol
         kor = ip + 2*nor

         r = 0.
         x = 0.

         DO i = 1 , 6
            r = real(r*aa + rs(i,kor))
            x = real(x*aa + xl(i,kol))
         ENDDO

         rc = 1.347
         IF ( Iprot==1 ) rc = 1.275

         ROSoq = r*facr
         Roc = rc*aa**(1./3.)
         Roso = r*facr*aa**(1./3.)

         XLAm = x*facl

         Aoc = 0.70
         Aoso = 0.70

         Vo = 49.6
         voo = Vo

         xappa = 0.86
         Vo = -Vo*(1.+isig*xappa*FLOAT(In-Iz)/aa)
         Voso = -0.25*XLAm*Vo*(hc/amc)**2/xmred**2

         RETURN
      ENDIF

      !---------------------------------------------------------

      IF ( Ichoic==4 ) THEN

         !        UNIVERSAL PARAMETERS

         Vo = 49.6
         voo = Vo
         xappa = 0.86

         IF ( Iprot==0 ) XLAm = 35.0*facl
         IF ( Iprot==1 ) XLAm = 36.0*facl

         IF ( Iprot==0 ) ROSoq = 1.31*facr
         IF ( Iprot==1 ) ROSoq = 1.32*facr

         IF ( Iprot==0 ) Roc = 1.347
         IF ( Iprot==1 ) Roc = 1.275

         Aoc  = 0.70
         Aoso = 0.70

         Vo = -Vo*(1.0+isig*xappa*FLOAT(In-Iz)/aa)
         Voso = -0.25*XLAm*Vo*(hc/amc)**2/xmred**2

         Roc = Roc*aa**(1./3.)
         Roso = ROSoq*aa**(1./3.)

         RETURN
      ENDIF

      !----------------------------------------------------------------

      IF ( Ichoic==5 ) THEN

         !         PARAMETERS FROM CARDS

         !READ (5,USERDF)

         ROSoq = urosoq
         xappa = uxappa
         Aoso = uaoso
         XLAm = uxlam
         Roc = uroc
         Aoc = uaoc

         Vo = uvo
         voo = Vo

         Vo = -Vo*(1.+isig*xappa*FLOAT(In-Iz)/aa)
         Voso = -0.25*XLAm*Vo*(hc/amc)**2/xmred**2
         Roc = Roc*aa**(1./3.)
         Roso = ROSoq*aa**(1./3.)

         RETURN
      ENDIF
      !-------------------------------------------------------------
      IF(Ichoic==6) THEN
         ! Modified by SGuo
         ! WS3.2 parameters
         Aoc = real(0.7842029057d0) !a0
         Roc = real(1.3840356614d0) !r0
         CSYM = 29.2876
         xappa = 1.4492
         XLAM0 = 26.3163

         Vo = real(47.478721388d0)
         voo = Vo

         ROSoq = Roc
         ASYM = CSYM * (1. - xappa * aa**(-1./3.) + (2. - ABS(FLOAT(In-Iz)/aa))/(2. + aa * ABS(FLOAT(In-Iz)/aa)))

         IF(Iprot==0) XLAm = XLAM0 * (1. + FLOAT(In)/aa)
         IF(Iprot==1) XLAm = XLAM0 * (1. + FLOAT(Iz)/aa)

         Vo = -(Vo + isig * ASYM * FLOAT(In - Iz)/aa)
         Voso = -0.25 * XLAm * Vo * (hc/amc)**2/xmred**2
         Roc = Roc * aa**(1./3.)
         Roso = Roc
         Aoso = Aoc
         RETURN
      ENDIF
   END SUBROUTINE PAR

   FUNCTION POTF(U,V)
      IMPLICIT NONE
      real :: POTF, U, V

      !          POTF  COMPUTES  THE NUCLEAR PART OF THE DEFORMED WOODS-SAXON
      !          POTENTIAL. (CALLED FROM ELMAT FOR THE MATRIX ELEMENT DEFIN.)

      POTF = VO/(1.+EXP(ROC/AOC*EL(U,V)))
   END FUNCTION POTF

   FUNCTION SORF(U,V)
      IMPLICIT NONE
      real :: SORF, U, uu, V, vv

      !          SORF COMPUTES SPIN-ORBIT PART OF THE WOODS-SAXON POTENTIAL
      !          (CALLED  FROM  ELMAT  FOR  THE MATRIX ELEMENT  DEFINITION)

      uu = U/(ROSo/ROC)
      vv = V/(ROSo/ROC)

      SORF = VOSo/(1.+EXP(ROSo/AOSo*EL(uu,vv)))
   END FUNCTION SORF

   FUNCTION EL(Uu,V)
      IMPLICIT NONE
      REAL ctmt0, delta, delta2, df, dod2, EL, eps
      REAL f, fff, od2, od20, odl,  &
      & pi, PLUs, pq, R20, rr
      REAL stmt0, teta, tetal, tetap, tetast, &
      & tetau, u, Uu, V, x
      INTEGER i, id, id0, iend, ier, iter

      DATA eps , iend , pi/1.E-6 , 50 , 3.1415926531/

      !          EL COMPUTES DISTANCE OF A POINT  (UU,V)  FROM THE NUCLEAR
      !          SURFACE. SURFACE IS EXPRESSED IN THE DIMENSIONLESS  UNITS
      !          COORDINATES. EL  IS DEFINED WITH THE MINIUS SIGN WHEN THE
      !          POINT  BELONGS  TO  THE  INTERIA  OF THE NUCLEAR SURFACE.

      delta = pi/NRAp
      delta2 = delta*0.5 + 0.001
      u = Uu

      IF ( V<0.0 .and. (iwrite/=0)) WRITE (12,99001) u , Uu , V

99001 FORMAT (1X,'U,UU,V=',3E12.5)

      IF ( (NSTep==2) .AND. (Uu>0.) ) u = -Uu

      R20 = u*u + V*V
      R01 = SQRT(R20)

      COT0 = u/R01
      SIT0 = V/R01
      od20 = 10000.

      id0 = 1

      DO i = 2 , NRAd

         ctmt0 = COTsp(i)*COT0 + SITsp(i)*SIT0
         stmt0 = SITsp(i)*COT0 - SIT0*COTsp(i)
         dod2 = DRDtsp(i)*(RSP(i)-R01*ctmt0) + RSP(i)*R01*stmt0

         id = -1

         IF ( dod2>0.000 ) id = +1
         IF ( i==(NRAp+1) ) id = -1

         IF ( id0==1 .AND. id==-1 ) THEN


            tetap = pi - delta*(i-1.5)

            !         N E W T O N   M E T H O D   AS THE ORDINARY TREATMENT

            CALL RTNI(teta,f,df,MINIM,tetap,eps,iend,ier)
            plus = dff
            IF ( ier/=0 .OR. PLUs<0. .OR. ABS(teta-tetap)>delta2 ) THEN

               tetast = teta
               tetal = pi - delta*(i-1)
               tetau = pi - delta*(i-2)

               iter = 200

               !         WHEN THE NEWTON METHOD FAILS BECAUSE OF THE DERIVATIVE
               !         SINGULARITIES  THE  ORDINARY  BISECTION METHOD IS USED

               CALL BISEC(tetal,tetau,teta,FUNDEL,eps,iter,ier)

               fff = FUNDEL(teta)

               IF ( ier/=0 ) THEN

                  IF ( ABS(tetast-teta)>=10.0*eps .and. iwrite/=0) THEN

                     WRITE (12,99002) u , V , tetal , tetau , teta , fff
                     !     write (12,55)U,V,TETA,TETAP,F,DF,DELTA
99002                FORMAT (1X,'U,V,TETA,TETAP,F,DF,DELTA=',10F12.6)

                     WRITE (12,99003) ier
99003                FORMAT (///,1X,'INCORRECT IER=',I4,' FROM EL',///)
                  ENDIF
               ENDIF

            ENDIF

            !     STOP 'ELFAILED'

            x = COS(teta)
            rr = RSUR(x)
            od2 = R20 + rr*(rr-2.*R01*(x*COT0+SQRT(1.-x*x)*SIT0))

            IF ( od2<od20 ) od20 = od2

            id0 = id
         ELSE


            id0 = id
         ENDIF

      ENDDO

      IF ( od20<0.0 .and. iwrite/=0) WRITE (12,99004) Uu , u , V , &
      & tetal, tetau, teta, fff, od20

99004 FORMAT (1X,'UU,U,V,TETAL,TETAU,TETA,FFF,OD20=',7E12.5,2X,E12.5)

      IF ( od20<0.0 ) od20 = 0.0000001

      odl = SQRT(od20)
      rr = RSUR(COT0)
      pq = -1.

      IF ( R01>rr ) pq = 1.

      EL = odl*pq*CMAle
   END FUNCTION EL

   FUNCTION FUNDEL(X)
      IMPLICIT NONE
      REAL cot, ctmt0, dr, drdt, FUNDEL, r, sit, stmt0, X

      !         CALLED BY BISEC FROM EL (USED WHEN THE NEWTONS METHOD
      !         FAILS TO FIND THE SOLUTION)

      cot = COS(X)
      sit = SIN(X)

      ctmt0 = cot*COT0 + sit*SIT0
      stmt0 = sit*COT0 - SIT0*cot

      r = RSUR(cot)
      dr = DRSUR(cot)
      drdt = -sit*dr

      FUNDEL = drdt*(r-R01*ctmt0) + r*R01*stmt0
   END FUNCTION FUNDEL

   SUBROUTINE MINIM(X,F,Df)
      IMPLICIT NONE
      REAL cot, ctmt0, ddr, ddrdt, Df, dr,  &
      & drdt, F, r, sit, stmt0, X

      !         MINIM COMPUTES THE THETA-ANGLE FOR THE DISTANCE FUNCTION
      !         (CALLED FROM EL)

      cot = COS(X)
      sit = SQRT(1.-cot*cot)
      ctmt0 = cot*COT0 + sit*SIT0
      stmt0 = sit*COT0 - SIT0*cot
      r = RSUR(cot)
      dr = DRSUR(cot)
      drdt = -sit*dr
      ddr = D2RSUR(cot)
      ddrdt = sit*sit*ddr - dr*cot
      F = drdt*(r-R01*ctmt0) + r*R01*stmt0
      Df = R01*(ctmt0*(r-ddrdt)+stmt0*2.*drdt) + drdt*drdt + r*ddrdt
      DFF = Df
   END SUBROUTINE MINIM

   FUNCTION COUL(U,V)
      IMPLICIT NONE
      REAL AZ, COUL, ct, gc, r, ro
      REAL st, U, V

      !          COUL COMPUTES THE COULOMB POTENTIAL OF THE UNIFORM CHARGE
      !          DISTRIBUTION.

      AZ=float(IZ)
      ro = ROC
      gc = 3./(4.*3.1415926)*(AZ-1.)*CMAle**2*(197.32891/137.03602)/ro
      r = SQRT(U*U+V*V)
      st = V/r
      ct = U/r

      COUL = real(-gc*XICOUL(st,ct,r))
   END FUNCTION COUL

   FUNCTION XICOUL(St,Ct,R)
      IMPLICIT NONE
      REAL Ct, R, St
      DOUBLE PRECISION  pi, xi,XICOUL
      !EXTERNAL FPXI1

      DATA pi/3.1415926536D0/

      !         CALLED FROM COUL (AUXILIARY ROUTINE)

      STU = St
      CTU = Ct
      RU = R
      RU2 = R*R

      CALL QG32(0.0D0,pi,FPXI1,xi)
      XICOUL = xi
   END FUNCTION XICOUL

   FUNCTION FPXI1(Teta2)
      IMPLICIT NONE
      REAL fpxi
      INTEGER i1, i2
      DOUBLE PRECISION e, q, ak, ck2, FPXI1, Teta2, st1, ct1, &
      & r1, dr1, ctt1, stt1, cts, ctr, r12, rr1, a2, b2, ak2, a

      !         CALLED FROM QG32 IN XI1 (AUXILIARY ROUTINE)

      st1 = DSIN(Teta2)
      ct1 = DCOS(Teta2)

      r1 = RSUR(SNGL(ct1))
      dr1 = -st1*DRSUR(SNGL(ct1))

      ctt1 = CTU*ct1
      stt1 = STU*st1
      cts = ctt1 - stt1
      ctr = ctt1 + stt1
      r12 = r1*r1
      rr1 = 2.*RU*r1
      a2 = RU2 + r12 - rr1*cts
      b2 = RU2 + r12 - rr1*ctr
      ak2 = a2 - b2

      IF ( DABS(a2)<1.0D-10 ) a2 = 1.0D-10

      ak2 = ak2/a2
      ak = DSQRT(ak2)
      a = DSQRT(a2)
      ck2 = 1.0D0 - ak2

      CALL CEL2(e,ak,1.0D0,ck2,i1)
      CALL CEL2(q,ak,1.0D0,1.0D0,i2)

      fpxi = real(q/a*(r1*st1*(RU2-r12)+dr1*(rr1*CTU-ct1*(RU2+r12))))
      FPXI1 = fpxi + a*e*(dr1*ct1-r1*st1)
   END FUNCTION FPXI1

   FUNCTION BC1(Cmale_)
      IMPLICIT NONE
      REAL bc, BC1, Cmale_ , res

      DOUBLE PRECISION  ra , rb , rc , um1 , um2 , xg , pi

      DATA pi/3.1415926536D0/

      !         BC1 CALCULATES CONSTANT BC (THE RATIO OF THE COULOMB ENERGY
      !         OF  THE DEFORMED UNIFORM CHARGE DISTRIBUTION TO THAT OF THE
      !         SPHERICAL  DISTRIBUTION  OF  THE  SAME  CHARGE  AND  VOLUME

      !         CALLED FROM INTRO

      xg = pi

      IF ( ABS(BETa3)<1.0E-6 .AND. ABS(BETa5)<1.0E-6 ) xg = pi/2.0D0

      um1 = xg/3.
      um2 = 2.*um1

      CALL QG16(0.0D0,um1,FPBC1,ra)
      CALL QG16(um1,um2,FPBC1,rb)
      CALL QG16(um2,xg,FPBC1,rc)

      res = real(ra + rb + rc)
      bc = real(res*3./4./pi*DBLE(Cmale_)**5)
      BC1 = bc
      IF ( DABS(xg-pi)<1.0D-6 ) BC1 = 0.5*bc
   END FUNCTION BC1

   FUNCTION FPBC1(Teta)
      IMPLICIT NONE
      REAL r
      DOUBLE PRECISION FPBC1, Teta, st, ct, rd

      !         CALLED FROM QG16 IN BC1

      st = DSIN(Teta)
      ct = DCOS(Teta)

      r = RSUR(SNGL(ct))

      rd = DBLE(r)

      FPBC1 = -rd**3*st*XI1(st,ct,rd)
   END FUNCTION FPBC1

   FUNCTION XI1(St,Ct,R)
      IMPLICIT NONE
      DOUBLE PRECISION XI1 , pi , xi ,  R , St , Ct
      !EXTERNAL FPXI1

      DATA pi/3.1415926536D0/

      !         CALLED FROM FPBC1

      STU = real(St)
      CTU = real(Ct)
      RU = real(R)
      RU2 = real(R*R)

      CALL QG32(0.0D0,pi,FPXI1,xi)
      XI1 = xi
   END FUNCTION XI1

   FUNCTION BS1(Cmale_)
      IMPLICIT NONE
      REAL bs, BS1, Cmale_ , res

      DOUBLE PRECISION  ra , rb , um1 , xd , xg

      !         EXTERNAL FPBS1

      !         BS1  CALCULATES CONSTANT BS (THE RATIO OF THE SURFACE ENERGY
      !         OF A DEFORMED BODY TO THAT OF THE SPHERICAL BODY OF THE SAME
      !         VOLUME

      !         CALLED FROM INTRO

      xg = +1.0D0

      IF ( ABS(BETa3)<1.0E-6 .AND. ABS(BETa5)<1.0E-6 ) xg = 0.0D0

      xd = -1.0D0
      um1 = (xg-xd)*0.5

      CALL QG32(xd,um1,FPBS1,ra)
      CALL QG32(um1,xg,FPBS1,rb)
      res = real(ra + rb)
      bs = real(res*DBLE(Cmale_)**2)
      BS1 = bs
      IF ( DABS(xg-1.0D0)<0.1D0 ) BS1 = 0.5*bs
   END FUNCTION BS1

   FUNCTION FPBS1(X)
      IMPLICIT NONE
      DOUBLE PRECISION FPBS1 , X , r , dr

      !         CALLED FROM QG32 IN BS1

      r = RSUR(SNGL(X))
      dr = DRSUR(SNGL(X))

      FPBS1 = r*DSQRT(r*r+dr*dr*(1.0D0-X*X))
   END FUNCTION FPBS1

   FUNCTION CEN(X)
      IMPLICIT NONE
      REAL r !, RSUR
      DOUBLE PRECISION CEN , X

      !         CALLED FROM QG16 IN CENTR

      r = RSUR(SNGL(X))
      CEN = X*DBLE(r)**4
   END FUNCTION CEN

   SUBROUTINE CENTR(Centr_,Cmale_)
      IMPLICIT NONE
      REAL Cmale_
      DOUBLE PRECISION cent
      !DOUBLE PRECISION CEN , cent
      REAL Centr_

      !         EXTERNAL CEN

      !         CENTR CALCULATES CONSTANT "UAVG" TO TAKE INTO ACCOUNT
      !         THE CONSTANCY  OF THE CNTRE OF MASS POSITION WHEN THE
      !         ODD-MULTIPOLE DEFORMATIONS ARE USED

      !         CALLED FROM INTRO

      CALL QG16(-1.0D0,+1.0D0,CEN,cent)
      Centr_ = real(cent*DBLE(Cmale_)**3*3./8.)
      RETURN
   END SUBROUTINE CENTR

   FUNCTION RSUR(X)
      IMPLICIT NONE
      REAL A0, A1, A2, A3, A4, A5, A6, RSUR, rsurr, X
      real B0, B1, B2, B3, B4, B5
      real C0, C1, C2, C3, C4

      !          RSUR COMPUTES R=R(COS(THETA)), X=COS(THETA)

      call STALE(A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,C0,C1,C2,C3,C4)

      rsurr = A0 + X*(A1+X*(A2+X*(A3+X*(A4+X*(A5+X*A6)))))

      IF ( rsurr<(-0.0000001) ) THEN

         if(iwrite/=0) WRITE (12,99001) rsurr
99001    FORMAT (///,1X,'RSUR = ',F10.5,///)
         STOP 'RSURNEGA'

      ENDIF

      RSUR = rsurr
   END FUNCTION RSUR

   FUNCTION DRSUR(X)
      IMPLICIT NONE
      REAL B0, B1, B2, B3, B4, B5, DRSUR, X
      real A0, A1, A2, A3, A4, A5, A6
      real C0, C1, C2, C3, C4

      !          DRSUR COMPUTES THE FIRST DERIVATIVE OF R(X)

      call STALE(A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,C0,C1,C2,C3,C4)

      DRSUR = B0 + X*(B1+X*(B2+X*(B3+X*(B4+X*B5))))
   END FUNCTION DRSUR

   FUNCTION D2RSUR(X)
      IMPLICIT NONE
      REAL C0, C1, C2, C3, C4, D2RSUR, X
      real A0, A1, A2, A3, A4, A5, A6
      real B0, B1, B2, B3, B4, B5

      !          D2RSUR COMPUTES THE CORRESPONDING SECOND DERIVATIVE

      call STALE(A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,C0,C1,C2,C3,C4)

      D2RSUR = C0 + X*(C1+X*(C2+X*(C3+X*C4)))
   END FUNCTION D2RSUR

   SUBROUTINE QUANTN(Dom,Iparr,Nmax,Homz,Hompr,Nmain, &
   & Nzmain,Nrmain,Lbmain,Ismain,Idim)
      IMPLICIT NONE
      REAL ebase1, Hompr, Homz
      INTEGER I, Idim, IMAx,  Iparr, isig, isigp, Ismain,  &
      & j, lamb, lambp, Lbmain, maxim, n
      INTEGER  Nmain, Nmax, nmaxpn, np, npr, nprp, &
      & Nrmain, nro, nrop, nz, Nzmain, nzp
      DIMENSION Nmain(Idim), Nzmain(Idim), Lbmain(Idim), Nrmain(Idim), &
      & Ismain(Idim), IMAx(NDROMX,NDNZMX)
      INTEGER Dom

      !         THIS  SUBROUTINE  GENERATES  THE QUANTUM NUMBERS
      !         OF THE HARMONIC OSCILLATOR BASIS AND DENUMERATES
      !         THE CORRESPONDING  BASIS  VECTORS  1, 2, 3,  ...

      nmaxpn = Nmax + 1

      I = 0

      DO isigp = 1 , 3 , 2

         isig = isigp - 2
         lamb = (Dom-isig)/2

         IF ( lamb<=Nmax ) THEN

            lambp = lamb + 1

            DO nprp = lambp , nmaxpn , 2

               npr = nprp - 1
               nro = (npr-lamb)/2
               nrop = nro + 1

               IF ( nrop>NDROMX ) THEN

                  if(iwrite/=0) WRITE (12,99001) nrop , NDROMX

99001             FORMAT (///,1X,'NROP =',I5,                           &
                  &'   EXCEEDS DIMENSION BOUND =',I5,            &
                  &'  IN QUANT,  STOP BY THE PROGRAM')

                  STOP 'NROQUANT'

               ENDIF

               DO np = nprp , nmaxpn

                  n = np - 1

                  IF ( Nstep/=1 ) THEN

                     IF ( (-1)**n/=(-1)**Iparr ) GOTO 10
                  ENDIF


                  nz = n - npr
                  nzp = nz + 1

                  IF ( nzp>NDNZMX ) THEN

                     if(iwrite/=0) WRITE (12,99002) nzp , NDNZMX

99002                FORMAT (///,1X,'NZP =',I5,                         &
                     &'   EXCEEDS DIMENSION BOUND =',I5,         &
                     &'  IN QUANT,  STOP BY THE PROGRAM')

                     STOP 'NZQUANT'

                  ENDIF

                  ebase1 = Homz*(nz+.5) + Hompr*(npr+1.)
                  !     write (12,9999)N,NZ,LAMB,NRO,EBASE1
                  !9999 FORMAT(1X,'N=',I2,'  NZ=',I2,'  LAMB=',I2,'  NRO=',I2,'  EBASE=',
                  !    *         F8.4)

                  IF ( ebase1<=Ecutb ) THEN

                     I = I + 1

                     Nmain(I) = n
                     Nzmain(I) = nz
                     Lbmain(I) = lamb
                     Nrmain(I) = nro
                     Ismain(I) = isig

                     IF ( isig<=0 ) THEN

                        IMIN(nrop,nzp) = I
                     ELSE

                        IMAX(nrop,nzp) = I
                        IPLU(NROP,NZP)=IMAX(NROP,NZP)
                     ENDIF
                  ENDIF

10             ENDDO
            ENDDO
         ENDIF
      ENDDO

      NMAXI=I
      IF ( I>Ndblok ) THEN

         if(iwrite/=0) WRITE (12,99003) Dom , I , Ndblok
99003    FORMAT (///,1X,'ACTUAL DIMENSION FOR THE OMEGA=',I2,'  BLOCK', &
         &'   I =',I3,' EXCEEDS DIMENSION BOUND =',I3)

         STOP 'DIMBOUND'

      ENDIF

      IF ( I>0 ) THEN

         maxim = -1000

         DO j = 1 , I
            IF ( maxim<Nmain(j) ) maxim = Nmain(j)
         ENDDO

         IF ( Nstep==1 ) THEN

            NDImom(Dom,1) = I
            NOSMAX(Dom,1) = maxim

         ELSE

            NDImom(Dom,Iparr) = I
            NOSMAX(Dom,Iparr) = maxim

         ENDIF

      ENDIF
   END SUBROUTINE QUANTN

   SUBROUTINE DIMDEF(Nmax,Homz,Hompr,Ndposi,Ndnega)
      IMPLICIT NONE
      REAL ebase1, Hompr, Homz
      INTEGER i, iparr, isig, isigp, lamb, lambp, n, &
      & Ndnega, Ndposi, Nmax, nmaxpn, np, npr,       &
      & nprp, nro, nrop, nz
      INTEGER nzp
      DIMENSION Ndposi(NDOMEG), Ndnega(NDOMEG)
      INTEGER dom

      !         THIS  SUBROUTINE  DEFINES THE SIZES OF THE OMEGA-BLOCKS
      !         AFTER CUTTING OFF THE HARMONIC OSCILLATOR BASIS

      nmaxpn = Nmax + 1

      DO dom = 1 , 2*Nmax + 1 , 2

         Ndnega(dom) = 0
         Ndposi(dom) = 0

         DO iparr = 1 , Nstep

            i = 0

            DO isigp = 1 , 3 , 2

               isig = isigp - 2
               lamb = (dom-isig)/2

               IF ( lamb<=Nmax ) THEN

                  lambp = lamb + 1

                  DO nprp = lambp , nmaxpn , 2

                     npr = nprp - 1
                     nro = (npr-lamb)/2
                     nrop = nro + 1

                     DO np = nprp , nmaxpn

                        n = np - 1

                        IF ( Nstep/=1 ) THEN
                           !
                           IF ( (-1)**n/=(-1)**iparr ) GOTO 2
                        ENDIF

                        nz = n - npr
                        nzp = nz + 1

                        ebase1 = Homz*(nz+.5) + Hompr*(npr+1.)

                        IF ( ebase1<=Ecutb ) i = i + 1

2                    ENDDO
                  ENDDO
               ENDIF
            ENDDO

            IF ( iparr/=1 ) Ndposi(dom) = i
            IF ( iparr==1 ) Ndnega(dom) = i

         ENDDO
      ENDDO
   END SUBROUTINE DIMDEF

   SUBROUTINE ELMAT(Hompr,Homz,Nmax,Cstret,Nw,ILS,AKAppa)
      IMPLICIT NONE
      REAL a1op, a1po, AKAppa, alamb, ambm, ancoef, anorm, &
      & ao, aoo, aop, apo, azcoef
      REAL b, b1op, b1po, boo, c, cantip, coulmb, cousor, &
      & cparal, Cstret, dfl, dfll, dfll1, dflpl, dher,    &
      & dhompr, dhomz
      REAL E, ebase1, ebase2, ENC, eta, etace, etaso,     &
      & fh, fl, fll
      REAL fll1, flpl, gfact, hc, hc2, her, Hompr, Homz,  &
      & pf, sqeta, q20, qqq20, sf
      REAL sqpr, sqz, suma, vsop1, vsop2, VVI, wh, wl,    &
      & x, xco, xh, xiro
      REAL xiz, xksi, xl, xlam0, xlam1, xmassn, xnnro,    &
      & xnrso1, xnrso2, xnz, xso, xso1, xx, xxl,          &
      & y, ynnro, zeta
      INTEGER i, ial1, ial2, ial3, ial4,    &
      & idim, ih, il, il1, il2, il3, il4, ilam,  &
      & ilamm, ILLim, ilm1
      INTEGER ilm2, ilpom, ILS, im, imm,  &
      & iorder, iparr, ipi, ipi1, ipo, ipop, ippar
      INTEGER  isig, isig1, isigp, isigp1, ismain,    &
      & j, JSIgn, k, kbas, KIJ, kik, kres, l, l1, l2, l3
      INTEGER l4, laa, lamb, lamb1, lambm, lambmp, lambp, &
      & lambp1, lambx, lbmain, licz, liczk, lpha, lphace, &
      & lphaso, m, n, n1
      INTEGER ndol, ndoll, ndsub, newi, ngor, nh
      INTEGER nho, nho1, nif, nl, nmain, Nmax, nmaxd,     &
      & nmaxpn, np, np1, npr, npr1, nprp, nprp1
      INTEGER nrmain, nro, nro1, nrom, nromp, nrop, nrop1, &
      & numer, Nw, nz, nz1, nzm, nzmain, nzmd,             &
      & nzmp, nznak, nzp, nzp1

      DOUBLE PRECISION fact, sfact, sq, xexp2, exp2, pi4
      INTEGER dom, domm
      CHARACTER*2 hplus, hblank, hminus, hplue, hblane, hminue, &
      & iznak(6)
      DIMENSION nmain(NDBLOK), nzmain(NDBLOK), lbmain(NDBLOK),  &
      & nrmain(NDBLOK), ismain(NDBLOK), q20(NDBLOK) ,           &
      & cparal(NDBLOK), cantip(NDBLOK), gfact(NDBLOK),          &
      & ancoef(NDBLOK,NDMAUX), azcoef(NDBLOK,NDMAUX)
      DIMENSION qqq20(NDENER), iorder(NDENER)
      DIMENSION sfact(NDFACT), xnz(NDFACT), sq(NDFACT), fact(64)
      DIMENSION fh(NDHERM,NDMAUX), her(NDINTG), xh(NDHERM),     &
      & wh(NDHERM), dher(NDINTG), xl(NDLAGR), wl(NDLAGR),       &
      & b(NDLAGR), c(NDLAGR), apo(NDLAGR), aop(NDLAGR), ao(NDLAGR)
      DIMENSION flpl(NDFLPL,NDLAGR,NDMAUX), &
      & dflpl(NDFLPL,NDLAGR,NDMAUX), fl(NDFLPL), dfl(NDFLPL)
      DIMENSION pf(NDLAGR,NDHERM), sf(NDLAGR,NDHERM),           &
      & xxl(NDINTG,NDLAGR)
      DIMENSION E(NDBLOK), VVI(NDBLOK), ENC(NDENER),            &
      & KIJ(NDENER), ILLim(NDENER)

      DATA hplus, hblank, hminus/'+ ' , '  ' , '- '/ , hplue,   &
      & hblane, hminue/'+E' , ' E' , '-E'/

      real, allocatable :: aqx(:,:)
      real, allocatable :: ZHC(:,:,:), ZHSo(:,:,:)
      allocate(aqx(NDBLOK,NDBLOK))
      allocate(ZHC(NDLAGR,NDINTG,NDINTG), ZHSo(NDLAGR,NDINTG,NDINTG))

      !         ILPUNL  -  SPECIFIES NUMBER  OF POINTS IN GAUSS-LAGUERRE
      !                    INTEGRATION FORMULAS
      !         ILPUNH  -  SPECIFIES NUMBER OF POINTS IN GAUSS-HERMITE
      !                    INTEGRATION FORMULAS
      !
      !         EPSL    -  ACCURACY IN THE GAUSS-LAGUERRE CALCULATIONS
      !         EPSH    -  ACCURACY IN THE GAUSS-HERMITE  CALCULATIONS
      !
      !
      !         MATRIX DECLARATION
      !
      !
      !         IF THE MAXIMUM VALUES FOR ILPUNL, ILPUNH AND NMAX ARE
      !         P, Q  AND R RESPECTIVELY, THEN DIMENSIONING SHOULD BE
      !         THE  FOLLOWING  ( MATRICES  ARE DIMENSIONED IN ELMAT)
      !
      !
      !
      !          FL((R+2)/2),DFL((R+2)/2),XL(P),WL(P), A(P),B(P),C(P)
      !          APO(P),AOP(P),AO(P),   XXL(R+2,P),   PF(P,Q),SF(P,Q)
      !          ZHC(P,R+2,R+2), ZHSO(P,R+2,R+2),  FH(Q,R+2),HER(R+2)
      !          XH(Q),WH(Q),FLPL((R+2)/2,P,R+1),DFLPL((R+2)/2,P,R+1)
      !          IMIN(P/2,R+1),IPLU(P/2,R+1)
      !
      !
      !
      !         T A B L E   O F   D I M E N S I O N S   (HAMILTONIAN)
      !
      !
      !         NMAX =  1         DIMENSION =  3
      !         NMAX =  2         DIMENSION =  6
      !         NMAX =  3         DIMENSION = 10
      !         NMAX =  4         DIMENSION = 15
      !         NMAX =  5         DIMENSION = 21
      !         NMAX =  6         DIMENSION = 28
      !         NMAX =  7         DIMENSION = 36
      !         NMAX =  8         DIMENSION = 45
      !         NMAX =  9         DIMENSION = 55
      !
      !         NMAX = 10         DIMENSION = 66          TOTAL = 286
      !         NMAX = 11         DIMENSION = 78          TOTAL = 364
      !         NMAX = 12         DIMENSION = 91          TOTAL = 455
      !         NMAX = 13         DIMENSION =105          TOTAL = 560
      !         NMAX = 14         DIMENSION =120          TOTAL = 680
      !         NMAX = 15         DIMENSION =136          TOTAL = 816
      !         NMAX = 16         DIMENSION =153          TOTAL = 969
      !         NMAX = 17         DIMENSION =171          TOTAL =1140
      !         NMAX = 18         DIMENSION =190          TOTAL =1330
      !         NMAX = 19         DIMENSION =210          TOTAL =1540
      !         NMAX = 20         DIMENSION =231          TOTAL =1771
      !
      !
      JSIgn = ichoic

      idim = NDBLOK

      fact(1) = 1.0D0
      DO i = 2 , 42
         fact(i) = (i-1)*fact(i-1)
      ENDDO

      DO i = 1 , NDENER
         iorder(i) = i
      ENDDO

      IF ( ILPunl>NDLAGR ) THEN
         if(iwrite/=0) WRITE (12,99001) ILPunl
99001    FORMAT (///,1X,'ILPUNL TOO LARGE =',I4,///)
         STOP 'NDLAGR'

      ELSEIF ( ILPunh<=NDHERM ) THEN

         IF ( Nmax>NDMAIN ) THEN

            !         CHECK THIS CASE WHEN INCREASING NDMAIN

            IF ( Nmax>NDMAIN+4 .OR. Nstep/=2 ) THEN

               IF ( IPRot==0 .and. iwrite/=0) WRITE (12,99002) Nmax
99002          FORMAT (///,1X,'PROTON  NMAX =',I5,'  TOO LARGE')

               IF ( IPRot==1 .and. iwrite/=0) WRITE (12,99003) Nmax
99003          FORMAT (///,1X,'NEUTRON NMAX =',I5,'  TOO LARGE')

               STOP 'NDMAX'
            ENDIF
         ENDIF

         numer = 0

         dhompr = 0.5*Hompr
         dhomz = 0.5*Homz
         xmassn = 938.90590
         hc = 197.32891
         alamb = hc/xmassn

         cousor = alamb**2*(hc/137.03602)

         iznak(1) = hplus
         iznak(2) = hblank
         iznak(3) = hminus
         iznak(4) = hplue
         iznak(5) = hblane
         iznak(6) = hminue

         vsop1 = AKAppa*xmassn/(hc**2)*SQRT(Homz*Hompr)
         vsop2 = AKAppa*xmassn/(hc**2)*Hompr

         nl = ILPunl
         nh = ILPunh
         nmaxpn = Nmax + 1
         nrom = Nmax/2
         nromp = nrom + 1
         lambm = Nmax
         lambmp = lambm + 1
         nzm = Nmax
         nzmp = nzm + 1

         DO i = 1 , idim
            DO j = 1 , idim
               ELM(i,j) = 1.E-10
            ENDDO
         ENDDO

         IF ( NPRint/=0 .and. iwrite/=0) WRITE (12,99014)

         !          ZEROES AND WEIGHTS FOR HERMITE INTEGRATIONS

         CALL HERMIT(nh,xh,wh,EPSh,fact)

         nzmd = nzmp + 1

         DO i = 1 , nh

            x = xh(i)

            IF ( Nstep==2 ) wh(i) = 2.0*wh(i)

            !         HERMITE FUNCTIONS AND DERIVATIVES

            CALL DHEP(x,nzmp,her,dher,NDINTG)
            !
            DO j = 1 , nzmd
               fh(i,j) = her(j)
            ENDDO
         ENDDO

         !         I  -  SPECIFIES POINT ,  J  -  SPECIFIES ORDER OF THE POLYN.

         sqpr = hc/SQRT(xmassn*Hompr)/Cstret
         sqz = hc/SQRT(xmassn*Homz)/Cstret

         hc2 = hc*hc
         xiz = hc2/xmassn/Homz
         xiro = hc2/xmassn/Hompr

         sfact(1) = 1.0D0
         sq(1) = 0.0D0

         exp2 = 1.0D0

         pi4 = 1.0D0/3.1415926531D0**(0.25D0)
         xexp2 = 1.0D0/DSQRT(2.0D0)

         DO i = 2 , NDFACT
            sq(i) = DSQRT(DBLE(FLOAT(i-1)))
            sfact(i) = sfact(i-1)*sq(i)
            exp2 = exp2*xexp2
            xnz(i) = SNGL(exp2*pi4/sfact(i))
         ENDDO

         xnz(1) = real(pi4)

         ambm = 0.
         lambx = 0

         DO nif = 1 , nl
            b(nif) = lambx + 2*nif - 1.
            c(nif) = (nif-1.)*(lambx+nif-1.)
         ENDDO

         !          ZEROES AND WEIGHTS FOR GAUSS-LAGUERRE INTEGRATIONS

         CALL LAGUER(nl,xl,wl,ambm,b,c,fact,EPSl)

         DO il = 1 , nl

            sqeta = sqpr*SQRT(xl(il))

            DO ih = 1 , nh

               IF ( xh(ih)<=0. .OR. Nstep/=2 ) THEN

                  !          CENTRE OF MASS TRANSFORMATION

                  zeta = sqz*xh(ih) + UAVg
                  !
                  pf(il,ih) = POTF(zeta,sqeta)
                  sf(il,ih) = SORF(zeta,sqeta)

                  IF ( IPRot/=0 ) THEN

                     coulmb = COUL(zeta,sqeta)

                     pf(il,ih) = pf(il,ih) + coulmb

                     IF ( ICOuso==1 ) sf(il,ih) = sf(il,ih)             &
                     & + cousor*coulmb
                  ENDIF

                  IF ( NPRint/=0 .and. iwrite/=0 ) THEN
                     WRITE (12,99004) il , ih , sqeta , zeta , pf(il,ih)&
                     & , sf(il,ih)
99004                FORMAT (1X,'FUNCT  ',2I3,4(1X,E16.9))
                  ENDIF
               ENDIF

            ENDDO
         ENDDO

         IF ( NPRint/=0 .and. iwrite/=0) WRITE (12,99014)

         DO ilam = 1 , lambmp
            ilamm = ilam - 1
            DO il = 1 , nl
               xx = xl(il)

               !          LAGUERRE FUNCTIONS AND DERIVATIVES STORED IN FL AND DFL

               CALL DLAP(xx,fl,dfl,nrom,ilamm,NDFLPL)

               DO im = 1 , nromp
                  imm = im
                  flpl(im,il,ilam) = fl(imm)
                  dflpl(im,il,ilam) = (ilamm-xx)*fl(imm) + 2.*dfl(imm)
               ENDDO
            ENDDO
         ENDDO

         nmaxd = nmaxpn + 1

         DO ilam = 1 , nmaxd
            ilpom = ilam - 2
            DO il = 1 , nl
               xxl(ilam,il) = xl(il)**ilpom
            ENDDO
         ENDDO

         DO n = 1 , nmaxd
            DO il = 1 , nl

               x = 0.0
               y = 0.0

               DO ih = 1 , nh
                  xksi = xh(ih)
                  IF ( xksi<=0.0 .OR. Nstep/=2 ) THEN
                     boo = fh(ih,n)*fh(ih,n)*wh(ih)
                     x = x + pf(il,ih)*boo
                     y = y + sf(il,ih)*boo
                  ENDIF
               ENDDO

               ZHC(il,n,n) = x
               ZHSo(il,n,n) = y

            ENDDO
         ENDDO

         IF ( Nstep/=2 ) THEN

            DO n = 2 , nmaxd
               m = n - 1
               DO il = 1 , nl

                  x = 0.0
                  y = 0.0

                  DO ih = 1 , nh

                     xksi = xh(ih)

                     IF ( xksi<=0.0 .OR. Nstep/=2 ) THEN

                        boo = fh(ih,n)*fh(ih,m)*wh(ih)
                        x = x + pf(il,ih)*boo
                        y = y + sf(il,ih)*boo
                     ENDIF

                  ENDDO

                  ZHC(il,n,m) = x
                  ZHC(il,m,n) = x
                  ZHSo(il,n,m) = y
                  ZHSo(il,m,n) = y

               ENDDO
            ENDDO
         ENDIF

         !          RECURRENCE RELATIONS WITH HERMITE POLYNOMIALS

         DO n = 3 , nmaxd , 2

            j = nmaxd - n + 1

            DO m = 1 , j

               k = n + m - 1

               DO il = 1 , nl

                  ZHC(il,k,m) = ZHC(il,k-1,m+1) - 2.0*(k-2)             &
                  & *ZHC(il,k-2,m)
                  ZHC(il,m,k) = ZHC(il,k,m)
                  ZHSo(il,k,m) = ZHSo(il,k-1,m+1) - 2.0*(k-2)           &
                  & *ZHSo(il,k-2,m)
                  ZHSo(il,m,k) = ZHSo(il,k,m)

                  IF ( m/=1 ) THEN

                     ZHC(il,k,m) = ZHC(il,k,m) + 2.0*(m-1)              &
                     & *ZHC(il,k-1,m-1)
                     ZHC(il,m,k) = ZHC(il,k,m)
                     ZHSo(il,k,m) = ZHSo(il,k,m) + 2.0*(m-1)            &
                     & *ZHSo(il,k-1,m-1)
                     ZHSo(il,m,k) = ZHSo(il,k,m)
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         !          RECURRENCE RELATIONS FOR ODD MULTIPOLES

         IF ( Nstep/=2 ) THEN

            DO n = 4 , nmaxd , 2

               j = nmaxd + 1 - n

               DO m = 1 , j

                  k = m + n - 1

                  DO il = 1 , nl

                     ZHC(il,k,m) = ZHC(il,k-1,m+1) - 2.0*(k-2)          &
                     & *ZHC(il,k-2,m)
                     ZHC(il,m,k) = ZHC(il,k,m)
                     ZHSo(il,k,m) = ZHSo(il,k-1,m+1) - 2.0*(k-2)        &
                     & *ZHSo(il,k-2,m)
                     ZHSo(il,m,k) = ZHSo(il,k,m)

                     IF ( m/=1 ) THEN

                        ZHC(il,k,m) = ZHC(il,k,m) + 2.0*(m-1)           &
                        & *ZHC(il,k-1,m-1)
                        ZHC(il,m,k) = ZHC(il,k,m)
                        ZHSo(il,k,m) = ZHSo(il,k,m) + 2.0*(m-1)         &
                        & *ZHSo(il,k-1,m-1)
                        ZHSo(il,m,k) = ZHSo(il,k,m)
                     ENDIF

                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         domm = 2*Nmax + 1
         ndol = 0
         ngor = 0

         DO dom = 1 , domm , 2
            DO iparr = 1 , Nstep

               ippar = (-1)**iparr

               CALL QUANTN(dom,iparr,Nmax,Homz,Hompr,nmain, &
               & nzmain,nrmain,lbmain,ismain,idim)

               DO isigp = 1 , 3 , 2

                  isig = isigp - 2
                  lamb = (dom-isig)/2

                  IF ( lamb<=Nmax ) THEN

                     lambp = lamb + 1

                     DO nprp = lambp , nmaxpn , 2

                        npr = nprp - 1
                        nro = (npr-lamb)/2
                        nrop = nro + 1

                        ipi = nrop + lamb

                        IF ( ipi>NDFACT ) STOP 'NDSFACT'

                        ynnro = real(sfact(nrop)/sfact(ipi))
                        ynnro = ynnro + ynnro

                        DO isigp1 = 1 , 3 , 2
                           isig1 = isigp1 - 2

                           IF ( isig1<=isig ) THEN

                              lamb1 = (dom-isig1)/2
                              laa = lamb - lamb1
                              lpha = (lamb+lamb1-1)/2 + 2

                              lphace = lamb + 2
                              lphaso = lamb + 1

                              IF ( lamb1<=Nmax ) THEN

                                 lambp1 = lamb1 + 1

                                 DO nprp1 = lambp1 , nmaxpn , 2

                                    npr1 = nprp1 - 1
                                    nro1 = (npr1-lamb1)/2
                                    nrop1 = nro1 + 1

                                    ipi1 = nrop1 + lamb1

                                    IF ( ipi1>NDFACT ) STOP 'NDSFACT'

                                    xnnro = real(sfact(nrop1)/sfact(ipi1))    &
                                    & *ynnro

                                    !         THE PRODUCT  OF  THE NORMALIZATION FACTORS

                                    xnrso1 = xnnro*vsop1
                                    xnrso2 = xnnro*vsop2

                                    !         INDEX  1  LABELS TERMS BEING NON-DIAGONAL IN LAMBDA

                                    xlam0 = FLOAT(lamb)
                                    xlam1 = FLOAT(lamb1)

                                    !          INTEGRATION OVER ETA

                                    DO il = 1 , nl

                                       eta = xl(il)
                                       fll = flpl(nrop,il,lambp)
                                       dfll = dflpl(nrop,il,lambp)

                                       IF ( laa/=0 ) THEN

                                          !         NON-DIAGONAL MATRIX ELEMENTS

                                          fll1 = flpl(nrop1,il,lambp1)
                                          dfll1 = dflpl(nrop1,il,lambp1)
                                          etaso = xxl(lpha,il)*xnrso1
                                          aoo = fll*fll1
                                          a1op = fll*dfll1
                                          a1po = fll1*dfll

                                          !         DFL  IS THE PRODUCT OF ETA  AND  DFL

                                          aop(il)                        &
                                          & = (xlam1*(aoo+a1op)-xlam0*   &
                                          & a1op)*etaso*wl(il)

                                          apo(il)                        &
                                          & = (xlam0*(aoo+a1po)-xlam1*   &
                                          & a1po)*etaso*wl(il)
                                       ELSE

                                          !         DIAGONAL MATRIX ELEMENTS

                                          fll1 = flpl(nrop1,il,lambp1)
                                          dfll1 = dflpl(nrop1,il,lambp1)
                                          etace = xxl(lphace,il)
                                          etaso = xxl(lphaso,il)
                                          aoo = fll1*fll
                                          a1po = fll1*dfll
                                          a1op = fll*dfll1
                                          apo(il) = -xlam0*(a1po+a1op)    &
                                          & *xnrso2*isig*etaso*wl(il)
                                          ao(il)                         &
                                          & = aoo*etace*0.5*xnnro*wl(il)
                                       ENDIF

                                    ENDDO


                                    DO np1 = nprp1 , nmaxpn

                                       n1 = np1 - 1

                                       IF ( Nstep/=1 ) THEN
                                          IF ( (-1)**n1/=ippar ) GOTO 4
                                       ENDIF

                                       nz1 = n1 - npr1

                                       ebase2 = Homz*(nz1+.5)           &
                                       & + Hompr*(npr1+1.)

                                       IF ( ebase2<=ECUtb ) THEN

                                          nho1 = npr1 + nz1
                                          nzp1 = nz1 + 1

                                          IF ( isig1<=0 ) THEN

                                             j = IMIN(nrop1,nzp1)
                                          ELSE

                                             j = IPLu(nrop1,nzp1)
                                          ENDIF

                                          DO np = nprp , nmaxpn
                                             n = np - 1

                                             IF ( Nstep/=1 ) THEN
                                                IF ( (-1)**n/=ippar ) GOTO 2
                                             ENDIF

                                             nz = n - npr

                                             ebase1 = Homz*(nz+.5)          &
                                             & + Hompr*(npr+1.)

                                             IF ( ebase1<=ECUtb ) THEN

                                                nho = npr + nz
                                                nzp = nz + 1

                                                IF ( isig<0 ) THEN
                                                   i = IMIN(nrop,nzp)
                                                ELSE
                                                   i = IPLu(nrop,nzp)
                                                ENDIF

                                                xso1 = 0.
                                                xso = 0.
                                                xco = 0.
                                                x = 0.

                                                IF ( j<=i ) THEN
                                                   IF ( IABS(n1-n)<=ILS ) THEN

                                                      IF ( laa/=0 ) THEN

                                                         !         INTEGRATION OF THE NON-DIAGONAL MATRIX ELEMENTS

                                                         DO il = 1 , nl

                                                            b1po = -0.5*ZHSo(il,nzp+1,nzp1)

                                                            IF ( nzp/=1 ) b1po = b1po +    &
                                                            & nz*ZHSo(il,nzp-1,nzp1)

                                                            b1op = -0.5*ZHSo(il,nzp,nzp1+1)

                                                            IF ( nzp1/=1 ) b1op = b1op +   &
                                                            & nz1*ZHSo(il,nzp,nzp1-1)

                                                            xso1 = xso1 +                  &
                                                            & (b1po*aop(il)+b1op*apo(il))

                                                         ENDDO

                                                         x = xso1
                                                      ELSE

                                                         !         INTEGRATION OF THE DIAGONAL MATRIX ELEMENTS


                                                         DO il = 1 , nl
                                                            xso = xso + apo(il)            &
                                                            & *ZHSo(il,nzp,nzp1)
                                                            xco = xco + ao(il)             &
                                                            & *ZHC(il,nzp,nzp1)
                                                         ENDDO

                                                         x = xso + xco
                                                      ENDIF
                                                   ENDIF

                                                   x = x*xnz(nzp)*xnz(nzp1)

                                                   !          KINETIC ENERGY OPERATOR

                                                   y = 1.0E-10
                                                   IF ( laa==0 ) THEN
                                                      IF ( nro==nro1 .AND. nz==nz1 ) &
                                                      & THEN
                                                         y = dhompr*FLOAT(npr+1)        &
                                                         & + dhomz*(nz+0.5)
                                                      ELSEIF ( IABS(nro-nro1)        &
                                                      & ==1 .AND. nz==nz1 ) THEN
                                                         ipo = nrop1 + lambp
                                                         ipop = nrop + lambp
                                                         IF ( nro>nro1 )                &
                                                         & y = real(dhompr*sq(nrop1+1)      &
                                                         & *sq(ipo))
                                                         IF ( nro1>nro )                &
                                                         & y = real(dhompr*sq(nrop+1)       &
                                                         & *sq(ipop))
                                                      ELSEIF ( nro==nro1 .AND.       &
                                                      & IABS(nz-nz1)==2 ) THEN
                                                         IF ( nz>nz1 )                  &
                                                         & y = real(-dhomz*sq(nzp1+1)       &
                                                         & *sq(nzp1+2)/2.)
                                                         IF ( nz1>nz )                  &
                                                         & y = real(-dhomz*sq(nzp+1)        &
                                                         & *sq(nzp+2)/2.)
                                                      ENDIF
                                                   ENDIF

                                                   !         S E T T I N G   T H E   H A M I L T O N I A N

                                                   suma = x + y
                                                   !
                                                   ELM(i,j) = suma
                                                   ELM(j,i) = suma

                                                   !          AUXILIARY PRINTS

                                                   numer = numer + 1

                                                   IF ( NPRint/=0 ) THEN
                                                      IF ( numer<=NPRint ) THEN
                                                         IF (((numer/28)*28)==numer .and. iwrite/=0)  &
                                                         & WRITE (12,99005)

99005                                                    FORMAT ('2',/,1X,              &
                                                         &'  I  J   LAMB SIG NRO  NZ   LAMP SGP NRP NZP    CENTRAL INTG     &
                                                         &S.O.-DIAG.     S.O.-NOND.     KIN.ENERGY'/)

                                                         anorm = xnz(nzp)*xnz(nzp1)

                                                         xso = xso*anorm
                                                         xco = xco*anorm

                                                         if(iwrite/=0) WRITE (12,99006) i , j , lamb ,&
                                                         & isig , nro , nz , lamb1 ,   &
                                                         & isig1 , nro1 , nz1 , xco ,  &
                                                         & xso , xso1 , y , n , n1 ,   &
                                                         & ELM(i,j) , ELM(j,i) , x
99006                                                    FORMAT (1X,2I3,3X,4I4,3X,4I4,  &
                                                         & 4X,4E16.9,/,1X,2I3,42X,     &
                                                         & 3E16.9)
                                                      ENDIF
                                                   ENDIF
                                                ENDIF
                                             ENDIF

                                             !          END OF THE MAIN LOOPS

2                                         ENDDO
                                       ENDIF
4                                   ENDDO
                                 ENDDO
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO

               i = NMAxi

               IF ( i/=0 ) THEN
                  IF ( i<=idim ) THEN

                     CALL QL(i,ELM,E,VVI,idim)

                     IF ( iquadr/=0 )                                     &
                     & CALL QQ20(i,lbmain,ismain,nrmain,nzmain,q20,  &
                     & idim,xiz,xiro,aqx,sq)

                     CALL ELIMIN(VVI,E,idim,i,dom,newi)

                     IF ( IDEcom==1 )                                   &
                     & CALL DECOMP(idim,i,nmain,nzmain,lbmain,   &
                     & dom,ancoef,azcoef,cparal,cantip,gfact,VVI,    &
                     & ippar,NSTep)

                     IF ( IDEcom==1 .AND. dom==domm .and. iwrite/=0) WRITE (12,99007)
99007                FORMAT ('1',///)

                     kres = newi
                     IF ( KOMple/=0 ) kres = i

                     IF ( JSIgn/=(-9) ) THEN

                        !===========================================================
                        !        *********   WRITE  WRITE  WRITE      *********
                        !===========================================================

                        if(iwrite/=0) WRITE (Nw,*) i , kres

                        IF ( kres/=0 .AND. Iarray==1 .and. iwrite/=0 ) WRITE (Nw,*)     &
                        & ((ELM(ilm1,ilm2),ilm1=1,i),ilm2=1,kres)

                        IF ( kres==0 .OR. Iarray/=1 .and. iwrite/=0) WRITE (Nw,*)      &
                        & NPArit
                     ENDIF

                     !          PRINTING OF THE FINAL RESULTS

                     ngor = ngor + i
                     nznak = 1

                     IF ( Nstep==1 ) nznak = 0

                     DO licz = 1 , i

                        liczk = ndol + licz

                        ENC(liczk) = VVI(licz)
                        qqq20(liczk) = q20(licz)
                        ILLim(liczk) = IFIX(E(licz))
                        KIJ(liczk) = (-1)**((iparr+1)*nznak)*dom

                     ENDDO

                     ndoll = ndol + 1
                     ndsub = ndol + kres
                     ndol = ndol + i

                     IF ( kres==0 ) THEN

                        !          *********   WRITE   *********

                        IF ( JSIgn/=(-9) .and. iwrite/=0 ) WRITE (Nw,*) NPArit

                        !===========================================================
                        !        *********   WRITE  WRITE  WRITE      *********
                        !===========================================================

                     ELSEIF ( IDEcom==1 ) THEN

                        IF ( JSIgn/=(-9) .and. iwrite/=0) WRITE (Nw,*)                 &
                        & (VVI(k),k=1,kres) , (KIJ(k),k=ndoll,ndsub) &
                        & , (q20(k),k=1,kres) ,                      &
                        & ((ancoef(k,kbas),k=1,kres),kbas=1,Nmax+1) ,&
                        & ((azcoef(k,kbas),k=1,kres),kbas=1,Nmax+1) ,&
                        & (cparal(k),k=1,kres) , (cantip(k),k=1,kres)&
                        & , (gfact(k),k=1,kres)

                     ELSE

                        IF ( JSIgn/=(-9) .and. iwrite/=0) WRITE (Nw,*)                 &
                        & (VVI(k),k=1,kres) , (KIJ(k),k=ndoll,ndsub)

                     ENDIF

                     IF ( JSIgn/=(-9) .and. iwrite/=0) WRITE (Nw,*) NPArit
                  ELSE

                     if(iwrite/=0) WRITE (12,99008) i
99008                FORMAT (//////,1X,                                 &
                     &'INCREASE THE DIMENSION OF THE HAMILTONIAN TO '&
                     & ,I3,'  AT LEAST')

                     STOP 'NDBLOK'
                  ENDIF
               ENDIF

            ENDDO
         ENDDO

         CALL ORD(ENC,KIJ,ILLim,iorder,ngor)

         if ( IPRot==0 ) SPL_N = enc ! store ENC array
         if ( IPRot==1 ) spl_P = enc ! store ENC array

         IF ( IPRot==0 .and. iwrite/=0) WRITE (12,99009)

99009    FORMAT (1X,129('*'),/,1X,'*',127(' '),'*',/,1X,'*',127(' '),   &
         & '*',/,1X,'*',8X,'SINGLE PARTICLE NUCLEAR SPECTRUM - ', &
         &'SHELL MODEL DEFORMED POTENTIAL WELL',35X,             &
         &'NEUTRONS      ','*',/,1X,'*',127X,'*',/,1X,'*',7X,    &
         &4(' NO)',3X,'  ENERGY',3X,' STATE ',5X),'*',/,1X,'*',  &
         & 127X,'*')
         IF ( IPRot==1 .and. iwrite/=0) WRITE (12,99010)

99010    FORMAT (1X,129('*'),/,1X,'*',127(' '),'*',/,1X,'*',127(' '),   &
         & '*',/,1X,'*',8X,'SINGLE PARTICLE NUCLEAR SPECTRUM - ', &
         &'SHELL MODEL DEFORMED POTENTIAL WELL',36X,             &
         &'PROTONS      ','*',/,1X,'*',127X,'*',/,1X,'*',7X,     &
         &4(' NO)',3X,'  ENERGY',3X,' STATE ',5X),'*',/,1X,'*',  &
         & 127X,'*')

         kik = ngor/4

         IF ( kik>IPRint ) kik = IPRint

         DO l = 1 , kik

            l1 = l
            l2 = kik + l
            l3 = 2*kik + l
            l4 = 3*kik + l

            il1 = 2 + ISIGN(1,KIJ(l1))*nznak + ILLim(l1)*3
            il2 = 2 + ISIGN(1,KIJ(l2))*nznak + ILLim(l2)*3
            il3 = 2 + ISIGN(1,KIJ(l3))*nznak + ILLim(l3)*3
            il4 = 2 + ISIGN(1,KIJ(l4))*nznak + ILLim(l4)*3

            ial1 = IABS(KIJ(l1))
            ial2 = IABS(KIJ(l2))
            ial3 = IABS(KIJ(l3))
            ial4 = IABS(KIJ(l4))

            if(iwrite/=0) WRITE (12,99015) l1 , ENC(l1) , ial1 , iznak(il1) , l2 ,    &
            & ENC(l2) , ial2 , iznak(il2) , l3 , ENC(l3) &
            & , ial3 , iznak(il3) , l4 , ENC(l4) , ial4 ,&
            & iznak(il4)

         ENDDO

         if(iwrite/=0) WRITE (12,99016)

         IF ( NCUtl/=0 ) CALL CUTEL(XEKat,Nmax)

         !         PRINTING THE QUADRUPOLE MOMENT TABLE

         IF ( IPRot==0 .AND. iquadr/=0 .and. iwrite/=0) WRITE (12,99011)

99011    FORMAT ('1',///,1X,129('*'),/,1X,'*',127(' '),'*',/,1X,'*',    &
         &127(' '),'*',/,1X,'*',8X,                              &
         &'SINGLE PARTICLE QUADRUPOLE MOMENTS - ',               &
         &'SHELL MODEL DEFORMED POTENTIAL WELL',33X,             &
         &'NEUTRONS      ','*',/,1X,'*',127X,'*',/,1X,'*',127X,  &
         &'*',/,1X,129('*'),/,1X,'*',127X,'*',/,1X,'*',127X,'*', &
         & /,1X,'*',7X,4(' NO)',3X,'  MOMENT',3X,' STATE ',5X),   &
         & '*',/,1X,'*',127X,'*')
         IF ( IPRot==1 .AND. iquadr/=0 .and. iwrite/=0) WRITE (12,99012)

99012    FORMAT ('1',///,1X,129('*'),/,1X,'*',127(' '),'*',/,1X,'*',    &
         &127(' '),'*',/,1X,'*',8X,                              &
         &'SINGLE PARTICLE QUADRUPOLE MOMENTS - ',               &
         &'SHELL MODEL DEFORMED POTENTIAL WELL',34X,             &
         &'PROTONS      ','*',/,1X,'*',127X,'*',/,1X,'*',127X,   &
         & '*',/,1X,129('*'),/,1X,'*',127X,'*',/,1X,'*',127X,'*', &
         & /,1X,'*',7X,4(' NO)',3X,'  MOMENT',3X,' STATE ',5X),   &
         & '*',/,1X,'*',127X,'*')

         IF ( iquadr/=0 ) THEN

            DO l = 1 , kik

               l1 = l
               l2 = kik + l
               l3 = 2*kik + l
               l4 = 3*kik + l

               il1 = 2 + ISIGN(1,KIJ(l1))*nznak + ILLim(l1)*3
               il2 = 2 + ISIGN(1,KIJ(l2))*nznak + ILLim(l2)*3
               il3 = 2 + ISIGN(1,KIJ(l3))*nznak + ILLim(l3)*3
               il4 = 2 + ISIGN(1,KIJ(l4))*nznak + ILLim(l4)*3

               ial1 = IABS(KIJ(l1))
               ial2 = IABS(KIJ(l2))
               ial3 = IABS(KIJ(l3))
               ial4 = IABS(KIJ(l4))

               if(iwrite/=0) WRITE (12,99015) l1 , qqq20(iorder(l1)) , ial1 , &
               & iznak(il1) , l2 , qqq20(iorder(l2)) ,   &
               & ial2 , iznak(il2) , l3 ,                &
               & qqq20(iorder(l3)) , ial3 , iznak(il3) , &
               & l4 , qqq20(iorder(l4)) , ial4 ,         &
               & iznak(il4)

            ENDDO

            if(iwrite/=0) WRITE (12,99016)
         ENDIF
      ELSE
         if(iwrite/=0) WRITE (12,99013) ILPunh
99013    FORMAT (///,1X,'ILPUNH TOO LARGE =',I4,///)
         STOP 'NDHERM'
      ENDIF

99014 FORMAT ('2',/)

99015 FORMAT (1X,'*',7X,4(I3,')',3X,F8.4,3X,I2,'/2',1X,A2,5X),'*')

99016 FORMAT (1X,'*',127X,'*',/,1X,129('*'))
   END SUBROUTINE ELMAT

   SUBROUTINE ELIMIN(Vvi,E,Idim,I,Dom,Newi)
      IMPLICIT NONE
      REAL E, gran, Vvi
      INTEGER I, ibiez, Idim, k, kli, klo, labda,  Newi
      INTEGER Dom
      DIMENSION Vvi(Idim), E(Idim)

      DO k = 1 , I
         E(k) = 0.0001
      ENDDO

      labda = (Dom+1)/2

      IF ( labda/=0 ) THEN
         gran = XEKat(labda)
      ELSE
         gran = 0.0
      ENDIF

      DO klo = 1 , I

         ibiez = klo
         IF ( Vvi(klo)>=gran ) GOTO 100

      ENDDO

      Newi = I

      RETURN

100   Newi = ibiez - 1

      DO kli = ibiez , I
         E(kli) = 1.0001
      ENDDO
   END SUBROUTINE ELIMIN

   SUBROUTINE EXTREM(Vekat,Iczy,Nmax,Rob,Cf)
      IMPLICIT NONE
      REAL AMN, AMP, Cf, hchv, potcc
      REAL Rob, Vekat, vge, vgehe, vjump, xmasa, HBArc
      INTEGER i, Iczy, k, kil, kilko, kis, kol, l, Nmax, nmaxt
      DIMENSION Vekat(Ndxekt), Iczy(Ndxekt), Rob(Ndplot), Cf(Ndplot)
      real, parameter :: hbarx = 197.32891
      real, parameter :: axn = 939.5527
      real, parameter :: axp = 938.2592

      !         ROUGH ESTIMATE OF THE EFFECTIVE BARRIERS VS. LAMBDA

      AMN = axn
      AMP = axp

      HBArc = hbarx

      nmaxt = Nmax + 1

      IF ( nmaxt<=NDLIM ) THEN

         vjump = .05
         xmasa = AMN

         IF ( IPRot==1 ) xmasa = AMP

         hchv = HBArc**2/2./xmasa/(CMAle*ROC)**2
         vgehe = 2.*Vmaxx
         kilko = IFIX(1.8*Vmaxx/vjump)

         DO i = 1 , nmaxt
            Iczy(i) = 0
            Vekat(i) = 0.
            Cf(i) = hchv*(FLOAT(i*i)-0.25)
         ENDDO

         DO kis = 1 , kilko

            vgehe = vgehe - vjump
            vge = vgehe**2

            potcc = POTF(Umaxx,vgehe)

            IF ( IPRot==1 ) potcc = potcc + COUL(Umaxx,vgehe)

            DO l = 1 , nmaxt

               IF ( Iczy(l)/=1 ) THEN

                  Rob(l) = potcc + Cf(l)/vge

                  IF ( Rob(l)>Vekat(l) ) THEN

                     Vekat(l) = Rob(l)
                  ELSE

                     Iczy(l) = 1
                  ENDIF
               ENDIF

            ENDDO

            k = 0

            IF ( nmaxt<=NDLIM ) THEN

               DO kol = 1 , nmaxt
                  k = k + Iczy(kol)
               ENDDO

               IF ( k==nmaxt ) GOTO 99999
            ENDIF

         ENDDO

         DO kil = 1 , nmaxt
            IF ( Iczy(kil)==0 ) Vekat(kil) = 0.
         ENDDO
         GOTO 99999
      ENDIF

      if(iwrite/=0) WRITE (12,99001)
99001 FORMAT (///,'STOP FROM EXTREM',///)

      STOP 'EXTREM'
99999 END SUBROUTINE EXTREM

   SUBROUTINE CUTEL(Vekat,Nmax)
      IMPLICIT NONE
      REAL umxxx, Vekat, vmxxx
      INTEGER kli, Nmax, nmaxt, nprint0
      DIMENSION Vekat(NDLIM)

      nmaxt = Nmax + 1

      if(iwrite/=0) WRITE (12,99001)
99001 FORMAT ('1',//)

      if(iwrite/=0) WRITE (12,99002) (kli,Vekat(kli),kli=1,8)
99002 FORMAT (1X,'MAXIMA OF THE EFFECTIVE BARRIERS VS. LAMBDA',//,1X,   &
      & 8('V(L=',I1,')=',F6.3,3X))

      nprint0 = 31

      IF ( Nstep==1 ) THEN

         if(iwrite/=0) WRITE (12,99003)
99003    FORMAT (//,1X,'DIMENSIONS OF THE OMEGA SUBBLOCKS VS. OMEGA',/)

         if(iwrite/=0) WRITE (12,99010) (kli,NDImom(kli,1),kli=1,nprint0,2)

         if(iwrite/=0) WRITE (12,99004)
99004    FORMAT (//,1X,'CORREPONDING MAXIMUM VALUES OF THE N-SHELLS',/)
         if(iwrite/=0) WRITE (12,99011) (kli,NOSMAX(kli,1),kli=1,nprint0,2)

      ELSE

         if(iwrite/=0) WRITE (12,99005)
99005    FORMAT (//,1X,'DIMENSIONS OF THE OMEGA SUBBLOCKS VS. OMEGA',   &
         & ' (PARITY=+)',/)

         if(iwrite/=0) WRITE (12,99010) (kli,NDImom(kli,2),kli=1,nprint0,2)

         if(iwrite/=0) WRITE (12,99006)
99006    FORMAT (//,1X,'CORREPONDING MAXIMUM VALUES OF THE N-SHELLS',   &
         & ' (PARITY=+)',/)

         if(iwrite/=0) WRITE (12,99011) (kli,NOSMAX(kli,2),kli=1,nprint0,2)

         if(iwrite/=0) WRITE (12,99007)
99007    FORMAT (//,1X,'DIMENSIONS OF THE OMEGA SUBBLOCKS VS. OMEGA',   &
         & ' (PARITY=-)',/)

         if(iwrite/=0) WRITE (12,99010) (kli,NDImom(kli,1),kli=1,nprint0,2)

         if(iwrite/=0) WRITE (12,99008)
99008    FORMAT (//,1X,'CORREPONDING MAXIMUM VALUES OF THE N-SHELLS',   &
         & ' (PARITY=-)',/)

         if(iwrite/=0) WRITE (12,99011) (kli,NOSMAX(kli,1),kli=1,nprint0,2)

      ENDIF

      umxxx = Umaxx*CMAle*ROC
      vmxxx = Vmaxx*CMAle*ROC

      if(iwrite/=0) WRITE (12,99009) Umaxx , umxxx , Vmaxx , vmxxx
99009 FORMAT (/,1X,70('='),1X,'UMAXX=',F6.3,' (',F6.3,' IN Z)',2X,      &
      & 'VMAXX=',F6.3,' (',F6.3,' IN RO)')
99010 FORMAT (1X,8('DIM(',I2,'/2)=',I3,3X),/,1X,                        &
      & 8('DIM(',I2,'/2)=',I3,3X))
99011 FORMAT (1X,8('NMX(',I2,'/2)=',I3,3X),/,1X,                        &
      & 8('NMX(',I2,'/2)=',I3,3X))
   END SUBROUTINE CUTEL

   SUBROUTINE DECOMP(Idim,Iact,Nmain,Nzmain,Lambda,Idom,An,Anz,  &
   & Cparal,Cantip,Gfacto,Energy,Ip,Ns)
      IMPLICIT NONE
      REAL An, antip, Anz, Cantip, Cparal, Energy, Gfacto, gk,   &
      & gl, gs, gsl, paral, spin
      INTEGER i, Iact, idd, Idim, Idom, Ip, j, k, &
      & lamb, Lambda, Nmain, Ns, nzact, Nzmain
      DIMENSION An(Idim,NDMAUX), Anz(Idim,NDMAUX)
      DIMENSION Energy(Idim), Nmain(Idim), Lambda(Idim), Nzmain(Idim), &
      & Cparal(Idim), Cantip(Idim), Gfacto(Idim)

      DO i = 1 , Iact

         DO j = 1 , NDMAUX
            Anz(i,j) = 0.00000
            An(i,j) = 0.00000
         ENDDO
      ENDDO

      !         CALCULATING THE NZ-Q.NUMBER DISTRIBUTION

      DO i = 1 , Iact
         DO j = 1 , Iact

            nzact = Nzmain(j)

            Anz(i,nzact+1) = Anz(i,nzact+1) + Elm(j,i)**2

         ENDDO
      ENDDO

      gl = 1.000
      gs = 5.585

      IF ( Iprot/=1 ) THEN

         gl = 0.0
         gs = -3.826
      ENDIF


      gs = 0.7*gs
      gsl = gs - gl

      IF ( Energy(1)<=EDEcgo .AND. Energy(Iact)>=EDEcdo ) THEN

         IF ( Ns==1 .and. iwrite/=0 ) WRITE (12,99001) Idom


99001    FORMAT (1X,'*',127X,'*',/,1X,'*',6X,'OMEGA =',I2,'/2',110X,'*',&
         & /,1X,'*',127X,'*',/,1X,'*',6X,'ENERGY',5X,'PRLL',2X,   &
         &'ANTI',4X,' GK  ',1X,' N=0 ',1X,' N=1 ',1X,' N=2 ',1X, &
         &' N=3 ',1X,' N=4 ',1X,' N=5 ',1X,' N=6 ',1X,' N=7 ',1X,&
         &' N=8 ',1X,' N=9 ',1X,'N=10 ',1X,'N=11 ',1X,'N=12 ',1X,&
         &'N=13 ',1X,'N=14 ',1X,'*',/,1X,'*',127X,'*')
      ENDIF

      IF ( Energy(1)<=EDEcgo .AND. Energy(Iact)>=EDEcdo ) THEN

         IF ( Ns==2 .and. iwrite/=0) WRITE (12,99002) Idom , Ip

99002    FORMAT (1X,'*',127X,'*',/,1X,'*',6X,'OMEGA =',I2,'/2',         &
         &'   PARITY =',I2,97X,'*',/,1X,'*',127X,'*',/,1X,'*',6X,&
         &'ENERGY',5X,'PRLL',2X,'ANTI',4X,' GK  ',1X,' N=0 ',1X, &
         &' N=1 ',1X,' N=2 ',1X,' N=3 ',1X,' N=4 ',1X,' N=5 ',1X,&
         &' N=6 ',1X,' N=7 ',1X,' N=8 ',1X,' N=9 ',1X,'N=10 ',1X,&
         &'N=11 ',1X,'N=12 ',1X,'N=13 ',1X,'N=14 ',1X,'*',/,1X,  &
         &'*',127X,'*')
      ENDIF

      DO i = 1 , Iact

         paral = 0.0
         antip = 0.0

         DO j = 1 , Iact

            !         PARALLEL COMPONENTS HERE

            lamb = Lambda(j)
            idd = 2*lamb + 1

            IF ( idd/=Idom ) THEN
               !         ... AND ANTIPARALLEL HERE
               antip = antip + Elm(j,i)**2
            ELSE
               paral = paral + Elm(j,i)**2
            ENDIF

            !         MAIN-SHELL DECOMPOSITION

            An(i,Nmain(j)+1) = An(i,Nmain(j)+1) + Elm(j,i)**2

         ENDDO

         spin = 0.5*FLOAT(Idom)
         gk = gl*spin + 0.5*gsl*(paral-antip)
         gk = gk/spin

         Cparal(i) = paral
         Cantip(i) = antip
         Gfacto(i) = gk

         IF ( Energy(i)<EDEcgo .AND. Energy(i)>EDEcdo .and. iwrite/=0) WRITE (12,99003)&
         & Energy(i) , paral , antip , gk , (An(i,k),k=1,15)

99003    FORMAT (1X,'*',5X,F8.4,3X,2(F5.3,1X),2X,F5.2,1X,15(F5.3,1X),1X,&
         &'*')

         IF ( Energy(i)<EDEcgo .AND. Energy(i)>EDEcdo .and. iwrite/=0) WRITE (12,99004)&
         & (Anz(i,k),k=1,15)
99004    FORMAT (1X,'*',' NZ>>',31X,15(F5.3,1X),1X,'*')

      ENDDO
   END SUBROUTINE DECOMP

   SUBROUTINE QQ20(I,Lamb,Isig,Nro,Nz,Q20,Idim,Xiz,Xiro,Aqx,Sq)
      IMPLICIT NONE
      REAL aq, aq1, Q20, Xiro, Xiz
      INTEGER I, ialfa, Idim, isid, Isig, isio, Lamb, lamd, lamg, &
      & m, m1, n, Nro, nrod, nrog, Nz, nzd, nzg, nzz
      DOUBLE PRECISION Sq
      DIMENSION Lamb(Idim), Isig(Idim), Sq(1), Q20(1), Nro(Idim), Nz(Idim)
      real, allocatable :: Aqx(:,:)
      allocate(Aqx(Idim, Idim))

      !         THE (Q20) MATRIX ELEMENTS IN THE HARMONIC OSCILLATOR BASIS

      DO m = 1 , I

         nzd = Nz(m)
         nrod = Nro(m)
         lamd = Lamb(m)
         isid = Isig(m)

         DO n = m , I

            Aqx(m,n) = 0.

            nrog = Nro(n)
            nzg = Nz(n)
            lamg = Lamb(n)
            isio = Isig(n)

            IF ( isid==isio ) THEN

               nzz = nzg - nzd

               IF ( nrog/=nrod ) THEN

                  IF ( nzz==0 ) THEN
                     IF ( nrog-nrod==1 ) THEN

                        !        NRO'=NRO+1 , NZ'=NZ

                        Aqx(m,n) = real(Xiro*Sq(nrog+1)*Sq(nrog+1+lamg))
                     ELSEIF ( nrod-nrog==1 ) THEN

                        !        NRO'=NRO-1 , NZ'=NZ

                        Aqx(m,n) = real(Xiro*Sq(nrod+1)*Sq(nrod+lamd+1))

                     ENDIF
                  ENDIF
               ELSEIF ( nzz==0 ) THEN

                  !        NZ'=NZ , NRO'=NRO

                  Aqx(m,n) = Xiz*(nzd+nzd+1.) - Xiro*(nrod+nrod+lamd+1.)
               ELSEIF ( IABS(nzz)==2 ) THEN
                  IF ( nzz==2 ) THEN

                     !        NZ'=NZ+2 , NRO'=NRO

                     Aqx(m,n) = real(Xiz*Sq(nzg)*Sq(nzg+1))
                  ELSE

                     !        NZ'=NZ-2 , NRO'=NRO

                     Aqx(m,n) = real(Xiz*Sq(nzd+1)*Sq(nzd))
                  ENDIF
               ENDIF
            ENDIF

         ENDDO

      ENDDO
      ! THE MATRIX ELEMENTS OF THE (Q20) OPERATOR IN THE WOODS-SAXON BASIS
      DO ialfa = 1 , I
         aq = 0.
         DO m = 1 , I
            aq1 = 0.
            m1 = m + 1
            IF ( m1<=I ) THEN
               DO n = m1 , I
                  aq1 = aq1 + Aqx(m,n)*Elm(n,ialfa)
               ENDDO
            ENDIF
            aq = aq + aq1*Elm(m,ialfa)*2. + Aqx(m,m)*Elm(m,ialfa)**2
         ENDDO
         Q20(ialfa) = aq
      ENDDO
   END SUBROUTINE QQ20

   FUNCTION FFFF(BBeta2)
      IMPLICIT NONE
      REAL b045, BBeta2, diffus, f1, f2, FFFF
      DATA diffus, b045/0.115 , 0.45/

      IF ( INCrea/=1 .OR. BBeta2<0.325 ) THEN
         FFFF = 1.0
         RETURN
      ENDIF
      f1 = (BBeta2-b045)/diffus
      f2 = (0.325-b045)/diffus
      FFFF = 1.0 + FINcla*(TANH(f1)-TANH(f2))
   END FUNCTION FFFF

   FUNCTION GGGG(BBeta2)
      IMPLICIT NONE
      REAL b045, BBeta2, diffus, f1, f2, GGGG
      DATA diffus , b045/0.115 , 0.45/

      IF ( INCrea/=1 .OR. BBeta2<0.325 ) THEN

         GGGG = 1.0
         RETURN

      ENDIF

      f1 = (BBeta2-b045)/diffus
      f2 = (0.325-b045)/diffus

      GGGG = 1.0 - FINcro*(TANH(f1)-TANH(f2))
   END FUNCTION GGGG

   SUBROUTINE SHAPE(P,Vd,Vg,N,Nn)
      IMPLICIT NONE
      REAL am, delmm, delp, delta, P, skala, space, spaceu, &
      & spacev, try, ubiez, Vd, Vg
      INTEGER i, indd, indeks, indg, j, k, lin, linp, ln, m,  &
      & mm, N, Nn
      DIMENSION P(N), Vd(N), Vg(N)
      CHARACTER*1 hblank, hstar, linia(101)
      DATA hblank, hstar/' ' , '*'/

      if(iwrite/=0) WRITE (12,99002)

      linp = 0

      CALL ORDD(P,Vd,Vg,N)

      am = 0.

      DO i = 1 , N
         IF ( Vg(i)>=am ) am = Vg(i)
      ENDDO

      spaceu = P(N) - P(1)
      spacev = 2.0*am
      space = AMAX1(spaceu,spacev)

      IF ( spaceu<=spacev ) linp = IFIX(FLOAT(Nn)/2.*(1.-spaceu/space))
      linp = linp + 1

      DO i = 1 , linp
         if(iwrite/=0) WRITE (12,99003)
      ENDDO

      lin = Nn - 2*(linp-1)

      m = N
      delta = 1./lin
      delp = 0.5/lin
      skala = 110./space/1.405

      !         NORMAL PRINT - 82.

      DO ln = 1 , lin

         mm = lin + 1 - ln
         m = m + 1

         DO j = 1 , 101
            linia(j) = hblank
         ENDDO

50       m = m - 1
         try = (P(m)-P(1))/spaceu
         delmm = mm*delta

         IF ( try<(delmm+delp) .AND. try>=(delmm-delp) ) THEN

            indeks = IFIX(Vg(m)*skala+1.5)
            indg = 50 + indeks
            indd = 52 - indeks

            linia(indg) = hstar
            linia(indd) = hstar

            GOTO 50
         ELSE

            ubiez = P(1) + spaceu*delmm
            if(iwrite/=0) WRITE (12,99001) ubiez , linia
99001       FORMAT (1X,5X,F7.3,4X,9X,'S',101(A1),'S')

         ENDIF

      ENDDO

      k = lin + linp

      DO i = k , Nn
         if(iwrite/=0) WRITE (12,99003)
      ENDDO

      if(iwrite/=0) WRITE (12,99002)

99002 FORMAT (26X,103('S'))
99003 FORMAT (26X,'S',50(' '),'X',50(' '),'S')
   END SUBROUTINE SHAPE

   SUBROUTINE QG32(Xl,Xu,FCT,Y)
      IMPLICIT NONE
      DOUBLE PRECISION a, b, c, Xl, Xu, Y
      DOUBLE PRECISION , external :: FCT

      !         THIS SUBROUTINE CALCULATES THE INTEGRAL OF THE EXTERNAL
      !         FUNCTION FCT USING A  32-POINT GAUSS QUADRATURE FORMULA
      !         (FROM  THE IBM SCIENTIFIC ROUTINE PACKAGE / IBM SYSTEM)

      a = .5*(Xu+Xl)
      b = Xu - Xl

      c = .4986319309D0*b
      Y = .0035093050D0*(FCT(a+c)+FCT(a-c))
      c = .4928057558D0*b
      Y = Y + .0081371974D0*(FCT(a+c)+FCT(a-c))
      c = .4823811278D0*b
      Y = Y + .0126960327D0*(FCT(a+c)+FCT(a-c))
      c = .4674530380D0*b
      Y = Y + .0171369315D0*(FCT(a+c)+FCT(a-c))
      c = .4481605779D0*b
      Y = Y + .0214179490D0*(FCT(a+c)+FCT(a-c))
      c = .4246838069D0*b
      Y = Y + .0254990296D0*(FCT(a+c)+FCT(a-c))
      c = .3972418980D0*b
      Y = Y + .0293420467D0*(FCT(a+c)+FCT(a-c))
      c = .3660910594D0*b
      Y = Y + .0329111114D0*(FCT(a+c)+FCT(a-c))
      c = .3315221335D0*b
      Y = Y + .0361728971D0*(FCT(a+c)+FCT(a-c))
      c = .2938578786D0*b
      Y = Y + .0390969479D0*(FCT(a+c)+FCT(a-c))
      c = .2534499545D0*b
      Y = Y + .0416559621D0*(FCT(a+c)+FCT(a-c))
      c = .2106756381D0*b
      Y = Y + .0438260465D0*(FCT(a+c)+FCT(a-c))
      c = .1659343011D0*b
      Y = Y + .0455869393D0*(FCT(a+c)+FCT(a-c))
      c = .1196436811D0*b
      Y = Y + .0469221995D0*(FCT(a+c)+FCT(a-c))
      c = .0722359808D0*b
      Y = Y + .0478193600D0*(FCT(a+c)+FCT(a-c))
      c = .0241538328D0*b
      Y = b*(Y+.0482700443D0*(FCT(a+c)+FCT(a-c)))
   END SUBROUTINE QG32

   SUBROUTINE QG16(Xl,Xu,FCT,Y)
      IMPLICIT NONE
      DOUBLE PRECISION a, b, c, FCT, Xl, Xu, Y

      !         THIS SUBROUTINE CALCULATES THE INTEGRAL OF AN EXTERNAL
      !         FUNCTION FCT USING A 16-POINT GAUSS QUADRATURE FORMULA
      !         (FROM THE IBM SCIENTIFIC ROUTINE PACKAGE)

      a = 0.5*(Xu+Xl)
      b = Xu - Xl

      c = .4947004675D0*b
      Y = .1357622970D-1*(FCT(a+c)+FCT(a-c))
      c = .4722875115D0*b
      Y = Y + .3112676197D-1*(FCT(a+c)+FCT(a-c))
      c = .4328156012D0*b
      Y = Y + .4757925584D-1*(FCT(a+c)+FCT(a-c))
      c = .3777022042D0*b
      Y = Y + .6231448563D-1*(FCT(a+c)+FCT(a-c))
      c = .3089381222D0*b
      Y = Y + .7479799441D-1*(FCT(a+c)+FCT(a-c))
      c = .2290083888D0*b
      Y = Y + .8457825970D-1*(FCT(a+c)+FCT(a-c))
      c = .1408017754D0*b
      Y = Y + .9130170752D-1*(FCT(a+c)+FCT(a-c))
      c = .4750625492D-1*b
      Y = b*(Y+.9472530523D-1*(FCT(a+c)+FCT(a-c)))
   END SUBROUTINE QG16

   SUBROUTINE RTNI(X,F,Derf,FCT,Xst,Eps,Iend,Ier)
      IMPLICIT NONE
      REAL a, Derf, dx, Eps, F, tdl, tdlf, X, Xst
      INTEGER i, Iend, Ier
      external FCT

      !         THIS SUBROUTINE CALCULATES THE ZERO OF THE EQUATION
      !
      !                         F(X)=0
      !
      !         USING THE NEWTONS METHOD
      !
      !
      !        X     ROOT OF EQUATION F(X)=0
      !        F     FUNCTION VALUE AT  X
      !        DERF  VALUE OF DERIVATIVE AT  X
      !
      !        FCT   EXTERNAL SUBROUTINE FCT(X,F,DERF).  IT COMPUTES
      !              FUNCTION VALUE F  AND DERIVATIVE  DERF  AT  X
      !
      !        XST   INITIAL GUESS OF  X
      !
      !        EPS   THE UPPER BOUND OF THE ERROR OF  0
      !
      !        IEND  MAX. NUMBER OF ITERATIONS
      !
      !        IER   RESULTANT ERROR PARAMETER
      !
      !              IER=0 - NO ERRORS
      !                  1 - NO CONVERGENCE AFTER IEND ITERATIONS
      !                  2 - AT ANY ITERATION DERF WAS EQUAL TO 0.
      !

      !        PREPARE ITERATION

      Ier = 0
      X = Xst
      tdl = X

      CALL FCT(tdl,F,Derf)

      tdlf = 100.*Eps

      !         START ITERATION LOOP

      DO i = 1 , Iend

         IF ( F==0 ) GOTO 100

         !        EQUATION IS NOT SATISFIED BY X

         IF ( Derf==0 ) GOTO 200

         !        ITERATION IS POSSIBLE

         dx = F/Derf
         X = X - dx
         tdl = X

         CALL FCT(tdl,F,Derf)

         !        TEST ON SATISFACTORY ACCURACY

         tdl = Eps
         a = ABS(X)

         IF ( a>1. ) tdl = tdl*a

         IF ( ABS(dx)<=tdl ) THEN
            IF ( ABS(F)<=tdlf ) GOTO 100
         ENDIF

      ENDDO

      !         END OF ITERATION LOOP

      Ier = 1

      !        NO CONVERGENCE AFTER IEND ITERATION STEPS, ERROR RETURN.

100   RETURN

      !        ERROR RETURN IN CASE OF ZERO DIVISOR

200   Ier = 2
   END SUBROUTINE RTNI

   SUBROUTINE BISEC(A,B,X,F,Eps,Iter,Ier)
      IMPLICIT NONE
      REAL A, abfa, abfb, abff, B, Eps, F, fa, fb, fx, X, xa, xb
      INTEGER Ier, it, Iter

      !         SOLVES THE EQUATION F(X)=0 USING A BISECTION METHOD

      xa = A
      xb = B
      it = 0

      fa = F(xa)
      fb = F(xb)
      !
      abfa = ABS(fa)
      X = A

      IF ( abfa>=Eps ) THEN

         abfb = ABS(fb)
         X = B

         IF ( abfb>=Eps ) THEN
            IF ( fa*fb>0.0 ) THEN

               if(iwrite/=0) WRITE (12,99001)
99001          FORMAT (///,1X,'BISEC OUT OF RANGE',///)

               STOP 'BISEC'
            ELSE

10             it = it + 1
               X = 0.5*(xa+xb)
               fx = F(X)
               abff = ABS(fx)

               IF ( abff>=Eps ) THEN
                  IF ( it>Iter ) THEN

                     Ier = 10
                     GOTO 99999
                  ELSE
                     IF ( fa*fx<0 ) THEN

                        !         FUNCTION CHANGES SIGN IN THE LEFT SUBINTERVAL

                        xb = X
                     ELSEIF ( fa*fx==0 ) THEN
                        GOTO 100
                     ELSE

                        !         FUNCTION CHANGES SIGN IN THE RIGHT SUBINTERVAL

                        xa = X
                     ENDIF
                     GOTO 10
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
100   Ier = 0
      RETURN
99999 END SUBROUTINE BISEC

   SUBROUTINE CEL2(Res,Ak,A0,B0,Ier)
      IMPLICIT NONE
      DOUBLE PRECISION A0, aa, aari, Ak, an, ari, B0, geo, Res, w
      INTEGER Ier


      !         THIS SUBROUTINE CALCULATES THE ELLIPTIC INTEGRALS OF THE
      !         FIRST AND OF THE SECOND KIND.
      !
      !         (FROM THE IBM SCIENTIFIC SUBROUTINR PACKAGE)

      Ier = 0

      geo = 1. - Ak*Ak

      IF ( geo<0 ) THEN

         Ier = 1

         RETURN
      ELSEIF ( geo==0 ) THEN

         IF ( B0<0 ) THEN

            Res = -1.D33

            RETURN
         ELSEIF ( B0==0 ) THEN

            Res = A0

            RETURN
         ELSE

            Res = 1.D33

            RETURN
         ENDIF
      ELSE

         geo = DSQRT(geo)

         ari = 1.
         aa = A0
         an = A0 + B0
         w = B0

50       w = w + aa*geo

         w = w + w
         aa = an
         aari = ari
         ari = geo + ari
         an = w/ari + an

         IF ( aari-geo<=1.D-9*aari ) THEN

            Res = .785398163397448D0*an/ari
         ELSE

            geo = DSQRT(geo*aari)

            geo = geo + geo

            GOTO 50
         ENDIF
      ENDIF
   END SUBROUTINE CEL2

   SUBROUTINE QL(N,A,E,V,Idim)
      IMPLICIT NONE
      REAL A, aav, b, c, E, f, g, h, hh, p, r, s, V
      INTEGER i, i1, iia, icr, Idim, ii1, ii2, j, j1, k, l, &
      & m, m1, N
      DIMENSION A(Idim,Idim), E(Idim), V(Idim)

      !         MATRIX DIAGONALISATION ROUTINE   (FORTRAN VERSION OF THE
      !         ORIGINAL ALGOL PROGRAM, COMPUTER CENTRE OF P.A.N. WARSAW

      IF ( N>=2 ) THEN

         DO icr = 2 , N

            i = 2 + N - icr
            l = i - 2
            f = A(i,i-1)
            g = 0.
            s = 0.

            IF ( l>=1 ) THEN

               DO k = 1 , l

                  c = A(i,k)**2
                  b = g + c
                  s = g - b + c + s
                  g = b

               ENDDO
            ENDIF

            g = g + s
            h = g + f*f

            IF ( g<=1.0D-38 ) THEN

               E(i) = f

               h = 0.
            ELSE

               l = l + 1

               IF ( f<0 ) THEN

                  E(i) = SQRT(h)
               ELSE

                  E(i) = -SQRT(h)
               ENDIF

               g = E(i)
               h = h - f*g
               A(i,i-1) = f - g
               f = 0.
               p = 0.

               IF ( l>=1 ) THEN

                  DO j = 1 , l

                     r = A(i,j)/h
                     A(j,i) = r
                     g = 0.
                     s = 0.

                     IF ( j>=1 ) THEN

                        DO k = 1 , j
                           c = A(j,k)*A(i,k)
                           b = g + c
                           s = g - b + c + s
                           g = b
                        ENDDO
                     ENDIF

                     j1 = j + 1

                     IF ( l>=j1 ) THEN

                        DO k = j1 , l
                           c = A(k,j)*A(i,k)
                           b = g + c
                           s = g - b + c + s
                           g = b
                        ENDDO
                     ENDIF

                     g = g + s
                     E(j) = g/h
                     c = g*r
                     b = f + c
                     p = f - b + c + p
                     f = b

                  ENDDO
               ENDIF

               hh = (f+p)/(2.*h)

               IF ( l>=1 ) THEN

                  DO j = 1 , l
                     !
                     f = A(i,j)
                     g = E(j) - hh*f
                     E(j) = g

                     IF ( j<1 ) GOTO 20

                     DO k = 1 , j
                        A(j,k) = A(j,k) - f*E(k) - g*A(i,k)
                     ENDDO

                  ENDDO
               ENDIF
            ENDIF

20          V(i) = h

         ENDDO
      ENDIF

      V(1) = 0.

      IF ( N>=1 ) THEN

         DO i = 1 , N

            l = i - 1

            IF ( V(i)/=0 ) THEN

               IF ( l>=1 ) THEN

                  DO j = 1 , l

                     g = 0.
                     s = 0.

                     IF ( l>=1 ) THEN

                        DO k = 1 , l
                           c = A(i,k)*A(k,j)
                           b = g + c
                           s = g - b + c + s
                           g = b
                        ENDDO
                     ENDIF

                     IF ( l<1 ) GOTO 40

                     DO k = 1 , l
                        A(k,j) = A(k,j) - g*A(k,i)
                     ENDDO

                  ENDDO
               ENDIF
            ENDIF

40          V(i) = A(i,i)
            A(i,i) = 1.

            IF ( l>=1 ) THEN

               DO j = 1 , l
                  A(i,j) = 0.
                  A(j,i) = 0.
               ENDDO
            ENDIF

         ENDDO
      ENDIF

      IF ( N>=2 ) THEN

         DO i = 2 , N
            ii1 = i - 1
            E(ii1) = E(i)
         ENDDO
      ENDIF

      E(N) = 0.

      k = N - 1

      IF ( N>=1 ) THEN

         DO l = 1 , N

            j = 0

60          IF ( k<l ) THEN

               m = N
            ELSE

               m = l

               IF ( ABS(E(m))>4.0E-9*(ABS(V(m))+ABS(V(m+1))) ) m = N
            ENDIF
            p = V(l)

            IF ( m/=l ) THEN

               IF ( j==30 ) GOTO 100

               j = j + 1

               g = (V(l+1)-p)/(2.0*E(l))
               r = SQRT(1.0+g**2)

               IF ( g<0 ) THEN

                  g = V(m) - p + E(l)/(g-r)
               ELSE

                  g = V(m) - p + E(l)/(g+r)
               ENDIF
               s = 1.

               c = 1.
               p = V(m)
               m1 = m - 1

               IF ( m1>=l ) THEN

                  DO icr = l , m1

                     i = l + m1 - icr
                     h = E(i)
                     f = s*h
                     b = c*h

                     IF ( ABS(f)<ABS(g) ) THEN
                        c = f/g

                        r = SQRT(c**2+1.)
                        ii1 = i + 1
                        E(ii1) = g*r
                        s = c/r
                        c = 1./r
                     ELSE

                        c = g/f

                        r = SQRT(c**2+1)
                        ii1 = i + 1
                        E(ii1) = f*r
                        s = 1./r

                        c = c/r
                     ENDIF

                     h = V(i)

                     f = c*h - s*b
                     g = c*b - s*p
                     r = h + p
                     p = c*f - s*g
                     g = s*f + c*g
                     ii2 = i + 1
                     V(ii2) = r - p

                     IF ( N>=1 ) THEN

                        DO iia = 1 , N
                           f = A(iia,i+1)
                           h = A(iia,i)
                           A(iia,i+1) = s*h + c*f
                           A(iia,i) = c*h - s*f
                        ENDDO
                     ENDIF

                  ENDDO
               ENDIF

               V(l) = p
               E(l) = g
               E(m) = 0.

               GOTO 60
            ENDIF

         ENDDO
      ENDIF

      IF ( N>=1 ) THEN

         DO i = 1 , N

            k = i
            p = V(i)
            i1 = i + 1

            IF ( N>=i1 ) THEN

               DO j = i1 , N

                  IF ( V(j)<p ) THEN

                     k = j

                     p = V(j)
                  ENDIF

               ENDDO
            ENDIF

            IF ( k/=i ) THEN

               aav = V(i)
               V(k) = aav
               V(i) = p

               IF ( N>=1 ) THEN

                  DO j = 1 , N
                     p = A(j,i)
                     A(j,i) = A(j,k)
                     A(j,k) = p
                  ENDDO
               ENDIF
            ENDIF

         ENDDO
      ENDIF
      GOTO 99999

100   if(iwrite/=0) WRITE (12,99001) N
99001 FORMAT (//,1X,'NO CONVERGENCE IN QL  N =',I3,///)
99999 END SUBROUTINE QL

   SUBROUTINE DHEP(X,N,Her,Dher,Ndim)
      IMPLICIT NONE
      REAL df, Dher, f, Her, X
      INTEGER i, N, Ndim
      DIMENSION Her(Ndim), Dher(Ndim)

      !         DHEP CALCULATES HERMITE POLYNOMIALS AT POINT X UP TO N-TH
      !         DEGREE. IT ALSO CALCULATES  DERIVATIVES  DHER(I) OF THESE
      !         POLYNOMIALS.
      !                            DHER(I)=2*I*HER(I-1)

      Her(1) = 1.
      Dher(1) = 0.

      IF ( N>0 ) THEN

         Her(2) = X + X
         Dher(2) = 2.

         IF ( N>1 ) THEN

            DO i = 2 , N

               f = X*Her(i) - FLOAT(i-1)*Her(i-1)
               Her(i+1) = f + f
               df = FLOAT(i)*Her(i)
               Dher(i+1) = df + df

            ENDDO
            GOTO 99999
         ENDIF
      ENDIF
      RETURN
99999 END SUBROUTINE DHEP

   SUBROUTINE HRECUR(Pn,Dpn,Pn1,X,Nn)
      IMPLICIT NONE
      DOUBLE PRECISION dp, dp1, Dpn, dq, fj, fj2, p, p1, Pn, &
      & Pn1, q, X
      INTEGER j, Nn

      !         AUXILIARY ROUTINE FROM THE BOOK BY STRODE

      p1 = 1.
      p = X
      dp1 = 0.
      dp = 1.

      DO j = 2 , Nn
         fj = j
         fj2 = (fj-1.)/2.
         q = X*p - fj2*p1
         dq = X*dp + p - fj2*dp1
         p1 = p
         p = q
         dp1 = dp
         dp = dq
      ENDDO
      Pn = p
      Dpn = dp
      Pn1 = p1
   END SUBROUTINE HRECUR

   SUBROUTINE HROOT(X,Nn,Dpn,Pn1,Eps)
      IMPLICIT NONE
      DOUBLE PRECISION d , dp , Dpn , Eps , p , Pn1 , X
      INTEGER iter , Nn

      !        IMPROVES THE APPROXIME ROOT X  IN ADDITION WE ALSO OBTAIN
      !                                         PN1=VALUE OF H(N-1) AT X
      !        FROM THE MONOGRAPH BY STROUD

      iter = 0
100   iter = iter + 1

      CALL HRECUR(p,dp,Pn1,X,Nn)

      d = p/dp
      X = X - d

      IF ( DABS(d)<=Eps ) THEN
         Dpn = dp
      ELSE
         IF ( iter<10 ) GOTO 100
         Dpn = dp
      ENDIF
   END SUBROUTINE HROOT

   SUBROUTINE HERMIT(Nn,X,A,Eps,Fact)
      IMPLICIT NONE
      REAL A, dpn, Eps, pn1, s, X, xt
      INTEGER i, n1, n2, ni, Nn
      DOUBLE PRECISION xtt, Fact, cc, eeps, ddpn, ppn1
      INTEGER fn
      DIMENSION X(Ndherm), A(Ndherm), Fact(1)

      !         HERMIT CALCULATES ZEROES X(I) OF THE N-TH ORDER HERMITE
      !         POLYNOMIAL.   THE  LARGEST  ZERO  IS  STORED  IN  X(1).
      !         IT CALCULATES  ALSO THE CORRESPONDING COEFFICIENTS A(I)
      !         OF THE  NN-TH ORDER  GAUSS-HERMITE  QUADRATURE  FORMULA
      !         OF THE DEGREE 2*NN-1 .
      !
      !         (TRIVIAL MODIFICATIONS OF THE ROUTINE FROM STROUD BOOK)

      eeps = Eps

      fn = Nn
      n1 = Nn - 1
      n2 = (Nn+1)/2

      cc = 1.7724538509D0*Fact(fn)/(2.**n1)
      s = (2.*fn+1)**.16667

      DO i = 1 , n2

         IF ( i<1 ) GOTO 100
         IF ( i==1 ) THEN

            !        LARGEST ZERO

            xt = s**3 - 1.85575/s
         ELSE

            IF ( i<2 ) GOTO 100
            IF ( i==2 ) THEN

               !        SECOND ZERO

               xt = xt - 1.14*fn**.426/xt
            ELSE

               IF ( i<3 ) GOTO 100
               IF ( i==3 ) THEN

                  !        THIRD ZERO

                  xt = 1.86*xt - .86*X(1)
               ELSE

                  IF ( i<4 ) GOTO 100
                  IF ( i==4 ) THEN

                     !        FOURTH ZERO

                     xt = 1.91*xt - .91*X(2)
                  ELSE

                     !        ALL OTHER ZEROS

                     xt = 2.*xt - X(i-2)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         xtt = xt

         CALL HROOT(xtt,Nn,ddpn,ppn1,eeps)

         dpn = SNGL(ddpn)
         pn1 = SNGL(ppn1)
         xt = SNGL(xtt)

         X(i) = xt
         A(i) = real(cc/dpn/pn1)
         ni = Nn - i + 1
         X(ni) = -xt
         A(ni) = A(i)

100   ENDDO
   END SUBROUTINE HERMIT

   SUBROUTINE DLAP(X,Xlag,Dxlag,N,Alpha,Ndim)
      IMPLICIT NONE
      REAL Dxlag, f, X, Xlag
      INTEGER i, N, Ndim
      INTEGER Alpha
      DIMENSION Xlag(Ndim), Dxlag(Ndim)

      !         DLAP CALCULATES LAGUERRE POLYNOMIALS AND THEIR DERIVATIVES
      !         UP THE THE N-TH ORDER

      Xlag(1) = 1.
      Dxlag(1) = 0.

      IF ( N>0 ) THEN

         Xlag(2) = 1. + Alpha - X
         Dxlag(2) = -X

         IF ( N>1 ) THEN

            DO i = 2 , N

               !         RECURRENCE DEFINITIONS

               f = (FLOAT(i+i+Alpha-1)-X)*Xlag(i) - FLOAT(i-1+Alpha)    &
               & *Xlag(i-1)
               Xlag(i+1) = f/FLOAT(i)
               Dxlag(i+1) = FLOAT(i)*Xlag(i+1) - FLOAT(i+Alpha)*Xlag(i)

            ENDDO
            GOTO 99999
         ENDIF
      ENDIF

      RETURN
99999 END SUBROUTINE DLAP

   SUBROUTINE LGROOT(X,Nn,Alf,Dpn,Pn1,B,C,Eps)
      IMPLICIT NONE
      REAL B, C
      INTEGER iter, Nn
      DOUBLE PRECISION d, p, dp, X, Dpn, Pn1, Alf, Eps
      DIMENSION B(Nn), C(Nn)

      !        IMPROVES  APPROXIMATE ROOT  X
      !        IN ADDITION WE ALSO OBTAIN
      !
      !                  DPN = DERIVATIVE OF P(N) AT X
      !                  PNI = VALUE OF P(N-1) AT  X

      iter = 0
100   iter = iter + 1

      CALL LGRECR(p,dp,Pn1,X,Nn,Alf,B,C)

      d = p/dp
      X = X - d

      IF ( DABS(d/X)<=Eps ) THEN

         Dpn = dp
      ELSE

         IF ( iter<10 ) GOTO 100
         Dpn = dp
      ENDIF
   END SUBROUTINE LGROOT

   SUBROUTINE LGRECR(Pn,Dpn,Pn1,X,Nn,Dalf,B,C)
      IMPLICIT NONE
      REAL B, C
      DOUBLE PRECISION Dalf, dp, dp1, Dpn, dq, p, p1, Pn, Pn1, q , X
      INTEGER j, Nn
      DIMENSION B(Nn), C(Nn)

      !         AUXILIARY ROUTINE - RECURRENCES (CALLED FROM LGROOT)

      p1 = 1.
      p = X - Dalf - 1.
      dp1 = 0.
      dp = 1.

      DO j = 2 , Nn
         q = (X-B(j))*p - C(j)*p1
         dq = (X-B(j))*dp + p - C(j)*dp1
         p1 = p
         p = q
         dp1 = dp
         dp = dq
      ENDDO

      Pn = p
      Dpn = dp
      Pn1 = p1

   END SUBROUTINE LGRECR

   SUBROUTINE LAGUER(Nn,X,A,Alf,B,C,Fact,Eps)
      IMPLICIT NONE
      REAL A, Alf, B, C, csa, csx, Eps, fi, fn, r1, r2, &
      & ratio, tsx, X, xt
      INTEGER i, ialfp, j,  Nn
      DOUBLE PRECISION cc, Fact, xtt, aalf, dpn, pn1, eeps
      DIMENSION X(Ndlagr), A(Ndlagr), B(Ndlagr), C(Ndlagr), Fact(1)

      !         LAGUER  CALCULATES  ZEROS AND WEIGHTS OF THE GAUSS-LAGUERRE
      !         INTEGRATION FORMULA
      !
      !         (TRIVIAL MODOFICATIONS OF THE ROUTINES FROM THE STROUD BOOK)

      fn = FLOAT(Nn)

      csx = 0.
      csa = 0.

      ialfp = IFIX(Alf) + 1
      cc = Fact(ialfp)
      tsx = fn*(fn+Alf)

      DO j = 2 , Nn
         cc = cc*C(j)
      ENDDO

      DO i = 1 , Nn

         IF ( i<1 ) THEN
         ELSEIF ( i==1 ) THEN

            !            SMALLEST ZERO

            xt = (1.+Alf)*(3.+.92*Alf)/(1.+2.4*fn+1.8*Alf)

         ELSEIF ( i<2 ) THEN
         ELSEIF ( i==2 ) THEN

            !            SECOND ZERO

            xt = xt + (15.+6.25*Alf)/(1.+.9*Alf+2.5*fn)
         ELSE

            !            ALL OTHER ZEROS

            fi = FLOAT(i-2)
            r1 = (1.+2.55*fi)/(1.9*fi)
            r2 = 1.26*fi*Alf/(1.+3.5*fi)
            ratio = (r1+r2)/(1.+.3*Alf)
            xt = xt + ratio*(xt-X(i-2))
         ENDIF

         xtt = xt
         aalf = Alf
         eeps = Eps

         CALL LGROOT(xtt,Nn,aalf,dpn,pn1,B,C,eeps)

         xt = real(xtt)

         X(i) = xt
         A(i) = real(cc/dpn/pn1)
         csx = csx + xt
         csa = csa + A(i)
      ENDDO

   END SUBROUTINE LAGUER

   SUBROUTINE WYKRES(A,B,C,IIa,Icontr,Ilo,Wmin,Wmax,Ibind)
      IMPLICIT NONE
      REAL A, aa, amax, amin, B, bb, bmax, bmin, C, cc, cmax,   &
      & cmin, da, diff, diffe, skala, skale, Wmax, Wmin, xmax, xmin
      INTEGER i, IIa, Ibind, Icontr, Ilo, ilog, ind, j, m1, m2, m3
      DIMENSION A(IIa), B(IIa), C(IIa), skala(11)
      CHARACTER*1 hbl, hx, ho, hplu, hsta, hdot, linia(100)
      DATA hbl, hx, ho, hplu, hsta, hdot/' ' , 'X' , 'O' , '+' ,   &
      &'*' , '.'/

      !         A,B,C   -  ARRAYS TO BE PLOTTED
      !         ICONTR  -  ACTUAL NUMBER OF THE ARRAYS (1.OR.2.OR.3)
      !         ILO     -  SCALE INDEX  (.EQ.0  -  LINEAR)
      !                                 (.EQ.1  -  LOGARYTHMIC)
      !                                 (.EQ.2  -  CHOOSE THE SCALE)
      !         IBIND   -  IF (.NE.0)  USE THE NORMALIZED SCALE BEING DEFINED
      !                                BY  BMIN  AND  BMAX

      ilog = Ilo
      amin = +1.0E+33
      bmin = +1.0E+33
      cmin = +1.0E+33
      amax = -1.0E+33
      bmax = -1.0E+33
      cmax = -1.0E+33
      ind = Icontr - 2
      IF ( ind<0 ) GOTO 100
      IF ( ind/=0 ) THEN
         DO i = 1 , IIa
            aa = A(i)
            IF ( amin>aa ) amin = aa
            IF ( amax<aa ) amax = aa
         ENDDO
      ENDIF
      DO i = 1 , IIa
         bb = B(i)
         IF ( bmin>bb ) bmin = bb
         IF ( bmax<bb ) bmax = bb
      ENDDO
100   DO i = 1 , IIa
         cc = C(i)
         IF ( cmin>cc ) cmin = cc
         IF ( cmax<cc ) cmax = cc
      ENDDO
      IF ( ind<0 ) THEN
         xmin = cmin
         xmax = cmax
      ELSEIF ( ind==0 ) THEN
         xmin = AMIN1(bmin,cmin)
         xmax = AMAX1(bmax,cmax)
      ELSE
         xmin = AMIN1(amin,bmin,cmin)
         xmax = AMAX1(amax,bmax,cmax)
      ENDIF
      IF ( Ibind/=0 ) THEN
         diffe = Wmax - Wmin
         IF ( diffe>=(xmax-xmin) ) THEN
            xmin = Wmin
            xmax = Wmax
         ENDIF
      ENDIF
      diff = xmax - xmin
      IF ( ilog==2 ) ilog = 0
      IF ( diff>9999. ) ilog = 1
      da = diff/100.
      IF ( ilog/=1 ) THEN
         DO j = 1 , 11
            skala(j) = xmin + (j-1)*da*10.
         ENDDO
         if(iwrite/=0) WRITE (12,99001) skala
99001    FORMAT (/,29X,102('S'),//,23X,11(2X,F8.2))
         if(iwrite/=0) WRITE (12,99002)
99002    FORMAT (3X,'A .....',2X,'B +++++',2X,'C *****',1X,102('S'))
         DO i = 1 , IIa
            DO j = 1 , 100
               linia(j) = hbl
            ENDDO
            m2 = 1
            m3 = 1
            m1 = IFIX((C(i)-xmin)/da)
            IF ( m1>100 ) m1 = 100
            IF ( m1<1 ) m1 = 1
            IF ( ind/=-1 ) THEN
               m2 = IFIX((B(i)-xmin)/da)
               IF ( m2>100 ) m2 = 100
               IF ( m2<1 ) m2 = 1
               IF ( ind/=0 ) THEN
                  m3 = IFIX((A(i)-xmin)/da)
                  IF ( m3>100 ) m3 = 100
                  IF ( m3<1 ) m3 = 1
               ENDIF
            ENDIF
            linia(m1) = hsta
            linia(m2) = hplu
            linia(m3) = hdot
            IF ( m1==m2 .OR. m1==m3 ) linia(m1) = hx
            IF ( m2==m3 ) linia(m2) = hx
            IF ( m1==m2 .AND. m1==m3 ) linia(m1) = ho
            if(iwrite/=0) WRITE (12,99003) A(i) , B(i) , C(i) , linia
99003       FORMAT (1X,3(1X,F8.2),' S',100A1,'S')
         ENDDO
         if(iwrite/=0) WRITE (12,99004)
99004    FORMAT (29X,102('S'))
      ELSE
         skale = ALOG10(diff+0.0001)
         da = skale/100.
         DO j = 1 , 11
            skala(j) = (j-1)*da*10.
         ENDDO
         if(iwrite/=0) WRITE (12,99005) (skala(j),j=1,10)
99005    FORMAT (/,32X,102('S'),//,29X,10(1X,E9.1),'    .')
         if(iwrite/=0) WRITE (12,99006)
99006    FORMAT (4X,'A .....',3X,'B +++++',3X,'C *****',1X,102('S'))
         DO i = 1 , IIa
            DO j = 1 , 100
               linia(j) = hbl
            ENDDO
            m2 = 1
            m3 = 1
            m1 = IFIX(ALOG10(C(i)-xmin+0.000001)/da)
            IF ( m1>100 ) m1 = 100
            IF ( m1<1 ) m1 = 1
            IF ( ind/=-1 ) THEN
               m2 = IFIX(ALOG10(B(i)-xmin+0.000001)/da)
               IF ( m2>100 ) m2 = 100
               IF ( m2<1 ) m2 = 1
               IF ( ind/=0 ) THEN
                  m3 = IFIX(ALOG10(A(i)-xmin+0.000001)/da)
                  IF ( m3>100 ) m3 = 100
                  IF ( m3<1 ) m3 = 1
               ENDIF
            ENDIF
            linia(m1) = hsta
            linia(m2) = hplu
            linia(m3) = hdot
            IF ( m1==m2 .OR. m1==m3 ) linia(m1) = hx
            IF ( m2==m3 ) linia(m2) = hx
            IF ( m1==m2 .AND. m1==m3 ) linia(m1) = ho
            if(iwrite/=0) WRITE (12,99007) A(i) , B(i) , C(i) , linia
99007       FORMAT (1X,3(1X,E9.1),' S',100A1,'S')
         ENDDO
         if(iwrite/=0) WRITE (12,99008)
99008    FORMAT (32X,102('S'))
      ENDIF
      if(iwrite/=0) WRITE (12,99009) amin , bmin , cmin , amax , bmax , cmax , xmin , &
      & xmax
99009 FORMAT (/,1X,'AMIN=',E12.5,'  BMIN=',E12.5,'  CMIN=',E12.5,/,1X,  &
      &'AMAX=',E12.5,'  BMAX=',E12.5,'  CMAX=',E12.5,/,1X,       &
      & 'XMIN=',E12.5,'  XMAX=',E12.5)
      IF ( ilog==1 .and. iwrite/=0) WRITE (12,99010)
99010 FORMAT (/,1X,'LOGARYTHMIC  SCALE  PRINTING')
      IF ( ilog/=1 .and. iwrite/=0) WRITE (12,99011)
99011 FORMAT (/,1X,'LINEAR  SCALE  PRINTING')
   END SUBROUTINE WYKRES

   SUBROUTINE ORD(A,K,L,M,N)
      IMPLICIT NONE
      REAL A, temp
      INTEGER i, ip1, j, K, L, lw, M, mw, N, n1, nemp
      DIMENSION A(N), K(N), M(N), L(N)

      !        ORDERING A SERIES OF REAL NUMBERS "A' IN INCREASING ORDER

      n1 = N - 1

      DO i = 1 , n1

         ip1 = i + 1

         DO j = ip1 , N

            IF ( A(i)>A(j) ) THEN

               temp = A(i)

               A(i) = A(j)
               A(j) = temp

               nemp = K(i)
               K(i) = K(j)
               K(j) = nemp

               lw = L(i)
               L(i) = L(j)
               L(j) = lw

               mw = M(i)
               M(i) = M(j)
               M(j) = mw
            ENDIF

         ENDDO
      ENDDO
   END SUBROUTINE ORD

   SUBROUTINE ORDD(A,P,Q,N)
      IMPLICIT NONE
      REAL A, P, pemp, Q, qw, temp
      INTEGER i, ip1, j, N, n1
      DIMENSION A(N), P(N), Q(N)

      !         SIMILAR TO ORD (CF. PRECEDING ROUTINE)

      n1 = N - 1

      DO i = 1 , n1

         ip1 = i + 1

         DO j = ip1 , N

            IF ( A(i)>A(j) ) THEN

               temp = A(i)
               A(i) = A(j)
               A(j) = temp

               pemp = P(i)
               P(i) = P(j)
               P(j) = pemp

               qw = Q(i)
               Q(i) = Q(j)
               Q(j) = qw
            ENDIF

         ENDDO
      ENDDO
   END SUBROUTINE ORDD

   SUBROUTINE STRUT(ma,mn,ENC,ESM)
      implicit none
      integer, parameter :: NBCS=272
      integer :: ITER
      integer :: i, j, k, ma, mn
      real(r8) :: HW0, GAMMA, ESP, ASN, EF
      real(r8) :: DN, AN, EN, X, PHI, X2, EX, EXG
      real(r8) :: ANY, ENY, DNY, ESR, AS, DNN, DF, EFST
      real(r8) :: SIPI=0.5641896_r8  !1/SQRT(PI)
      real(r8) :: ESM
      real(r8) :: E(NBCS), ENC(NBCS)

      DO J = 1, NBCS
         E(J) = ENC(J)
      ENDDO

      HW0 = 41./real(ma)**0.3333

      K = NBCS
      GAMMA = HW0*1.2
      ESP = 0.d0
      DO I = 1,mn/2
         ESP = ESP + E(I)
      END DO
      ESP = ESP*2.
      IF(MOD(mn,2).EQ.1)ESP = ESP + E(mn/2+1)
      ASN = real(mn)/2
      EF = 0.5*(E(mn/2)+E(mn/2+1))

      ITER=0
      do while (.true.)
         ITER = ITER+1
         DN = 0.
         AN = 0.
         EN = 0.
         DO I = 1, K
            X = (EF-E(I))/GAMMA
            PHI = ERF(X)
            X2 = X*X
            EX = SIPI*EXP(-X2)
            EXG = EX/GAMMA
            ANY = (1.+PHI)/2.+(1./12.*X2*X2-2./3.*X2+19./16.)*X*EX
            ENY = (1./12.*X2*X2*X2-5./8.*X2*X2+15./16.*X2-5./32.)* &
            &     EX*GAMMA+E(I)*ANY
            DNY = (35./16.-35./8.*X2+7./4.*X2*X2-1./6.*X2*X2*X2)*EXG
            DN = DN+DNY
            AN = AN+ANY
            EN = EN+ENY
         enddo
         ESR = EN
         AS = AN
         DNN = ASN-AS
         DF = DNN/DN
         EFST = EF
         IF(abs(DNN) .LT. 1.E-3) exit
         EF = EF+DF
      enddo

      ESM = EN*2
      ESM = ESP-ESM

   END SUBROUTINE STRUT


   function get_e_shell(mz,mn,bt2,bt3,bt4) result(esm)
      implicit none
      integer, intent(in) :: mz, mn
      real, intent(in) :: bt2, bt3, bt4
      integer :: ma
      real(r8) :: esm_n, esm_p
      real :: esm

      ma = mz+mn

      iz = mz
      in = mn
      beta1 = 0.
      beta2 = bt2
      beta3 = bt3
      beta4 = bt4

      ! calc single particle energy levles
      call WSBETA

      ! calc Eshell of neutrons
      call STRUT(ma, mn, spl_n, esm_n)

      ! calc Eshell of protons
      call STRUT(ma, mz, spl_p, esm_p)

      ! calc final shell correction, with parameters WS3.2
      ESM = real(0.7274 * (esm_n+esm_p))

   end function get_e_shell



end module mod_shell_correction



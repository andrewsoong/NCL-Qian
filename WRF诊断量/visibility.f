C NCLFORTSTART
      SUBROUTINE CALVIS(IM,JM,QV,QC,QR,QI,QS,TT,PP,VIS)
C
C   This routine computes horizontal visibility at the
C   surface or lowest model layer, from qc, qr, qi, and qs.  
C   qv--water vapor mixing ratio (kg/kg)
C   qc--cloud water mixing ratio (kg/kg)
C   qr--rain water mixing ratio  (kg/kg)
C   qi--cloud ice mixing ratio   (kg/kg)
C   qs--snow mixing ratio        (kg/kg)
C   tt--temperature              (k)
C   pp--pressure                 (Pa)
C
C   If iice=0:
C      qprc=qr     qrain=qr and qclw=qc if T>0C
C      qcld=qc          =0          =0  if T<0C
C                  qsnow=qs and qclice=qc  if T<0C
C                       =0            =0   if T>0C
C   If iice=1:
C      qprc=qr+qs   qrain=qr and qclw=qc
C      qcld=qc+qi   qsnow=qs and qclice=qc
C
C   Independent of the above definitions, the scheme can use different
C   assumptions of the state of hydrometeors:
C        meth='d': qprc is all frozen if T<0, liquid if T>0
C        meth='b': Bocchieri scheme used to determine whether qprc
C           is rain or snow. A temperature assumption is used to
C           determine whether qcld is liquid or frozen.
C        meth='r': Uses the four mixing ratios qrain, qsnow, qclw,
C           and qclice
C
C   The routine uses the following
C   expressions for extinction coefficient, beta (in km**-1),
C   with C being the mass concentration (in g/m**3):
C
C      cloud water:  beta = 144.7 * C ** (0.8800)
C      rain water:   beta =  2.24 * C ** (0.7500)
C      cloud ice:    beta = 327.8 * C ** (1.0000)
C      snow:         beta = 10.36 * C ** (0.7776)
C
C   These expressions were obtained from the following sources:
C
C      for cloud water: from Kunkel (1984)
C      for rainwater: from M-P dist'n, with No=8e6 m**-4 and
C         rho_w=1000 kg/m**3
C      for cloud ice: assume randomly oriented plates which follow
C         mass-diameter relationship from Rutledge and Hobbs (1983)
C      for snow: from Stallabrass (1985), assuming beta = -ln(.02)/vis
C
C   The extinction coefficient for each water species present is
C   calculated, and then all applicable betas are summed to yield
C   a single beta. Then the following relationship is used to
C   determine visibility (in km), where epsilon is the threshhold
C   of contrast, usually taken to be .02:
C
C      vis = -ln(epsilon)/beta      [found in Kunkel (1984)]
C
C------------------------------------------------------------------
C      implicit none

      INTEGER IM,JM
      REAL QV(IM,JM),QC(IM,JM),QI(IM,JM),QR(IM,JM),QS(IM,JM)
     1,    TT(IM,JM),PP(IM,JM),VIS(IM,JM)
C NCLEND
      REAL CELKEL,TICE,COEFLC,COEFLP,COEFFC,COEFFP,EXPONLC,EXPONLP,
     &     EXPONFC,EXPONFP,CONST1,RHOICE,RHOWAT
      REAL QPRC,QCLD,QRAIN,QSNOW,QCLW,QCLICE,TV,H1,D608,RHOAIR,
     &     RD,VOVERMD,CONCLC,CONCLP,CONCFC,CONCFP,BETAV
      INTEGER I,J
C
      CHARACTER METH*1
C------------------------------------------------------------------
C------------------------------------------------------------------
      D608=0.608
      H1=1.0
      RD=287.04
      CELKEL=273.15
      TICE=CELKEL-10.
      COEFLC=144.7
      COEFLP=2.24
      COEFFC=327.8
      COEFFP=10.36
      EXPONLC=0.8800
      EXPONLP=0.7500
      EXPONFC=1.0000
      EXPONFP=0.7776
      CONST1=-LOG(.02)
      RHOICE=917.
      RHOWAT=1000.
C
      DO J=1,JM
      DO I=1,IM
C       IF(IICE.EQ.0)THEN
C         QPRC=QR
C         QCLD=QC
C         IF(TT.LT.CELKEL)THEN
C           QRAIN=0.
C           QSNOW=QPRC
C           QCLW=0.
C           QCLICE=QCLD
C         ELSE
C           QRAIN=QPRC
C           QSNOW=0.
C           QCLW=QCLD
C           QCLICE=0.
C         ENDIF
C       ELSE
          QPRC=QR(I,J)+QS(I,J)
          QCLD=QC(I,J)+QI(I,J)
          QRAIN=QR(I,J)
          QSNOW=QS(I,J)
          QCLW=QC(I,J)
          QCLICE=QI(I,J)
C       ENDIF
C       TV=VIRTUAL(TT,QV)
        TV=TT(I,J)*(H1+D608*QV(I,J))
        RHOAIR=PP(I,J)/(RD*TV)
C       IF(METH.EQ.'D')THEN
C         IF(TT.LT.CELKEL)THEN
C           VOVERMD=(1.+QV)/RHOAIR+(QPRC+QCLD)/RHOICE
C           CONCLC = 0.
C           CONCLP = 0.
C           CONCFC = QCLD/VOVERMD*1000.
C           CONCFP = QPRC/VOVERMD*1000.
C         ELSE
C           VOVERMD=(1.+QV)/RHOAIR+(QPRC+QCLD)/RHOWAT
C           CONCLC = QCLD/VOVERMD*1000.
C           CONCLP = QPRC/VOVERMD*1000.
C           CONCFC = 0.
C           CONCFP = 0.
C         ENDIF
C       ELSEIF(METH.EQ.'B')THEN
C         IF(TT.LT.TICE)THEN
C           VOVERMD=(1.+QV)/RHOAIR+(QPRC+QCLD)/RHOICE
C           CONCLC = 0.
C           CONCLP = 0.
C           CONCFC = QCLD/VOVERMD*1000.
C           CONCFP = QPRC/VOVERMD*1000.
C         ELSEIF(PRSNOW.GE.50.)THEN
C           VOVERMD=(1.+QV)/RHOAIR+QPRC/RHOICE+QCLD/RHOWAT
C           CONCLC = QCLD/VOVERMD*1000.
C           CONCLP = 0.
C           CONCFC = 0.
C           CONCFP = QPRC/VOVERMD*1000.
C         ELSE
C           VOVERMD=(1.+QV)/RHOAIR+(QPRC+QCLD)/RHOWAT
C           CONCLC = QCLD/VOVERMD*1000.
C           CONCLP = QPRC/VOVERMD*1000.
C           CONCFC = 0.
C           CONCFP = 0.
C         ENDIF
C       ELSEIF(METH.EQ.'R')THEN
          VOVERMD=(1.+QV(I,J))/RHOAIR+(QCLW+QRAIN)/RHOWAT+
     1            (QCLICE+QSNOW)/RHOICE
          CONCLC=AMAX1(0., QCLW/VOVERMD*1000.)
          CONCLP=AMAX1(0., QRAIN/VOVERMD*1000.)
          CONCFC=AMAX1(0., QCLICE/VOVERMD*1000.)
          CONCFP=AMAX1(0., QSNOW/VOVERMD*1000.)
C       ENDIF
        BETAV=COEFFC*CONCFC**EXPONFC+COEFFP*CONCFP**EXPONFP
     1       +COEFLC*CONCLC**EXPONLC+COEFLP*CONCLP**EXPONLP
     2       +1.E-10
C CHANGED GSM 3-10-00 -->  no point in distinguishing values
C       above 20 km, so make that value the max (prev max was 80)
C        VIS(I,J)=1.E3*MIN(20.,CONST1/BETAV)   ! max of 20km
C Chuang: Per Geoff, the max visibility was changed to be 
C         cosistent with visibility ceiling in obs
C change max to be consistent with obs
      VIS(I,J)=1.E3*MIN(24.135,CONST1/BETAV)   
      ENDDO
      ENDDO
C
      RETURN
      END


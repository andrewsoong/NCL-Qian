C NCLFORTSTART
      SUBROUTINE CAL_BRN2D(NLONS,NLATS,NLVLS,
     &             P,H,T,TD,U,V,CAPE,BRN)
C    INPUT   P ======> PRESSURE   (HPA) 
C            H ======> GEOPOTENTIAL HEIGHT  (GPM)
C            T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT      (CELSIUS DEGREE)
C            U ======> EAST-WEST WIND (M/S)
C            V ======> NORTH-SOUTH WIND (M/S)
C            LEV ====> LEVEL
C            CAPE ===> EQUIVALENT POTENTIAL TEMPERATURE (J/KG)
      
      INTEGER  NLONS,NLATS,NLVLS
      REAL P(NLONS,NLATS,NLVLS), H(NLONS,NLATS,NLVLS), 
     &     T(NLONS,NLATS,NLVLS), TD(NLONS,NLATS,NLVLS), 
     &     U(NLONS,NLATS,NLVLS), V(NLONS,NLATS,NLVLS),
     &     CAPE(NLONS,NLATS), BRN(NLONS,NLATS)
C NCLEND
      real WD(NLVLS), WS(NLVLS)
      do j=1,NLATS
        do i=1,NLONS
          call FORM_TRAN_UV_WDS(U(i,j,:),V(i,j,:),NLVLS,WD,WS)
          call CAL_BRN(P(i,j,:),H(i,j,:),T(i,j,:),TD(i,j,:),
     &                 WD(:),WS(:),NLVLS,CAPE(i,j),BRN(i,j))
        end do
      end do

      END SUBROUTINE CAL_BRN2D

C NCLFORTSTART
      SUBROUTINE CAL_SI2D(NLONS,NLATS,NLVLS,
     &             P,T,TD,SI)
C    INPUT   P ======> PRESSURE   (HPA) 
C            T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT      (CELSIUS DEGREE)
      INTEGER  NLONS,NLATS,NLVLS
      REAL P(NLONS,NLATS,NLVLS), T(NLONS,NLATS,NLVLS), 
     &     TD(NLONS,NLATS,NLVLS), SI(NLONS,NLATS)
C NCLEND
      real WD(NLVLS), WS(NLVLS)

      do j=1,NLATS
        do i=1,NLONS
          call CAL_SI(P(i,j,:),T(i,j,:),TD(i,j,:),
     &                 NLVLS,SI(i,j))
        end do
      end do

      END SUBROUTINE CAL_SI2D

C NCLFORTSTART
      SUBROUTINE CAL_EHI2D(NLONS,NLATS,NLVLS,
     &             P,H,U,V,CAPE,EHI)
C     INPUT  P ======> PRESSURE   (HPA) 
C            H ========> GEOPOTENTIAL HEIGHT (GPM)
C            U ======> EAST-WEST WIND (M/S)
C            V ======> NORTH-SOUTH WIND (M/S)
C            CAPE =====> CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
C     OUTPUT EHI ======> ENERGETIC HELICITY INDEX
      INTEGER  NLONS,NLATS,NLVLS
      REAL P(NLONS,NLATS,NLVLS), H(NLONS,NLATS,NLVLS), 
     &     U(NLONS,NLATS,NLVLS), V(NLONS,NLATS,NLVLS),
     &     CAPE(NLONS,NLATS), EHI(NLONS,NLATS)
C NCLEND
      real WD(NLVLS), WS(NLVLS)

      do j=1,NLATS
        do i=1,NLONS
          call FORM_TRAN_UV_WDS(U(i,j,:),V(i,j,:),NLVLS,WD,WS)
          call CAL_EHI(P(i,j,:),H(i,j,:),WD(:),WS(:),
     &                 NLVLS,CAPE(i,j),EHI(i,j))
        end do
      end do

      END SUBROUTINE CAL_EHI2D


      SUBROUTINE CAL_ES(T,LEV,ES)
C     THIS SUBROUTINE CALCULATE SATURATE VAPOR PRESSURE
C     NOTICE THE FORMULA ARE DIFFERENT WHEN TEMPERATURE 
C      IS ABOVE OR BELOW ZERO DEGREE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            LEV ====> LEVEL
C    OUTPUT  ES  ====> SATURATE VAPOR PRESSURE (HPA)
      INTEGER LEV
      DIMENSION T(LEV),ES(LEV)
      DO K=1,LEV
          IF(T(K).GT.-150.AND.T(K).LT.60.)THEN
              IF(T(K).GE.0)ES(K)=6.112*10**(7.5*T(K)/(237.3+T(K)))
              IF(T(K).LT.0)ES(K)=6.112*10**(9.5*T(K)/(265.5+T(K)))
          ELSE
              ES(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_E(TD,LEV,E)
C     THIS SUBROUTINE CALCULATE VAPOR PRESURE
C     NOTICE THE FORMULA ARE DIFFERENT WHEN DEW POINT TEMPERATURE 
C      IS ABOVE OR BELOW ZERO DEGREE
C    INPUT   TD =====> DEW POINT      (CELSIUS DEGREE)
C            LEV ====> LEVEL
C    OUTPUT  E  ====> VAPOR PRESSURE (HPA)
      INTEGER LEV
      DIMENSION TD(LEV),E(LEV)
      DO K=1,LEV
          IF(TD(K).GT.-150.AND.TD(K).LT.60.)THEN
              IF(TD(K).GE.0)E(K)=6.112*10**(7.5*TD(K)/(237.3+TD(K)))
              IF(TD(K).LT.0)E(K)=6.112*10**(9.5*TD(K)/(265.5+TD(K)))
          ELSE
              E(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_Q(TD,P,LEV,Q)
C     THIS SUBROUTINE CALCULATE SPECIFIC HUMIDITY AND SATURATE SPECFIC HUMIDITY
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_E
C    INPUT   TD ======> DEW POINT TEMPERATURE (CELSIUS DEGREE)
C            P  =====> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  Q   ====> SPECIFIC HUMIDITY
C    TEMP ARRAY
C            E  =====> VAPOR PRESSURE   (HPA)
      INTEGER LEV
      DIMENSION TD(LEV),P(LEV),Q(LEV)
      DIMENSION E(LEV)
      CALL CAL_E(TD,LEV,E)
      DO K=1,LEV
            IF(E(K).LT.100.)THEN
                Q(K)=0.622*E(K)/(P(K)-0.378*E(K))
          ELSE
              Q(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_QS(T,P,LEV,QS)
C     THIS SUBROUTINE CALCULATE SPECIFIC HUMIDITY AND SATURATE SPECFIC HUMIDITY
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_ES
C    INPUT   TD =====> DEW POINT      (CELSIUS DEGREE)
C            P  =====> PRESSURE   (HPA)
C            LEV ====> LEVEL
C            QS  ====> SATURATE SPECIFIC HUMIDITY
C    TEMP ARRAY
C            ES =====> SATURATED VAPOR PRESSURE    (HPA)
      INTEGER LEV
      DIMENSION T(LEV),P(LEV),QS(LEV)
      DIMENSION ES(LEV)
      CALL CAL_ES(T,LEV,ES)
      DO K=1,LEV
            IF(ES(K).LT.100.)THEN
                QS(K)=0.622*ES(K)/(P(K)-0.378*ES(K))
          ELSE
              QS(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END
      
      SUBROUTINE CAL_R(TD,P,LEV,R)              
C     THIS SUBROUTINE CALCULATE MIXING RATION
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_E
C    INPUT   TD =====> DEW POINT      (CELSIUS DEGREE)
C            P  =====> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  R   ====> MIXING RATION
C    TEMP ARRAY
C            E  =====> VAPOR PRESSURE   (HPA)
      INTEGER LEV
      DIMENSION TD(LEV),P(LEV),R(LEV)
      DIMENSION E(LEV)
      CALL CAL_E(TD,LEV,E)
      DO K=1,LEV
            IF(E(K).LT.100.)THEN
                R(K)=0.622*E(K)/(P(K)-E(K))
          ELSE
              R(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_RH(T,TD,LEV,RH)
C     THIS SUBROUTINE CALCULATE RELATIVE HUMIDITY
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_E AND CAL_ES
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT      (CELSIUS DEGREE)
C            LEV ====> LEVEL
C    OUTPUT  RH  ====> RELATIVE HUMIDITY
C    TEMP ARRAY
C            E  =====> VAPOR PRESSURE     (HPA)
C            ES =====> SATURATED VAPOR PRESSURE   (HPA)
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),RH(LEV)
      DIMENSION E(LEV),ES(LEV)
      CALL CAL_E(TD,LEV,E)
      CALL CAL_ES(T,LEV,ES)
      DO K=1,LEV
          IF(E(K).LT.100.AND.ES(K).LT.100)THEN
              RH(K)=E(K)/ES(K)
          ELSE
              RH(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_D(T,TD,LEV,D)
C     THIS SUBROUTINE CALCULATE SATURATED DIFFICIENCY
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_E AND CAL_ES
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT      (CELSIUS DEGREE)
C            LEV ====> LEVEL
C    OUTPUT  D   ====> SATURATED VAPOR DIFFICIENCY  (HPA)
C    TEMP ARRAY
C            E  =====> VAPOR PRESSURE    (HPA)
C            ES =====> SATURATED VAPOR PRESSURE  (HPA)
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),D(LEV)
      DIMENSION E(LEV),ES(LEV)
      CALL CAL_E(TD,LEV,E)
      CALL CAL_ES(T,LEV,ES)
      DO K=1,LEV
          IF(E(K).LT.100.AND.ES(K).LT.100)THEN
              D(K)=ES(K)-E(K)
          ELSE
              D(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_TTD(T,TD,LEV,TTD)
C     THIS SUBROUTINE CALCULATE DEW POINT DIFFICIENCY
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT      (CELSIUS DEGREE)
C            LEV ====> LEVEL
C    OUTPUT  TTD   ====> DEW DIFFICIENCY (CELSIUS DEGREE)
C    TEMP ARRAY
C            E  =====> VAPOR PRESSURE  (HPA)
C            ES =====> SATURATED VAPOR PRESSURE  (HPA)
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),TTD(LEV)
      DO K=1,LEV
          IF(T(K).LT.70.AND.TD(K).LT.70)THEN
              TTD(K)=T(K)-TD(K)
          ELSE
              TTD(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END
   
      SUBROUTINE CAL_TV(T,TD,P,LEV,TV)
C     THIS SUBROUTINE CALCULATE VIRTUAL TEMPERATURE
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_Q
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT      (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  TV  ====> VIRTUAL TEMPERATURE   (KELVIN)
C    TEMP ARRAY
C            Q  =====> SPECFIC HUMIDITY
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),TV(LEV),P(LEV)
      DIMENSION Q(LEV)
      CALL CAL_Q(TD,P,LEV,Q)
      DO K=1,LEV
          IF(T(K).LT.70.AND.Q(K).LT.70)THEN
              TV(K)=(1+0.61*Q(K))*(T(K)+273.15)
          ELSE
              TV(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_TH(T,P,LEV,TH)
C     THIS SUBROUTINE CALCULATE POTENTIAL TEMPERATURE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  TH  ====> POTENTIAL TEMPERATURE   (KELVIN)
      INTEGER LEV
      DIMENSION T(LEV),P(LEV),TH(LEV)
      DO K=1,LEV
          IF(T(K).LT.70)THEN
              TH(K)=(T(K)+273.15)*(1000./P(K))**(287.04/1005.)
          ELSE
              TH(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_THV(T,TD,P,LEV,THV)
C     THIS SUBROUTINE CALCULATE VIRTUAL POTENTIAL TEMPERATURE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  THV ====> POTENTIAL TEMPERATURE   (KELVIN)
C    TEMP ARRAY
C            Q ======> SPECIFIC HUMIDITY
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),P(LEV),THV(LEV)
      DIMENSION Q(LEV)
      CALL CAL_Q(TD,P,LEV,Q)
      DO K=1,LEV
          IF(T(K).LT.70.AND.Q(K).LT.1)THEN
              THV(K)=(T(K)+273.15)*(1+0.61*Q(K))
     &               *(1000./P(K))**(287.04/1005.)
          ELSE
              THV(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_THE(T,TD,P,LEV,THE)
C     THIS SUBROUTINE CALCULATE EQUIVALENT POTENTIAL TEMPERATURE
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_E, CAL_R AND CAL_RH
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  THE ====> EQUIVALENT POTENTIAL TEMPERATURE   (KELVIN)
C    TEMP ARRAY
C            R ======> MIXING RATIO
C            RH =====> RELATIVE HUMIDITY
C            E ======> VAPOR PRESSURE
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),P(LEV),THE(LEV)
      DIMENSION R(LEV),E(LEV),RH(LEV)
      REAL RD,RV,CPD
      RD=287.04
      RV=461.50
      CPD=1005.
      CALL CAL_E(TD,LEV,E)
      CALL CAL_R(TD,P,LEV,R)
      CALL CAL_RH(T,TD,LEV,RH)
      DO K=1,LEV
C          P0=P(K)
C            T0=T(K)
C            TD0=TD(K)
C          CALL CAL_TC_PC(P0,T0,TD0,PC0,TC0,THSE0)
C          THE(K)=THSE0
          IF(T(K).LT.70.AND.R(K).LT.1)THEN
              THE(K)=(T(K)+273.15)*(1000./(P(K)-E(K)))**(RD/CPD)
     1               *RH(K)**(R(K)*RV/CPD)
     2               *EXP((2500000.-2368*T(K))*R(K)/(CPD*(T(K)+273.15)))
          ELSE
              THE(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_IL(T,TD,P,LEV,IL)
C     THIS SUBROUTINE CALCULATE INDEX OF CONDITIONAL STABILITY
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_THE, CAL_TE AND CAL_THES
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  IL =====> INDEX OF CONDITIONAL STABILITY
C    TEMP ARRAY
C            THE ====> EQUIVALENT OF POTENTIAL TEMPERAURE
C            THES ===> SATURATED EQUIVALENT OF POTENTIAL TEMPERATURE
C            THESD ==> ENVIRONMENT OF SATURATE EQUIVALENT POTENTIAL TEMPERATURE
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),P(LEV),THE(LEV)
     1,THES(LEV),THESD(110)
      REAL IL,THE0,THES500
      IF(P(1).LT.500.)THEN
          IL=99999.9
          RETURN
      ELSE
            CALL CAL_THE(T,TD,P,LEV,THE)
          THE0=THE(1)
          CALL CAL_THES(T,P,LEV,THES)
          CALL CAL_TE(P,THES,THESD,LEV)
          P_L=INT(P(1)/10.)*10.
          DO K=1,70
              IF(P_L.EQ.500)THEN
                  THES500=THESD(K)
                  IL=THES500-THE0
                  RETURN
              ENDIF
              P_L=P_L-10.
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_ILC(T,TD,P,LEV,ILC)
C     THIS SUBROUTINE CALCULATE INDEX OF CONDITIONAL COVECTIVE STABILITY
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_THE, CAL_THES AND CAL_TE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  ILC ====> INDEX OF CONDITIONAL CONVECTIVE STABILITY
C    TEMP ARRAY
C            THE ====> EQUIVALENT OF POTENTIAL TEMPERAURE
C            THED ===> ENVIRONMENT EQUIVALENT POTENTIAL TEMPERAURE
C            THES ===> SATURATED EQUIVALENT OF POTENTIAL TEMPERATURE
C            THESD ==> ENVIRONMENT OF SATURATE EQUIVALENT POTENTIAL TEMPERATURE
      INTEGER LEV
      DIMENSION T(LEV),TD(LEV),P(LEV),THE(LEV)
     1,THES(LEV),THED(110),THESD(110)
      REAL ILC,THE850,THE500,THES500,THE0
      IF(P(1).LT.850.)THEN
          ILC=99999.9
          RETURN
      ELSE
            CALL CAL_THE(T,TD,P,LEV,THE)
          THE0=THE(1)
           CALL CAL_TE(P,THE,THED,LEV)
          CALL CAL_THES(T,P,LEV,THES)
          CALL CAL_TE(P,THES,THESD,LEV)
          P_L=INT(P(1)/10.)*10.
          DO K=1,70
                IF(P_L.EQ.850)THE850=THED(K)
              IF(P_L.EQ.500)THEN
                  THES500=THESD(K)
                  THE500=THED(K)
                  ILC=THES500-THE0+THE500-THE850
                  RETURN
              ENDIF
              P_L=P_L-10.
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_THES(T,P,LEV,THES)
C     THIS SUBROUTINE CALCULATE SATURATE EQUIVALENT POTENTIAL TEMPERATURE
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_ES AND CAL_R
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  THES ====> SATURATE EQUIVALENT POTENTIAL TEMPERATURE   (KELVIN)
C    TEMP ARRAY
C            RS ======> SATURATE MIXING RATIO
C            ES ======> SATURATE VAPOR PRESSURE
      INTEGER LEV
      DIMENSION T(LEV),P(LEV),THES(LEV)
      DIMENSION RS(LEV),ES(LEV)
      REAL RD,CPD
      RD=287.04
      CPD=1005.
      CALL CAL_ES(T,LEV,ES)
      CALL CAL_R(T,P,LEV,RS)
      DO K=1,LEV
          IF(T(K).LT.70.AND.RS(K).LT.1)THEN
              THES(K)=(T(K)+273.15)*(1000./(P(K)-ES(K)))**(RD/CPD)
     1              *EXP((2500000.-2368*T(K))*RS(K)/(CPD*(T(K)+273.15)))
          ELSE
              THES(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_TSIG(P,T,TD,H,LEV,TSIG)
C     THIS SUBROUTINE CALCULATE TOTAL TEMPRATURE AT EACH LEVEL
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_Q
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            H ======> GEOPOTENTIAL HEIGHT  (GPM)
C            TD =====> DEW POINT  (CELSIUS)
C            LEV ====> LEVEL
C    OUTPUT  TSIG ===> TOTAL TEMPERATURE  (CELSIUS)
C    TEMP ARRAY
C            Q ======> SPECIFIC HUMIDITY
      INTEGER LEV
      DIMENSION T(LEV),P(LEV),TD(LEV),H(LEV),TSIG(LEV)
      DIMENSION Q(LEV)
      REAL CPD
      CPD=1005.
      CALL CAL_Q(TD,P,LEV,Q)
      DO K=1,LEV
          TSIG(K)=T(K)+(2500000.-2368*T(K))*Q(K)/CPD+9.8*H(K)/CPD
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_TSIGS(P,T,H,LEV,TSIGS)
C     THIS SUBROUTINE CALCULATE TOTAL TEMPRATURE AT EACH LEVEL
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_Q
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            H ======> GEOPOTENTIAL HEIGHT  (GPM)
C            LEV ====> LEVEL
C    OUTPUT  TSIGS ==> SATURATE TOTAL TEMPERATURE  (CELSIUS)
C    TEMP ARRAY
C            QS =====> SATURATE SPECIFIC HUMIDITY
      INTEGER LEV
      DIMENSION T(LEV),P(LEV),H(LEV),TSIGS(LEV)
      DIMENSION QS(LEV)
      REAL CPD
      CPD=1005.
      CALL CAL_Q(T,P,LEV,QS)
      DO K=1,LEV
          TSIGS(K)=T(K)+(2500000.-2368*T(K))*QS(K)/CPD+9.8*H(K)/CPD
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_F(P,T,LEV,F)
C     THIS SUBROUTINE CALCULATE FUNCTION OF CONDENSATION
C     INPUT  P ========> PRESSURE     (HPA)
C            T ========> TEMPERATURE  (CELSIUS DEGREE)
C            LEV ======> NUMBER OF LAYER
C     OUTPUT F ========> FUNCTION OF CONDENSATION       (1/HPA)
C     TEMP ARRAY
C            QS =======> SATURATED SPECIFIC HUMIDITY
      DIMENSION P(LEV),T(LEV),F(LEV)
      DIMENSION QS(LEV)
      INTEGER LEV
      REAL CPD,RV,RD,LV
      CPD=1005.
      RV=461.
      RD=287.
      CALL CAL_QS(T,P,LEV,QS)
      DO K=1,LEV
          TBIG=T(K)+273.15
          LV=2500000.-2368*T(K)
          F(K)=QS(K)*TBIG/P(K)*(LV*RD-CPD*RV*TBIG)
     1         /(CPD*RV*TBIG**2+QS(K)*LV**2)
      ENDDO
      RETURN
      END

       SUBROUTINE CAL_PW(TD,P,LEV,PW)
C     THIS SUBROUTINE CALCULATE TOTAL QUALITY WATER
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_Q
C    INPUT   TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  QW =====> POTENTIAL QUALITY WATER  (CM)
C    TEMP ARRAY
C            Q ======> SPECIFIC HUMIDITY
      INTEGER LEV
      DIMENSION TD(LEV),P(LEV)
      DIMENSION Q(LEV)
      REAL PW,G
      G=9.8
      CALL CAL_Q(TD,P,LEV,Q)
      PW=0.0
      IF(Q(1).GT.1.OR.Q(2).GT.1)THEN
          PW=99999.9
          RETURN
      ENDIF
      DO K=1,LEV-1
          IF(Q(K).LT.1.0.AND.Q(K+1).LT.1.0)THEN
              PW=PW+(Q(K)+Q(K+1))*(P(K)-P(K+1))/.98
          ELSE
              RETURN
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_THSE(P,T,TD,LEV,THSE)
C     THIS SUBROUTINE CALCULATE PSEUDO-EQUIVALENT POTENTIAL TEMPERATURE
C      AND IT IS USED WITH CAL_TC_PC,CAL_TA,CAL_TE,CAL_ENGY,
C      IN WHICH SUBROUTINE CAL_TA AND CAL_TC_PC ARE USED WITH CAL_THSE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  THSE ===> PSEUDO EQUIVALENT POTENTIAL TEMPERATURE   (KELVIN)
C    TEMP ARRAY
      DIMENSION P(LEV),T(LEV),TD(LEV),THSE(LEV)
      REAL A,B,THSE,TCM,PCM,ECM
      DO K=1,LEV
            IF(T(K).LT.80.AND.TD(K).LT.80)THEN
                TCM=T(K)-(T(K)-TD(K))*0.976/(0.976-0.000833
     1            *(237.3+TD(K))**2/(273.15+TD(K)))
              PCM=P(K)*((273.15+TCM)/(273.15+T(K)))**3.5
              IF(TCM.GE.0)A=7.5*TCM/(237.3+TCM)
              IF(TCM.LT.0)A=9.5*TCM/(265.5+TCM)
              ECM=6.112*10.**A
              B=1000.0/(PCM-ECM)
      WRITE(*,*)PCM
              THSE(K)=(273.16+TCM)*B**0.286*EXP((2.5*10.**6-2368.*TCM)
     1                *0.622*ECM/(1004.*(273.16+TCM)*(PCM-ECM)))
          ELSE
              THSE(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_PLCL_TLCL_THSE(P,T,TD,LEV,PLCL,TLCL,THSE)
C     THIS SUBROUTINE CALCULATE LIFTING CONDENSATION LAYER PRESSURE AND TEMPERATURE
C      AND IT IS USED WITH CAL_E,CAL_R,
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  PLCL ===> LIFTING CONDENSATION LAYER (HPA)
C    TEMP ARRAY
C            R ======> MIXING RATIO
C            E ======> VAPOR PRESSURE    (HPA)
      DIMENSION P(LEV),T(LEV),TD(LEV),PLCL(LEV),TLCL(LEV),THSE(LEV)
      DIMENSION E(LEV),R(LEV)
      REAL CPD,CPV,RD
      RD=287.04
      CPD=1005.
      CPV=1850.
      CALL CAL_E(TD,LEV,E)
      CALL CAL_R(TD,P,LEV,R)
      DO K=1,LEV
            IF(E(K).LT.300.)THEN
              TLCL(K)=2840./(3.5*LOG(T(K)+273.15)-LOG(E(K))-4.805)+55
              PLCL(K)=P(K)*EXP((CPD+R(K)*CPV)/(RD*(1+R(K)/0.622))
     1                *LOG(TLCL(K)/(T(K)+273.15)))
              TCM=TLCL(K)-273.15
              PCM=PLCL(K)
              IF(TCM.GE.0)A=7.5*TCM/(237.3+TCM)
              IF(TCM.LT.0)A=9.5*TCM/(265.5+TCM)
              ECM=6.112*10.**A
              B=1000.0/(PCM-ECM)
              THSE(K)=(273.16+TCM)*B**0.286*EXP((2.5*10.**6-2368.*TCM)
     1                *0.622*ECM/(1004.*(273.16+TCM)*(PCM-ECM)))
          ELSE
              TLCL(K)=99999.9
              PLCL(K)=99999.9
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_BCAPE(P,T,TD,LEV,BCAPE)
C     THIS SUBROUTINE CALCULATE THE BEST COVECTIVE AVAILABLE POTENTIAL ENERGY
C      THIS PROGRAM IS USED WITH SUB(CAL_TC_PC,CAL_TA,CAL_TE,CAL_ENGY,
C      IN WHICH SUBROUTINE CAL_TA AND CAL_TC_PC ARE USED WITH CAL_THSE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT
C            BCAPE ==> BEST CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
C     TEMP ARRAY
C    INPUT   TP =====> TEMPERATURE FROM ABOVE MOIST UNSTABLE LAYER
C                      (CELSIUS DEGREE)
C            TDP ====> DEW POINT FROM ABOVE MOIST UNSTABLE LAYER
C                      (CELSIUS DEGREE)
C            PP =====> PRESSURE FROM ABOVE MOIST UNSTABLE LAYER  (HPA)
      DIMENSION P(LEV),T(LEV),TD(LEV),PP(LEV),TP(LEV),TDP(LEV)
      DIMENSION TE(120),TA(120)
      REAL BCAPE,TC,PC,THSE,PC0,TC0,TD0,THSE0
      INTEGER KP
      KP=0
      THSE=220.
      DO K=1,LEV
          P0=P(K)
            T0=T(K)
            TD0=TD(K)
          CALL CAL_TC_PC(P0,T0,TD0,PC0,TC0,THSE0)
          IF(P(1)-P0.LT.300.AND.THSE0.GT.THSE)THEN
              THSE=THSE0
              PC=PC0
              TC=TC0
              KC=K
          ENDIF
      ENDDO
      KL=LEV-KC+1
      DO K=1,KL
          PP(K)=P(K+KC-1)
          TP(K)=T(K+KC-1)
          TDP(K)=TD(K+KC-1)
      ENDDO
      P0=PP(1)
      T0=TP(1)
      CALL CAL_TE(PP,TP,TE,KL)
      CALL CAL_TA(P0,T0,PC,TC,THSE,TA)
      CALL CAL_VIR_ENGY(PP,TP,TDP,KL,TE,TA,BCAPE)
C      CALL CAL_ENGY(P0,TE,TA,BCAPE)
      RETURN
      END

      SUBROUTINE CAL_BLI(P,T,TD,LEV,BLI)
C     THIS SUBROUTINE CALCULATE THE BEST LIFING INDEX
C      NOTICE THIS PROGRAM IS USED WITH CAL_TC_PC,CAL_TA,CAL_TE,COM_BLI,
C      IN WHICH SUBROUTINE CAL_TA AND CAL_TC_PC ARE USED WITH CAL_THSE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT
C            BLI ====> BEST LIFTING INDEX   (CELSIUS DEGREE)
C    TEMP ARRAY
C            TA =====> PARCEL TEMPERATURE  (CELSIUS DEGREE)
C            TE =====> ENVIRONMENTAL TEMPERATURE (CELSIUS DEGREE)
      DIMENSION P(LEV),T(LEV),TD(LEV),PP(LEV),TP(LEV),TDP(LEV)
      DIMENSION TE(120),TA(120)
      REAL BLI,PC0,TC0,TD0,THSE0,BLI0
      BLI=99999.9
      DO KC=1,LEV
            IF(P(KC).GE.700)THEN
              KL=LEV-KC+1
                  DO K=1,KL
                        PP(K)=P(K+KC-1)
                        TP(K)=T(K+KC-1)
                        TDP(K)=TD(K+KC-1)
                  ENDDO
              P0=PP(1)
                T0=TP(1)
                TD0=TDP(1)
              CALL CAL_TE(PP,TP,TE,KL)
                  DO K=1,KL
                        PP(K)=0.
                  ENDDO
              CALL CAL_TC_PC(P0,T0,TD0,PC0,TC0,THSE0)
              CALL CAL_TA(P0,T0,PC0,TC0,THSE0,TA)
              CALL COM_BLI(P0,TA,TE,BLI0)
              IF(BLI0.LE.BLI)BLI=BLI0
          ELSE
              RETURN
            ENDIF
      ENDDO
      RETURN
      END
      
      SUBROUTINE COM_BLI(P0,TA,TE,BLI0)
C     THIS SUBROUTINE CALCULATE THE LIFTING INDEX BY USING PARCEL 
C     AND ENVIRONMENTAL TEMPERATURE
C     INPUT  T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C     OUTPUT
C            BLI0 ===> LIFTING INDEX (CELSIUS DEGREE)
C     TEMP ARRAY
C     INPUT   TA ====> PARCEL TEMPERATURE    (CELSIUS DEGREE)
C            TE =====> ENVIRONMENTAL TEMPERATURE (CELSIUS DEGREE) 
      DIMENSION TE(120),TA(120)
      REAL P0,BLI0,P_L
      BLI0=99999.9
      P_L=INT(P0/10.)*10.
      DO K=1,70
          IF(P_L.EQ.500.)THEN
              TE500=TE(K)
              TA500=TA(K)
              BLI0=TE500-TA500
              P_L=P_L-10.
              RETURN
          ELSEIF(P_L.GT.500)THEN
              P_L=P_L-10.
          ENDIF
      ENDDO
      RETURN
      END

       SUBROUTINE CAL_TT(P,T,TD,LEV,TT)
C     THIS SUBROUTINE CALCULATE TOTAL INDEX
C      NOTICE THIS PROGRAM IS USED WITH CAL_TE AND CAL_TDE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT
C            TT ====> TOTOAL INDEX   (CELSIUS DEGREE)
C    TEMP ARRAY
C            TE =====> ENVIRONMENTAL TEMPERATURE (CELSIUS DEGREE)
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION TE(120),TDE(120)
      REAL TT,T850,TD850,T500
      TT=99999.9
      IF(P(1).LT.850)RETURN
      CALL CAL_TE(P,T,TE,LEV)
      CALL CAL_TDE(P,TD,TDE,LEV)
C      WRITE(*,*)(TE(K),K=1,55)
C      WRITE(*,*)(TDE(K),K=1,55)
C      PAUSE 999
      P_L=INT(P(1)/10.)*10.
      DO K=1,70
          IF(P_L.EQ.850)THEN
              TD850=TDE(K)
              T850=TE(K)
          ENDIF
          IF(P_L.EQ.500)THEN
                T500=TE(K)
              TT=T850-2*T500+TD850
              RETURN
          ENDIF
          P_L=P_L-10.
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_CCL(TD,T,P,LEV,CCL,TCON)
C     THIS SUBROUTINE CALCULATE CONVECTIVE CONDENSATION LEVEL
C     INPUT TD =========> DEW POINT  (CELSIUS)
C           P ==========> PRESSURE   (HPA)
C           T ==========> TEMPERATURE  (CELSIUS)
C           LEV ========> LEVEL
C     OUTPUT CCL =======> CONVECTIVE CONDENSATION LEVEL (HPA)
C            TCON ======> TEMPERATURE OF CONVECTIVE (CELSIUS)
C     TEMP ARRAY
C          TE ==========> ENVIRONMENT TEMPERATURE (CELSIUS)
C          TDS =========> SATURATED DEW POINT OF PARCEL AT EACH LEVEL
      DIMENSION TD(LEV),T(LEV),P(LEV)
      DIMENSION TE(110),TDS(110)
      INTEGER LEV
      REAL CCL
      CALL CAL_TE(P,T,TE,LEV)
      TD0=TD(1)
      CALL CAL_TDS(TD0,P,LEV,TDS)
      KLFC=0
      IF(TDS(1).GT.TE(1).OR.(TDS(1).EQ.TE(1).AND.TDS(2).GT.TE(2)))THEN
          CCL=P(1)
          TCON=T(1)
      ELSEIF((TDS(1).EQ.TE(1).AND.TDS(2).LE.TE(2)).OR
     1   .TDS(2).LE.TE(2))THEN
          P_L=INT(P(1)/10.)*10
          P_T=P_L-10
          DO K=2,90
              IF(TDS(K).GT.TE(K))THEN
                  KLFC=K
                  P_B=P_T+10
                  DELT1=TE(K-1)-TDS(K-1)
                  DELT2=TDS(K)-TE(K)
                  CCL=(P_B*DELT2+P_T*DELT1)/(DELT1+DELT2)
                  CCT=TE(K-1)-DELT1/(DELT1+DELT2)*(TE(K-1)-TE(K))
                  THCL=(CCT+273.15)*(1000./CCL)**(287/1005.)
                  TCON=THCL*(1000/P(1))**(-287./1005.)-273.15
                  GOTO 290
               ENDIF
          P_T=P_T-10
          ENDDO
         IF(KLFC.EQ.0)CCL=99999.9
      ENDIF
290      CONTINUE
      RETURN
      END

      SUBROUTINE CAL_TDS(TD0,P,LEV,TDS)
C     THIS SUBROUTINE CALCULATE LAPSE RATE OF DEW POINT
C     INPUT TDO ======> DEW POINT OF THE FIRST PRESSURE LEVEL (CELSIUS)
C           P ========> PRESSURE  (HPA)
C           LEV ======> NUMBER OF LEVEL
C     OUTPUT TDS =====> DEW POINT ALONG THE SPECIFIC Q LINE
      DIMENSION P(LEV),TDS(110)
      INTEGER LEV
      REAL C
      C=461/2500000.
      TDS(1)=TD0
C      WRITE(*,*)' TD0 = ',TD0
      P1=INT(P(1)/10.)*10
      P2=P1-10.
      DO K=2,110
         TDS(K)=1/(C*ALOG(P1/P2)+1/(TDS(K-1)+273.15))-273.15
          P1=P1-10.
          P2=P2-10.
          IF(P2.LT.300.)RETURN
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_TDS0(TD0,P,LEV,TDS)
C     THIS SUBROUTINE CALCULATE LAPSE RATE OF DEW POINT USING ANOTHER WAY
C     INPUT TDO ======> DEW POINT OF THE FIRST PRESSURE LEVEL (CELSIUS)
C           P ========> PRESSURE  (HPA)
C           LEV ======> NUMBER OF LEVEL
C     OUTPUT TDS =====> DEW POINT ALONG THE SPECIFIC Q LINE
      DIMENSION P(LEV),TDS(110)
      INTEGER LEV
      REAL C
      C=461/2500000.
      e0=6.112*10**(7.5*td0/(237.3+td0))
      q00=e0/(p(1)-0.378*e0)
      P1=INT(P(1)/10.)*10.
      DO K=1,110
          e1=q00*p1/(1+0.378*q00)
          c1=alog10(e1/6.112)
          tds(k)=c1*237.3/(7.5-c1)
          P1=P1-10.
          IF(P1.LT.300.)RETURN
      ENDDO
      RETURN
      END

       SUBROUTINE CAL_KI(P,T,TD,LEV,KI)
C     THIS SUBROUTINE CALCULATE K INDEX
C      NOTICE THIS PROGRAM IS USED WITH CAL_TE AND CAL_TDE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT
C            KI ====> BEST LIFTING INDEX   (CELSIUS DEGREE)
C    TEMP ARRAY
C            TE =====> ENVIRONMENTAL TEMPERATURE (CELSIUS DEGREE)
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION TE(120),TDE(120)
      REAL KI,T850,TD850,T700,TD700,T500
      KI=99999.9
      IF(P(1).LT.850)RETURN
      CALL CAL_TE(P,T,TE,LEV)
      CALL CAL_TDE(P,TD,TDE,LEV)
      P_L=INT(P(1)/10.)*10.
      DO K=1,70
          IF(P_L.EQ.850)THEN
              TD850=TDE(K)
              T850=TE(K)
          ENDIF
          IF(P_L.EQ.700)THEN
              TD700=TDE(K)
              T700=TE(K)
          ENDIF
          IF(P_L.EQ.500)THEN
                T500=TE(K)
              KI=T850-T500+TD850-T700+TD700
              RETURN
          ENDIF
          P_L=P_L-10.
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_MK(P,T,TD,LEV,MK)
C     THIS SUBROUTINE CALCULATE MODIFIED K INDEX
C      NOTICE THIS PROGRAM IS USED WITH CAL_TE AND CAL_TDE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT
C            KI ====> BEST LIFTING INDEX   (CELSIUS DEGREE)
C    TEMP ARRAY
C            TE =====> ENVIRONMENTAL TEMPERATURE (CELSIUS DEGREE)
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION TE(120),TDE(120)
      REAL MK,T850,TD850,T700,TD700,T500
      MK=99999.9
      IF(P(1).LT.850)THEN
            RETURN
      ELSE
          T0=T(1)
          TD0=TD(1)
          CALL CAL_TE(P,T,TE,LEV)
          CALL CAL_TDE(P,TD,TDE,LEV)
          P_L=INT(P(1)/10.)*10.
          DO K=1,70
              IF(P_L.EQ.850)THEN
                  TD850=TDE(K)
                  T850=TE(K)
              ENDIF
              IF(P_L.EQ.700)THEN
                  TD700=TDE(K)
                  T700=TE(K)
              ENDIF
              IF(P_L.EQ.500)THEN
                    T500=TE(K)
                  MK=(T0+T850)/2.-T500+(TD0+TD850)/2.-T700+TD700
                  RETURN
              ENDIF
              P_L=P_L-10.
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_DCI(P,T,TD,LEV,LI,DCI)
C     THIS SUBROUTINE CALCULATE DEEP CONVECTIN INDEX
C      NOTICE THIS PROGRAM IS USED WITH CAL_TE AND CAL_TDE
C     AND IT SHOULD BE SUBMITTED AFTER CAL_CAPE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C            LI LIFTING INDEX
C    OUTPUT
C            DCI ====> DEEP CONVECTION INDEX
C    TEMP ARRAY
C            TE =====> ENVIRONMENTAL TEMPERATURE (CELSIUS DEGREE)
C            TDE ====> DEW POINT
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION TE(120),TDE(120)
      REAL DCI,LI,T850,TD850
      DCI=99999.9
      IF(P(1).LT.850.)THEN
            RETURN
      ELSE
          CALL CAL_TE(P,T,TE,LEV)
          CALL CAL_TDE(P,TD,TDE,LEV)
          P_L=INT(P(1)/10.)*10.
          DO K=1,30
              IF(P_L.EQ.850)THEN
                  TD850=TDE(K)
                  T850=TE(K)
                  DCI=T850+TD850-LI
              RETURN
              ENDIF
              P_L=P_L-10.
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_MDCI(P,T,TD,LEV,LI,MDCI)
C     THIS SUBROUTINE CALCULATE MODIFIED DEEP CONVECTION INDEX
C      NOTICE THIS PROGRAM IS USED WITH CAL_TE AND CAL_TDE,
C     AND IT SHOULD BE SUBMITED AFTER CAL_CAPE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C            LI =====> LIFTING INDEX
C    OUTPUT
C            KI ====> BEST LIFTING INDEX   (CELSIUS DEGREE)
C    TEMP ARRAY
C            TE ====> ENVIRONMENTAL TEMPERATURE (CELSIUS DEGREE)
C            TDE ===> DEW POINT (CELSIUS DEGREE)
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION TE(120),TDE(120)
      REAL MDCI,LI,T850,TD850,T0,TD0
      DCI=99999.9      
      IF(P(1).LT.850.)THEN
            RETURN
      ELSE
          T0=T(1)
          TD0=TD(1)
          CALL CAL_TE(P,T,TE,LEV)
          CALL CAL_TDE(P,TD,TDE,LEV)
          P_L=INT(P(1)/10.)*10.
          DO K=1,30
              IF(P_L.EQ.850)THEN
                  TD850=TDE(K)
                  T850=TE(K)
                  MDCI=(T0+T850)/2.+(TD0+TD850)/2.-LI
              RETURN
              ENDIF
              P_L=P_L-10.
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_MDPI(P,THE,LEV,MDPI)
C     THIS SUBROUTINE CALCULATE MICRO-DOWNBURST POTENTIAL INDEX
C     INPUT P =======> PRESSURE
C           THE =====> EQUIVALENT POTENTIAL TEMPERATURE
C           LEV =====> LEVEL
C     OUTPUT
C           MDPI ====> MICRO-DOWNBURST POTENTIAL INDEX
      DIMENSION P(LEV),THE(LEV)
      INTEGER LEV
      REAL MDPI
      THEB=200.
      THET=500.
      IF(P(1).LT.850)THEN
          MDPI=99999.9
          RETURN
      ELSE
            DO K=1,LEV
              IF(THE(K).GT.THEB.AND.(P(1)-P(K)).LE.150.)THEN
                  THEB=THE(K)
              ENDIF
               IF(THE(K).LT.THET.AND.P(K).LE.650.AND.P(K).GE.500.)THEN
                  THET=THE(K)
              ENDIF
              IF(P(K).LE.500)THEN
C                        WRITE(*,*)'THEB = ',THEB,' THET = ',THET
                      MDPI=(THEB-THET)/20.
                  RETURN
               ENDIF
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_IC(P,T,TD,LEV,IC)
C     THIS SUBROUTINE CALCULATE STABILITY OF CONVECTION
C      NOTICE THIS PROGRAM IS USED WITH CAL_TE AND CAL_TDE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT
C            IC ====> INDEX OF CONVECTION  (KELVIN)
C    TEMP ARRAY
C            THE =====> EQUIVELENT OF POTENTIAL TEMPERATURE (KELVIN)
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION THE(LEV)
      REAL IC,THE850,THE500
      IF(P(1).LT.850)THEN
          IC=99999.9
          RETURN
      ELSE
          CALL CAL_THE(T,TD,P,LEV,THE)
          DO K=1,LEV
              IF(P(K).EQ.850)THEN
                  THE850=THE(K)
              ENDIF
              IF(P(K).EQ.500)THEN
                  THE500=THE(K)
                  IC=THE500-THE850
                  RETURN
              ENDIF
            ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_BIC(P,T,TD,LEV,BIC)
C     THIS SUBROUTINE CALCULATE BEST STABILITY INDEX OF CONVECTION
C      NOTICE THIS PROGRAM IS USED WITH CAL_TE AND CAL_TDE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT
C            BIC ====> BEST INDEX OF CONVECTIVE INDEX  (KELVIN)
C    TEMP ARRAY
C            THE =====> EQUIVLENT POTENTIAL TEMPERATURE (KELVIN)
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION THE(LEV)
      REAL BIC,THMAX,THMIN
      THMAX=220
      THMIN=500
      IF(P(1).LT.850.)THEN
          BIC=99999.9
          RETURN
      ELSE
          CALL CAL_THE(T,TD,P,LEV,THE)
          DO K=1,LEV
              IF(P(1)-P(K).LT.150.AND.THE(K).GT.THMAX)THEN
                  THMAX=THE(K)
              ENDIF
              IF(P(K).GE.500.AND.P(K).LE.650.AND.THE(K).LT.THMIN)THEN
                  THMIN=THE(K)
              ENDIF
              IF(P(K).LE.500)THEN
                  BIC=THMIN-THMAX
                  RETURN
               ENDIF
            ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_TW(T,TD,P,LEV,TW)
C     THIS SUBROUTINE CALCULATE WET-BULB TEMPERATURE
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY CAL_DHZHF
C     INPUT   T ========> TEMPERATURE   (CELSIUS DEGREE)
C             TD =======> DEW POINT     (CELSIUS DEGREE)
C             P ========> PRESSURE      (HPA)
C             LEV ======> LEVEL
C     OUTPUT      TW =======> WET-BULB TEMPERATURE (CELSIUS DEGREE)
C     TEMP ARRAY
C             Q ========> SPECIFIC HUMIDITY
      DIMENSION T(LEV),TD(LEV),P(LEV)
      DIMENSION TW(LEV),Q(LEV)
      INTEGER LEV
      REAL CPD,TWW
      CPD=1005.
      DO K=1,LEV
          Q(K)=0.622*6.112*10**(7.5*TD(K)/(237.3+TD(K)))/P(K)
      ENDDO
      DO K=1,LEV
          TOTAL0=T(K)+(2500000.-2368*TD(K))*Q(K)/CPD
          PK=P(K)
          T0=T(K)
          TD0=TD(K)
          CALL CAL_DHZHF(TOTAL0,T0,TD0,PK,TWW)
          TW(K)=TWW
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_DHZHF(TOTAL0,T0,TD0,PK,TW)
C     THIS SUBROUTINE CALCULATE A PROPER TEMPERATUE BY USING
C     ITERATIVE METHOD
C     NOTICE IT IS SUPPORTED BY CAL_QL
C     INPUT   TOTAL0 ==> TARGET FOR ITERATION
C             T0 ======> TEMPERATURE    (CELSIUS DEGREE)
C             TD0 =====> DEW POINT      (CELSIUS DEGREE)
C             PK ======> PRESSURE       (HPA)
C     OUTPUT  TW ======> WET-BULB TEMPERATURE (CELSIUS DEGREE)
      REAL TOTAL0,T0,TD0,PK,TW,LQ
      CPD=1005.
      TM0=T0
      TM1=TD0
101      TM=(TM0+TM1)/2.
      LQ=(2500000.-2368*TM)*6.112*10**(7.5*TM/(237.3+TM))
      TOTALM=TM+0.622*LQ/PK/CPD
      IF(ABS(TOTALM-TOTAL0).GT.0.001)THEN
          IF(TOTALM.GT.TOTAL0)THEN
              TM0=TM
              TM1=TM1
          ELSE
              TM0=TM0
              TM1=TM
          ENDIF
          GOTO 101
      ELSE
          TW=TM
      ENDIF
      RETURN
      END
     
      SUBROUTINE CAL_HT(P,T,H,LEV,TX,HT)
C     THIS SUBROUTINE CALCULATE HEIGHT FOR THE SPECIAL TEMPERATURE
C     INPUT P ========> PRESSURE  (HPA)
C           T ========> TEMPERATURE  (CELSIUS DEGREE)
C           H ========> GEOPOTENTIAL HEIGHT
C           LEV ======> LEVEL
C           TX =======> SPECIAL TEMPERATURE (CELSIUS DEGREE)
C     OUTPUT
C           HT =======> GEOPOTENTIAL HEIGHT OF THE SPECIAL TEMPERATUE
C                        LAYER
      DIMENSION P(LEV),T(LEV),H(LEV)
      INTEGER LEV
      REAL TX,PF,HT
      HT=99999.9
      IF(T(1).LE.TX)THEN
          IF(T(1).EQ.TX)HT=H(1)
          IF(T(1).LT.TX)HT=99999.9
          RETURN
      ELSE
c     lyd lev==>lev-1
          DO K=1,LEV-1
              IF(T(K+1).LT.TX)THEN
                  PF=EXP((T(K)*LOG(P(K+1))-T(K+1)*LOG(P(K))
     1                     -TX*LOG(P(K+1)/P(K)))/(T(K)-T(K+1)))
                  HT=H(K)+287.03*(2*273.15+T(K)+TX)/2./9.8*ALOG(P(K)/PF)
                  RETURN
                ENDIF
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CHANGE_DATA(FILE1,FILE2)
      CHARACTER FILE1*12,FILE2*12
      DIMENSION IA(6),RA(6)
      INTEGER NUM,IA,NSTN
      REAL RLAT,RLNG,HTST
      OPEN(7,FILE=FILE1)
      OPEN(8,FILE=FILE2)
      DO KSTN=1,999
          READ(7,'(I5,3F10.2,I5)',END=11,ERR=11)NSTN,RLAT,RLNG,HTST,NUM
          NUM=NUM/6
          WRITE(8,'(2I8)')NSTN,NUM
          DO K=1,NUM
              READ(7,'(6I7)')(IA(I),I=1,6)
              RA(1)=IA(1)
              IF(IA(2).EQ.9999)RA(2)=99999.0
              IF(IA(2).NE.9999)RA(2)=IA(2)
              RA(3)=IA(3)/10.
              IF(IA(4).NE.9999)RA(4)=IA(4)/10.
              IF(IA(4).EQ.9999)RA(4)=99999.0
              IF(IA(5).NE.9999)RA(5)=IA(5)
              IF(IA(5).EQ.9999)RA(5)=99999.0
              IF(IA(6).NE.9999)RA(6)=IA(6)
              IF(IA(6).EQ.9999)RA(6)=99999.0
              IF(K.EQ.1.AND.IA(2).EQ.9999)RA(2)=HTST
              WRITE(8,'(6F10.2)')(RA(I),I=1,6)
          ENDDO
      ENDDO
11    CLOSE(7)
      CLOSE(8)
      RETURN
      END

      SUBROUTINE WRITE_DATA(UNIT,VARC,VAR,LEV)
      DIMENSION VAR(LEV)
      INTEGER LEV,UNIT
      CHARACTER VARC*2
      WRITE(UNIT,'(A2,60F10.5)')VARC,(VAR(K),K=1,LEV)
C      WRITE(*,'(5F10.5)')(VAR(K),K=1,LEV)
      RETURN
      END

       SUBROUTINE CAL_RH_TD(T,RH,LEV,TD)
C     THIS SUBROUTINE CALCULATE DEW POINT BY USING TEMPERATURE AND RELATIVE HUMIDITY
C     INPUT   T ========> TEMPERATURE (CELSIUS DEGREE)
C             RH =======> RELATIVE HUMIDITY  (E/ES)
C             LEV ======> LEVEL
C     OUTPUT  TD =======> DEW POINT  (CELSIUS DEGREE)
      DIMENSION T(LEV),RH(LEV),TD(LEV)
      INTEGER LEV
      DO K=1,LEV
          IF(RH(K).GT.1.0)RH(K)=1.0
          IF(RH(K).LT.0.05)RH(K)=0.05
          RMA=ALOG(RH(K))*461.*273.15/4.18/597./1000.+T(K)/
     1      (273.15+T(K))
          TD(K)=273.15*RMA/(1-RMA)
      ENDDO
      RETURN
      END
        SUBROUTINE CAL_CAPE(P,T,TD,LEV,CAPE,PLFC,PE,CIN,LI,QL,PC,TC)
C     THIS SUBROUTINE CALCULATE CONVECTIVE AVAILABLE POTENTIAL ENERGY
C      AND IT IS USED WITH CAL_TC_PC,CAL_TA,CAL_TE,CAL_ENGY,
C      IN WHICH SUBROUTINE CAL_TA AND CAL_TC_PC ARE USED WITH CAL_THSE
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT  CAPE ===> EQUIVALENT POTENTIAL TEMPERATURE   (J/KG)
C            PLFC ===> FREE CONVECTION LEVEL
C            PE =====> EMBALANCE LEVEL
C            CIN ====> CONVECTIVE INHIBITION
C            LI =====> LIFTING INDEX            
C    TEMP ARRAY
C                TE =====> ENVIRONMENT TEMPERATURE   (KELVIN)
C            TA =====> PARCEL TEMPERATURE        (KELVIN)
C            TC
      DIMENSION P(LEV),T(LEV),TD(LEV),QL(LEV)
      DIMENSION TE(120),TA(120)
      REAL CAPE,TC,PC,THSE,PC0,TC0,TD0,THSE0
      REAL PLFC,PE,LI,CIN
      INTEGER LEV
C    ASSUMING THE BOTTOM POINT AS THE LIFTING LAYER
      THSE=220.
      PC=P(1)
      TC=T(1)
      P0=P(1)
      T0=T(1)
      TD0=TD(1)
      CALL CAL_TC_PC(P0,T0,TD0,PC0,TC0,THSE0)
      PC=PC0
      TC=TC0
      CALL CAL_TE(P,T,TE,LEV)
C      WRITE(*,*)'PC =',PC,' TC =',TC,' THSE =',THSE0
      CALL CAL_TA(P0,T0,PC,TC,THSE0,TA)
      WRITE(*,'(10F8.3)')(TE(I),I=1,70)
      WRITE(*,'(10F8.3)')(TA(I),I=1,70)
C      PAUSE 110
      CALL CAL_PLFC_AND_PE(P,TE,TA,LEV,PLFC,PE)
C      CALL CAL_ENGY(P0,TE,TA,CAPE)
      CALL CAL_VIR_ENGY(P,T,TD,LEV,TE,TA,CAPE)
      CALL CAL_CIN(P0,TE,TA,CIN)
      CALL CAL_LI(P0,TE,TA,LI)
      CALL CAL_QL(PC,TC,P,TA,LEV,QL)
      RETURN
      END

      SUBROUTINE CAL_SI0(T850,TD850,T500,SI)
C    CALCULATE SHAWALTER INDEX
C    INPUT T850 =====> TEMPERATURE AT 850HPA (CELSIUS DEGREE)
C          TD850 ====> DEW POINT AT 850HPA  (CELSIUS DEGREE)
C          T500 =====> TEMPERATURE AT 500HPA (CELSIUS DEGREE)
C    OUTPUT
C          SI  ======> SHAWALTER INDEX  (CELSIUS DEGREE)
      REAL SI,PC,TC,TH_0
      CALL CAL_TC_PC(850.,T850,TD850,PC,TC,TH_0)
      CALL CAL_TA_FOR_SI(T850,PC,TC,TH_0,T500,SI)
      RETURN
      END

      SUBROUTINE CAL_TC_PC(P,T,TD,PC,TC,THSE)
C    CALCULATE LIFTING CONDENSATION LAYER AND LIFTING CONDENSATION TEMPERATURE
C    AND ALSO PSEUDO-POTENTIAL TEMPERATURE AT CONDESATION POINT
C    INPUT P ========> LIFTING PRESSURE (HPA)
C          T ========> LIFTING TEMPERATURE (CELSIUS DEGREE)
C          TD =======> DEW POINT (CELSIUS DEGREE)
C    OUTPUT
C          PC =======> CONDESATION PRESSURE
C          TC =======> CONDENSATION TEMPERATURE
C          THSE =====> PSEUSO EQUIVALENT POTENTION TEMPERATURE (KELVIN)
      REAL P,T,TD,PC,TC,THSE
      TC=T-(T-TD)*0.976/(0.976-0.000833*(237.3+TD)**2/(273.15+TD))
      PC=P*((273.15+TC)/(273.15+T))**3.5
      CALL CAL_THSE0(PC,TC,THSE)
      RETURN
      END

      SUBROUTINE CAL_THSE1(P,T,THSE)
C      本子程序考虑 0 度以下，假相当位温计算公式不同

C     INPUT P =======> PRESSURE
C           T =======> TEMPERATURE
C     OUTPUT 
C           THSE ====> PSEUDO EQUIVALENT POTENTIAL TEMPERATURE
      REAL E,A,B,P,T,THSE
      IF(T.GE.0)A=7.5*T/(237.3+T)
      IF(T.LT.0)A=9.5*T/(265.5+T)
      E=6.11*10.**A
      B=1000.0/(P-E)
      IF(T.GE.0) THSE=(273.16+T)*B**0.286*EXP((597.3+0.5*T)*0.622*E
     1      /(0.24*(273.16+T)*(P-E)))
      IF(T.LT.0) THSE=(273.16+T)*B**0.286*EXP((667.0+0.5*T)*0.622*E
     1      /(0.24*(273.16+T)*(P-E)))
      RETURN
      END

      SUBROUTINE CAL_THSE0(P,T,THSE)
C      计算假相当位温，不考虑 0 度以下假相当位温计算公式的差别

C     INPUT P =======> PRESSURE
C           T =======> TEMPERATURE
C     OUTPUT 
C           THSE ====> PSEUDO EQUIVALENT POTENTIAL TEMPERATURE
      REAL E,A,B,P,T,THSE
      A=7.5*T/(237.3+T)
      E=6.11*10.**A
      B=1000.0/(P-E)
      THSE=(273.16+T)*B**0.286*EXP((2.5*10.**6-2368.*T)
     1*0.622*E/(1004.*(273.16+T)*(P-E)))
      RETURN
      END

      SUBROUTINE CAL_TA(PP0,TP0,PC,TC,THSE,TA)
C      计算状态曲线各整层高度（间隔 10HPA）上的温度

C     INPUT PP0 ======> PRESSURE AT THE LIFTING LEVEL (HPA)
C           TP0 ======> TEMPERATURE AT THE LIFTING LEVEL (CELSIUS DEGREE)
C           PC =======> LIFTING CONDENSATION LEVEL (HPA)
C           TC =======> LIFTING CONDESATION TEMPERATURE (CELSIUS DEGREE)
C           THESE ====> PSEUDO EQUIVALENT POTENTIAL TEMPERATURE (KELVIN)
C     OUTPUT ARRAY
C           TA =======> PARCEL TEMPERATURE AT EACH LEVEL (CELSIUS DEGREE)
      REAL PP0,TP0,PC,TC,THSE,THSE1,THSE2,P_L,T1,T2,RD,CPD
      INTEGER KPC
      DIMENSION TA(110)
      RD=287.04
      CPD=1005.
      DO K=1,110
            TA(K)=-130.
      ENDDO
C ************** DRY ADIABATIC PROCESS BEGIN
      KPC=INT(PP0/10)-INT(PC/10)
      P_L=INT(PP0/10.)*10.
      THD=(TP0+273.15)*(1000./PP0)**(RD/CPD)
      DO I=1,KPC
            TA(I)=THD*(P_L/1000.)**(RD/CPD)-273.15
            P_L=P_L-10.
      ENDDO
C      DO I=1,KPC
C           TA(I)=TP0+(TC-TP0)/ALOG(PC/PP0)*ALOG(P_L/PP0)
C            P_L=P_L-10.
C      ENDDO
C ************** DRY ADIABATIC PROCESS END
      IN=KPC
      T1=TC
      P1=PC
C      WRITE(*,*)' IN = ',IN,' TC = ',TC,' PC = ',PC
      KZERO=1
C      KZERO IS A JUDGEMENT FOR TEMPERATURE
C *************** MOIST ADIABATIC PROCESS BEGIN
10    CALL CAL_THSE0(P_L,T1,THSE1)
      CALL CAL_THSE0(P_L,T1-0.1,THSE2)
      T2=T1-(THSE1-THSE)/(THSE1-THSE2)*0.1
C      NOTICE DIFFERENCE OF THSE WHEN TEMPERERTURE OBOVE OR BELEOW ZERO
      IF(T2.LT.0.AND.KZERO.EQ.1)THEN
      CALL CAL_THSE0(P_L,T2,THSE)
      KZERO=0
      ENDIF       
      IN=IN+1
      TA(IN)=T2
      T1=T2
      P_L=P_L-10.
C      WRITE(*,*)' T1 = ',T1
C      PAUSE
      IF(P_L.GE.90.)GOTO 10
C *************** MOIST ADIABATIC PROCESS END
      RETURN
      END

      SUBROUTINE CAL_SI(P,T,TD,LEV,SI)
C    THIS SUBROUTINE CALCULATE SHAWALTER INDEX
C    NOTICE IT IS SUPPORTED BY CAL_TE,CAL_TDE,AND CAL_SI0 AND THERE
C    CORRESPONDING SUBROUTINE
C    CALCULATE SHAWALTER INDEX
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C    OUTPUT SI ======> SHAWALTER INDEX
      INTEGER LEV
      REAL SI
      DIMENSION P(LEV),T(LEV),TD(LEV)
      DIMENSION TES(110),TDE(110)
      CALL CAL_TE(P,T,TES,LEV)
      DO K=1,LEV
          IF(TD(K).GT.90)   TD(K)=T(K)-30.
      ENDDO
      CALL CAL_TDE(P,TD,TDE,LEV)
      NP850=1
      NP500=1
      P_L=INT(P(1)/10.)*10.
1     CONTINUE
      IF(P_L.GT.850.)THEN
          P_L=P_L-10.
            NP850=NP850+1
          GOTO 1
      ENDIF
      NP500=NP850+35
      T850=TES(NP850)
      TD850=TDE(NP850)
      T500=TES(NP500)
      CALL CAL_SI0(T850,TD850,T500,SI)
      RETURN
      END

      SUBROUTINE CAL_TA_FOR_SI(T850,PC,TC,THSE,T500,SI)
C      计算850HPA对应点绝热抬升所对应的状态曲线及沙氏指数
      REAL PP0,TP0,PC,TC,THSE,THSE1,THSE2,P_L,T1,T2
      INTEGER KPC
      PP0=850.
      TP0=T850
      KPC=INT(PP0/10)-INT(PC/10)
      P_L=INT(PP0/10.)*10.
C     干绝热部分

      DO I=1,KPC
            T1=TP0+(TC-TP0)/ALOG(PC/PP0)*ALOG(P_L/PP0)
            P_L=P_L-10.
            IF(I.EQ.35)THEN
                  SI=T500-T1
                  RETURN
            ENDIF
      ENDDO
C     湿绝热部分

      T1=TC
      P1=PC
      KZERO=1
C      KZERO IS A JUDGEMENT FOR TEMPERATURE
10      CALL CAL_THSE0(P_L,T1,THSE1)
      CALL CAL_THSE0(P_L,T1-0.1,THSE2)
      T2=T1-(THSE1-THSE)/(THSE1-THSE2)*0.1
C      NOTICE DIFFERENCE OF THSE WHEN TEMPERERTURE ABOVE OR BELEOW ZERO DEGREE
      IF(T2.LT.0.AND.KZERO.EQ.1)THEN
      CALL CAL_THSE0(P_L,T2,THSE)
      KZERO=0
      ENDIF
      T1=T2
      P_L=P_L-10.
      IF(P_L.GE.500.)GOTO 10
      SI=T500-T2
      RETURN
      END
      
      SUBROUTINE CAL_TE(PP,TP,TE,KMAX)
C      ********** CALCULATE TEMPERTURE AT EACH LEVEL *********
C      计算层结每隔 10HPA 的温度（对数气压线性内插）
C      PARAMETER (KMAX=30)
      DIMENSION PP(KMAX),TP(KMAX),TE(110)
      REAL P_L
      INTEGER I_NL,IN
      DO K=1,110
            TE(K)=-60.
      ENDDO
C      WRITE(*,*)(PP(K),K=1,KMAX)
C      WRITE(*,*)(TP(K),K=1,KMAX)
c      PAUSE 31
      I_NL=KMAX
      I=1
      IN=1
101            P_L=INT(PP(I)/10.)*10.0
            P_T=INT(PP(I+1)/10.)*10.+10.
10      TE(IN)=TP(I)+(TP(I)-TP(I+1))/ALOG(PP(I)/PP(I+1))*ALOG(P_L/PP(I))
      IN=IN+1
      P_L=P_L-10.
      IF(P_T.GT.PP(I))IN=IN-1
      IF(P_L.GE.P_T)GOTO 10
            I=I+1
C            IF(PP(I+1).GT.50.OR.(I+1).LT.(I_NL))THEN
            IF((I+1).LE.(I_NL))THEN
                  GOTO 101
            ELSE
                  GOTO 20
            ENDIF
20      CONTINUE
      IF((PP(I_NL).GE.100.).AND.
     1            (INT(PP(I_NL)/10.)*10.).EQ.(PP(I_NL)))
     2      TE(IN)=TP(I_NL)
      RETURN
      END

      SUBROUTINE CAL_TDE(PP,TP,TE,KMAX)
C      ********** CALCULATE TEMPERTURE AT EACH LEVEL *********
C      计算层结每隔 10HPA 的温度（对数气压线性内插）
C      PARAMETER (KMAX=30)
      DIMENSION PP(KMAX),TP(KMAX),TE(110)
      REAL P_L
      INTEGER I_NL,IN
      DO K=1,110
            TE(K)=-120.
      ENDDO
      I_NL=KMAX
      I=1
      IN=1
101            P_L=INT(PP(I)/10.)*10.0
            P_T=INT(PP(I+1)/10.)*10.+10.
10    TE(IN)=TP(I)+(TP(I)-TP(I+1))/ALOG(PP(I)/PP(I+1))*ALOG(P_L/PP(I))
      IN=IN+1
      P_L=P_L-10.
      IF(P_T.GT.PP(I))IN=IN-1
      IF(P_L.GE.P_T)GOTO 10
      I=I+1
C     IF(PP(I+1).GT.50.OR.(I+1).LT.(I_NL))THEN
      IF((I+1).LE.(I_NL))THEN
        GOTO 101
      ELSE
        GOTO 20
      ENDIF
20    CONTINUE

      IF((PP(I_NL).GE.100.).AND.
     1            (INT(PP(I_NL)/10.)*10.).EQ.(PP(I_NL)))
     2  TE(IN)=TP(I_NL)
      RETURN
      END

      SUBROUTINE CAL_CIN(PP0,TE,TA,CIN)
C      CALCULATE CONVECTIVE INHIBITION
C     INPUT PP0 ======> BOTTOM PRESSURE
C           TE =======> ENVIRONMENTAL TEMPERATURE
C           TA =======> PARCEL TEMPERATURE
C     OUTPUT
C           CIN ======> CONVECTIVE INHIBITION (J/KG)
      REAL PP0,R,P_L,CIN
      DIMENSION TE(110),TA(110)
      CIN=0.
      R=287.04
      P_L=INT(PP0/10.)*10
      P_T=P_L-10.
      IF(TA(1).GT.TE(1).OR.TA(2).GT.TE(2))THEN
          CIN=0.0
          RETURN
      ELSE
          DO K=1,90
                IF(TA(K).LE.TE(K).AND.TA(K+1).LT.TE(K+1).AND.P_L.GE.300)
     1            THEN
                  CIN=CIN-R*(TA(K)+TA(K+1)-TE(K)-TE(K+1))
     1               /2.*ALOG(P_L/P_T)
              ELSEIF(TA(K+2).GT.TE(K+2))THEN
                  RETURN
              ELSEIF(P_L.LT.300.)THEN
                  CIN=99999.9
                  RETURN
              ENDIF
                P_L=P_T
                P_T=P_T-10.
          ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_ENGY(PP0,TE,TA,CAPE)
C      计算对流有效位能（能量积分、不考虑虚位温）
      REAL PP0,CAPE,R,P_L
      DIMENSION TE(110),TA(110)
            CAPE=0.
      R=287.3
      P_L=INT(PP0/10.)*10.
      P_T=P_L-10.

      DO K=1,100
            IF(TA(K).GT.TE(K).AND.TA(K+1).GT.TE(K+1))
     1      CAPE=CAPE+R*(TA(K)+TA(K+1)-TE(K)-TE(K+1))/2.*ALOG(P_L/P_T)
            P_L=P_T
            P_T=P_T-10.
            IF(P_T.LT.150)GOTO 10
C      WRITE(*,*)TA(K),TE(K)
C      PAUSE
      ENDDO
10      CONTINUE
      RETURN
      END

      SUBROUTINE CAL_QL(PC,TC,P,TA,LEV,QL)
C      计算绝热比含水量
      REAL PC,TC,QL,QS0
      DIMENSION TA(110),P(LEV),QL(LEV)
      INTEGER LEV
      QS0=0.622*6.112*10**(7.5*TC/(237.3+TC))/PC
      DO K=1,LEV
            IF(P(K).GE.PC)THEN
                QL(K)=0.0
            ELSE
                IP=INT(P(1)/10.)-INT(P(K)/10.)+1
              PL=INT(P(K)/10.)*10.
              QL(K)=QS0-0.622*6.112*10**(7.5*TA(IP)/(237.3+TA(IP)))/PL
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE CAL_LI(PP0,TE,TA,LI)
C      CALCULATE LIFTING INDEX
C     INPUT PP0 ========> BOTTON PRESSURE
C           TE =========> ENVIRONMENTAL TEMPERATURE
C           TA =========> PARCEL TEMPERATURE 
C     OUTPUT
C           LI =========> LIFTING INDEX
      REAL PP0,LI,P_L
      DIMENSION TE(110),TA(110)
       NP=1
       P_L=INT(PP0/10.)*10.
1     CONTINUE
      IF(P_L.GT.500.)THEN
          P_L=P_L-10.
          NP=NP+1
      GOTO 1
      ENDIF
      LI=TE(NP)-TA(NP)
C      WRITE(*,*)' NP === ',NP
C      WRITE(*,*)' LI === ',LI
      RETURN
      END

      SUBROUTINE CAL_VIR_ENGY(P,T,TD,LEV,TE,TA,CAPEV)
C      计算对流有效位能（能量积分、考虑虚位温）
C    NOTICE IT IS SUPPORTED BY CAL_TDE,CAL_Q AND CAL_QS
C    INPUT   T ======> TEMPERATURE (CELSIUS DEGREE)
C            TD =====> DEW POINT (CELSIUS DEGREE)
C            P ======> PRESSURE   (HPA)
C            LEV ====> LEVEL
C            TE =====> ENVIRONMENT TEMPERATURE (CELSIUS DEGREE)
C            TA =====> PARCEL TEMPERATURE (CELSIUS DEGREE)
C    OUTPUT
C            CAPEV ==> CONVECTIVE AVAILABLE POTENTIAL ENERGY (J/KG)
C    TEMP ARRAY
C            PDS ====> PRESSURE (HPA)
C            QE =====> SPECIFIC HUMICITY
C            QSE ====> SATURATE SPECIFIC HUMIDITY
C            TDE ====> DEW POINT (CELSIUS DEGREE)
      REAL CAPEV,R,P_L
      DIMENSION TE(110),TA(110),TDE(110),P(LEV),T(LEV),TD(LEV)
      DIMENSION PDS(110),QE(110),QSE(110)
      DO K=1,LEV
          IF(TD(K).GT.90)TD(K)=T(K)-30.
      ENDDO
      CALL CAL_TDE(P,TD,TDE,LEV)
C      WRITE(*,'(10F8.2)')(TDE(K),K=1,90)
      NP=1
      P_L=INT(P(1)/10.)*10.
1     CONTINUE
      IF(P_L.GT.100.)THEN
          PDS(NP)=P_L
          P_L=P_L-10.
          NP=NP+1
      GOTO 1
      ENDIF
      NP=NP-1
      CALL CAL_Q(TDE,PDS,NP,QE)
      CALL CAL_QS(TA,PDS,NP,QSE)
C      WRITE(*,*)' NP = ',NP
C      WRITE(*,'(10F8.4)')(QE(I),I=80,110)
C      WRITE(*,'(10F8.4)')(QSE(I),I=80,110)
C      PAUSE 110
      CAPEV=0.
      R=287.04
      P_L=INT(P(1)/10.)*10
      P_T=P_L-10.
      DO K=1,NP
            IF(TA(K).GT.TE(K).AND.TA(K+1).GT.TE(K+1))
     1      CAPEV=CAPEV+R*((1+0.61*QSE(K))*(TA(K)+273.15)
     2           +(1+0.61*QSE(K+1))*(TA(K+1)+273.15)
     3           -(1+0.61*QE(K))*(TE(K)+273.15)
     4           -(1+0.61*QE(K+1))*(TE(K+1)+273.15))/2.*ALOG(P_L/P_T)
C     1      CAPEV=CAPEV+R*(TA(K)+TA(K+1)-TE(K)-TE(K+1))/2.*ALOG(P_L/P_T)
            P_L=P_T
            P_T=P_T-10.
            IF(P_T.LT.P(LEV))GOTO 10
      ENDDO
10      CONTINUE
      RETURN
      END

      SUBROUTINE CAL_PLFC_AND_PE(P,TE,TA,LEV,PLFC,PE)
C     THIS SUBROUTINE CALCULATE LAYER OF FREE CONVECTION AND EMBALANCE LAYER
C     NOTICE THIS SUBROUTINE IS NOT INDEPENDENT!!!
C       INPUT  TE =====> ENVIRONMENT TEMPERATURE   (KELVIN)
C            TA =====> PARCEL TEMPERATURE        (KELVIN)
C            P  =====> PRESSURE        (HPA)
C            LEV ====> NUMBER OF LAYER
C     OUTPUT PLFC ===> FREE CONVECTION LAYER  (HPA)
C            PE =====> EMBALANCE PRESSURE  (HPA)
      REAL PLFC,PE
      DIMENSION TE(110),TA(110),P(LEV)
      INTEGER LEV,KLFC
      REAL DELT1,DELT2
C      WRITE(*,'(10F8.2)')(TA(K),K=1,10),(TE(K),K=1,10)
      KLFC=0
C      IF(TA(1).GE.TE(1).AND.TA(2).GT.TE(2))THEN
      IF(TA(1).GT.TE(1).OR.(TA(1).EQ.TE(1).AND.TA(2).GT.TE(2)))THEN
          PLFC=P(1)
      ELSEIF((TA(1).EQ.TE(1).AND.TA(2).LE.TE(2)).OR.TA(2).LE.TE(2))THEN
          P_L=INT(P(1)/10.)*10
          P_T=P_L-10
          DO K=2,90
              IF(TA(K).GT.TE(K))THEN
                  KLFC=K
                  P_B=P_T+10
                  DELT1=TE(K-1)-TA(K-1)
                  DELT2=TA(K)-TE(K)
                  PLFC=(P_B*DELT2+P_T*DELT1)/(DELT1+DELT2)
                  GOTO 290
               ENDIF
          P_T=P_T-10
          ENDDO
         IF(KLFC.EQ.0)PLFC=99999.9
      ENDIF
290      CONTINUE
       NB=INT(P(1)/10.)
      NT=INT((P(LEV)-0.001)/10.)+1
      NE=NB-NT+1
      KPE=NE
      P_T=INT(P(LEV)/10.)*10+10.
      IF(TE(NE).LT.TA(NE))THEN
          PE=99999.9
      ELSE
          DO K=NE-1,1,-1
              IF(TE(K).LT.TA(K))THEN
                  KPE=K
                  P_B=P_T+10
C    LYD K-1 ==>K+1
                  DELT1=TE(K+1)-TA(K+1)
                  DELT2=TA(K)-TE(K)
                  PE=(P_B*DELT2+P_T*DELT1)/(DELT1+DELT2)
                  GOTO 390
              ENDIF
          P_T=P_T+10
          ENDDO
          IF(KPE.EQ.NE)PE=99999.9
      ENDIF
390   CONTINUE
      RETURN
      END
      SUBROUTINE CAL_SRH(H0,P,H,WD,WS,LEV,RSH)
C    THIS SUBROUTINE CALCULATE STORM_RELATIVE ENVIRONMENTAL HELICITY(0-H0)
C    OUTPUT  RSH  ====> STORM_RELATIVE ENVIRONMENTAL HELICITY (M2/S2)
C    TEMP ARRAY
C            C  =====>STORM MOVE SPEED  
C            DC  =====>STORM MOVE DIRECTION 

      DIMENSION P(LEV),WD(LEV),WS(LEV),U(LEV),V(LEV),P1(LEV),H1(LEV),
     *          H(LEV)
      REAL H0
      PI=3.14159
      RSH=99999.9
      N1=0
      DO I=1,LEV 
        IF(P(I).GE.700.AND.WD(I).NE.99999.AND.WS(I).NE.99999) THEN
           N1=N1+1
        END IF
      END DO

      IF(N1.LE.2)RETURN

      M=0
      DO I=1,LEV 
         IF(WD(I).NE.99999.AND.WS(I).NE.99999)THEN
             M=M+1
             U(M)=WS(I)*SIN((WD(I)-180)*PI/180)
             V(M)=WS(I)*COS((WD(I)-180)*PI/180)
             P1(M)=P(I)
             H1(M)=H(I)
         END IF
      END DO

      LEV1=M
      CX=0
      CY=0
      N=0
      DO I=1,LEV1
          IF(P1(I).LE.850.AND.P1(I).GE.400)THEN
              CX=CX+U(I)
              CY=CY+V(I)
              N=N+1
            END IF
      END DO

      CX=CX/N
      CY=CY/N
      
      IF(ABS(CY).LT.0.00001)THEN
          IF(CX.GT.0) DC=2*PI-PI/2
          IF(CX.LT.0) DC=PI/2
      ELSE IF(ABS(CX).LT.0.00001.AND.CY.NE.0)THEN
          IF(CY.GT.0) DC=PI
          IF(CY.LT.0) DC=0
      ELSE
           DC=ATAN(CX/CY)+PI
      END IF

        DC=DC+PI*40./180.
        C=0.75*SQRT(CX**2+CY**2)
        CX=C*SIN((DC-180)*PI/180)
        CY=C*COS((DC-180)*PI/180)

C-------------CALCULATE U0,V0 ON H0---------
      H0=H0+H(1)
      CALL CAL_H0_UV(H0,H1,U,V,LEV1,U0,V0)

      RSH=0
        DO I=2,LEV1
           IF(H1(I).LT.H0)THEN
               KK=I 
             RSH=RSH+(U(I)-CX)*(V(I-1)-CY)-(U(I-1)-CX)*(V(I)-CY)
           END IF
        END DO
      RSH=RSH+(U0-CX)*(V(KK)-CY)-(U(KK)-CX)*(V0-CY)
      RETURN
      END 

C-----------------
      SUBROUTINE CAL_H0_UV(H0,H,U,V,LEV,U0,V0)
C    THIS SUBROUTINE CALCULATE WIND SPEED(U0,V0) ON GIVEN HIGHT H0
      DIMENSION H(LEV),U(LEV),V(LEV)
      INTEGER LEV
C    LYD LEV ==>LEV-1
      DO K=1,LEV-1
         IF(H(K+1).GT.H0.AND.H(K).LE.H0)THEN
            CC=(H(K+1)-H0)/(H(K+1)-H(K))
            U0=U(K+1)+CC*(U(K)-U(K+1))
            V0=V(K+1)+CC*(V(K)-V(K+1))
         END IF
      END DO

      RETURN
      END 

C**************************************************

      SUBROUTINE CAL_BRN(P,H,T,TD,WD,WS,LEV,CAPE,BRN)
C     THIS SUBROUTINE CALCULATE BULK RICHARSON NUMBER
C     OUTPUT SHR  
C     TEMP ARRAY
C             R  ====> THE AIR DENSITY
C
      DIMENSION P(LEV),WD(LEV),WS(LEV),H(LEV),U(LEV),V(LEV),T(LEV),
     *          TD(LEV),R(LEV),H1(LEV),T1(LEV),TD1(LEV),P1(LEV),
     *          WS1(LEV)
      INTEGER LEV,LEV1
      REAL CAPE
      PI=3.14159
      M=0
      DO I=1,LEV
        IF(WD(I).NE.99999.AND.WS(I).NE.99999)THEN
                 M=M+1
             U(M)=WS(I)*SIN((WD(I)-180)*PI/180)
             V(M)=WS(I)*COS((WD(I)-180)*PI/180)
             P1(M)=P(I)
               H1(M)=H(I)
             T1(M)=T(I)
             TD1(M)=TD(I)
             WS1(M)=WS(I)
             R(M)=P1(M)/(287*(273.15+T1(M))) 
        END IF
      END DO

      LEV1=M
C-----------CALCULATE U,V,P,T,R ON 6KM HEIGHT------
      H0=6000+H(1) 
      CALL CAL_H0_UVPT(H0,H1,T1,TD1,P1,U,V,LEV1,P6,T6,TD6,U6,V6)
      R6=P6/(287*(273.15+T6))
      WS6=SQRT(U6**2+V6**2)
C---------------------------------------------------

      A=R(2)*WS1(2)*(H1(2)-H1(1)) 
      B=R(2)*(H1(2)-H1(1))
      DO I=3,LEV1
         IF(H1(I).LT.H0)THEN
            KK=I
            A=A+R(I)*WS1(I)*(H1(I)-H1(I-1))
            B=B+R(I)*(H1(I)-H1(I-1))
         END IF
      END DO

      A=A+R6*WS6*(H0-H1(KK))
      B=B+R6*(H0-H1(KK))

C-------------CALCULATE U,V ON 0.5KM HEIGHT---------
      H0=500+H(1)
      CALL CAL_H0_UV(H0,H1,U,V,LEV1,U5,V5)
C---------------------------------------------------

      SHR1=A/B-SQRT((U5+U(1))**2+(V5+V(1))**2)
      IF(ABS(SHR1).GT.1.0)THEN
          BRN=CAPE/(2.*SHR1**2)
      ELSE
          BRN=CAPE/2.
      ENDIF
      RETURN
      END

      SUBROUTINE CAL_H0_UVPT(H0,H,T,TD,P,U,V,LEV,P0,T0,TD0,U0,V0)
C    THIS SUBROUTINE CALCULATE U0,V0,P0,T0,TD0 ON GIVEN HIGHT H0
      DIMENSION H(LEV),T(LEV),TD(LEV),P(LEV),U(LEV),V(LEV)
      INTEGER LEV 
c     LYD LEV=> LEV-1
      DO K=1,LEV-1
               IF(H(K+1).GT.H0.AND.H(K).LE.H0)THEN
              CC=(H(K+1)-H0)/(H(K+1)-H(K))
              T0=T(K+1)+CC*(T(K)-T(K+1))             !!DT~DH
              TD0=T0-(T(K+1)-TD(K+1))
                  C1=(T(K)+T(K+1)+2*273.15)/(T0+T(K+1)+2*273.15)
                  C2=C1*CC
              P0=P(K+1)*(P(K)/P(K+1))**C2                !! DLNP~H/TP
              U0=U(K+1)+CC*(U(K)-U(K+1))               !!DU~DH
              V0=V(K+1)+CC*(V(K)-V(K+1))               !!DV~DH
         END IF
      END DO
      RETURN
      END 


C******************************************************
      SUBROUTINE CAL_SSI(H0,H,WD,WS,LEV,CAPE,SSI)
C     THIS SUBROUTINE CALCULATE INDEX OF STORM
C     OUTPUT SSI 
C     TEMP ARRAY
C             R  ====> THE AIR DENSITY

      DIMENSION WD(LEV),WS(LEV),H(LEV),U(LEV),V(LEV),
     *          H1(LEV)
      INTEGER LEV,LEV1
      REAL SHR2
      PI=3.14159
      M=0
      DO I=1,LEV
        IF(WD(I).NE.99999.AND.WS(I).NE.99999.)THEN
                 M=M+1
             U(M)=WS(I)*SIN((WD(I)-180)*PI/180)
             V(M)=WS(I)*COS((WD(I)-180)*PI/180)
               H1(M)=H(I)
        END IF
      END DO

      LEV1=M

      H01=H0+H(1) 
      CALL CAL_H0_UV(H01,H1,U,V,LEV1,U0,V0)
c      SHR2=SQRT((U0-U(1))**2+(V0-V(1))**2)/(H0/1000)
      SHR2=SQRT((U0-U(1))**2+(V0-V(1))**2)/H0
      SSI=100*(2+0.276*ALOG(SHR2)+2.011*CAPE*0.0001)
      RETURN
      END

      SUBROUTINE CAL_SWISS00(P,T,TD,H01,H02,H,WD,WS,LEV,SI,SWISS00)
C     THIS SUROUTIINE CALCULATE SWISS STORM INDEX AT 00 UTC
      DIMENSION P(LEV),T(LEV),TD(LEV),H(LEV),WD(LEV),WS(LEV)
      DIMENSION TE(110),TDE(110)
      INTEGER LEV
      REAL H01,H02,SI,SWISS00,WSH3_6,TTD600
      SWISS=99999.9
      IF(P(1).LT.850.)THEN
            RETURN
      ELSE
      CALL CAL_SHR3(H01,H02,H,WD,WS,LEV,WSH3_6)
C     WSH3_6 (M/S PER 3KM)
      WSH3_6=WSH3_6*3
          CALL CAL_TE(P,T,TE,LEV)
          CALL CAL_TDE(P,TD,TDE,LEV)
          P_L=INT(P(1)/10.)*10.
          DO K=1,50
              IF(P_L.EQ.600)THEN
                  TD600=TDE(K)
                  T600=TE(K)
                  TTD600=T600-TD600
                  SWISS00=SI+0.4*WSH3_6+0.1*TTD600
              RETURN
              ENDIF
              P_L=P_L-10.
          ENDDO
      ENDIF      
      RETURN
      END

      SUBROUTINE CAL_SWISS12(P,T,TD,H0,H,WD,WS,LEV,LI,SWISS12)
C     THIS SUROUTIINE CALCULATE SWISS STORM INDEX AT 12 UTC
      DIMENSION P(LEV),T(LEV),TD(LEV),H(LEV),WD(LEV),WS(LEV)
      DIMENSION TE(110),TDE(110)
      INTEGER LEV
      REAL H0,LI,SWISS12,WSH0_3,TTD650
      SWISS=99999.9
      IF(P(1).LT.650.)THEN
            RETURN
      ELSE
      CALL CAL_SHR2(H0,H,WD,WS,LEV,WSH0_3)
C     WSH3_6 (M/S PER 3KM)
      WSH0_3=WSH0_3*3
          CALL CAL_TE(P,T,TE,LEV)
          CALL CAL_TDE(P,TD,TDE,LEV)
          P_L=INT(P(1)/10.)*10.
          DO K=1,50
              IF(P_L.EQ.650)THEN
                  TD650=TDE(K)
                  T650=TE(K)
                  TTD650=T650-TD650
                  SWISS12=LI+0.3*WSH0_3+0.3*TTD650
              RETURN
              ENDIF
              P_L=P_L-10.
          ENDDO
      ENDIF      
      RETURN
      END


C******************************************************
      SUBROUTINE CAL_SHR2(H0,H,WD,WS,LEV,SHR2)
C     THIS SUBROUTINE CALCULATE VERTICAL SHEAR(0-H0)
C     OUTPUT SHR (M/S PER 1KM)
C     TEMP ARRAY
C             R  ====> THE AIR DENSITY

      DIMENSION WD(LEV),WS(LEV),H(LEV),U(LEV),V(LEV),
     *          H1(LEV)
      INTEGER LEV,LEV1
      REAL SHR2
      PI=3.14159
      M=0
      DO I=1,LEV
        IF(WD(I).NE.99999.AND.WS(I).NE.99999.)THEN
                 M=M+1
             U(M)=WS(I)*SIN((WD(I)-180)*PI/180)
             V(M)=WS(I)*COS((WD(I)-180)*PI/180)
               H1(M)=H(I)
        END IF
      END DO

      LEV1=M

      H01=H0+H(1) 
      CALL CAL_H0_UV(H01,H1,U,V,LEV1,U0,V0)
      SHR2=SQRT((U0-U(1))**2+(V0-V(1))**2)/(H0/1000)
      RETURN
      END

C*****************************************************

      SUBROUTINE CAL_SHR3(H01,H02,H,WD,WS,LEV,SHR3)
C     THIS SUBROUTINE CALCULATE VERTICAL SHEAR FROM H01 TO H02
C     OUTPUT SHR3 (M/S PER 1KM )

      DIMENSION WD(LEV),WS(LEV),H(LEV),U(LEV),V(LEV),H1(LEV)
      INTEGER LEV,LEV1
      PI=3.14159
      M=0
      DO I=1,LEV
        IF(WD(I).NE.99999.AND.WS(I).NE.99999)THEN
                 M=M+1
             U(M)=WS(I)*SIN((WD(I)-180)*PI/180)
             V(M)=WS(I)*COS((WD(I)-180)*PI/180)
               H1(M)=H(I)
        END IF
      END DO

      LEV1=M

      H01=H01+H(1) 
      CALL CAL_H0_UV(H01,H1,U,V,LEV1,U1,V1)
      H02=H02+H(1) 
      CALL CAL_H0_UV(H02,H1,U,V,LEV1,U2,V2)
      SHR3=SQRT((U1-U2)**2+(V1-V2)**2)/(ABS(H02-H01)/1000)
      RETURN
      END


C***************************************************

      SUBROUTINE CAL_SWEAT(H,T,TD,WD,WS,LEV,SWEAT)
C    THIS SUBROUTINE CALCULATE SWEAT INDEX
C    OUTPUT  SWEAT   

      DIMENSION WD(LEV),WS(LEV),H(LEV),U(LEV),V(LEV),T(LEV),
     *          TD(LEV),H1(LEV),T1(LEV),TD1(LEV)
      INTEGER LEV,LEV1

      SWEAT=0.0

      PI=3.14159
      M=0
      DO I=1,LEV
        IF(WD(I).NE.99999.AND.WS(I).NE.99999)THEN
                 M=M+1
             U(M)=WS(I)*SIN((WD(I)-180)*PI/180)
             V(M)=WS(I)*COS((WD(I)-180)*PI/180)
               H1(M)=H(I)
             T1(M)=T(I)
               TD1(M)=TD(I)
        END IF
      END DO

      LEV1=M

      H0=900+H(1)
      CALL CAL_H_TUV(H0,H1,T1,TD1,U,V,LEV1,T01,TD01,DD01,FF01) 
      H0=5000+H(1)
      CALL CAL_H_TUV(H0,H1,T1,TD1,U,V,LEV1,T02,TD02,DD02,FF02) 

C     1M/S=1.94 SEA MILE/HOUR
      FF01=1.94*FF01
      FF02=1.94*FF02

      TT=T01+TD01-2*T02

      FA=0
      IF(DD01.GE.130.AND.DD01.LE.250.AND.
     *   DD02.GE.210.AND.DD02.LE.310.AND.
     *   DD02.GT.DD01.AND.
     *   FF01.GT.15.AND.
     *   FF02.GT.15)THEN
        FA=SIN(DD02-DD01)+0.2 
      END IF

       D=TD01
      IF(TD01.LT.0)  D=0
      IF(TT.LT.49)  TT=49

      SWEAT=12*D+20*(TT-49)+2*FF01+FF02+125*FA
      RETURN
      END

      SUBROUTINE CAL_H_TUV(H0,H,T,TD,U,V,LEV,T0,TD0,DD0,FF0)
C     THIS SUBROUTINE CALCULATE TEMPERATURE AND WIND SPEED ON GIVEN HEIGHT 
      DIMENSION H(LEV),T(LEV),TD(LEV),U(LEV),V(LEV)
      INTEGER LEV 
      PI=3.14159
C    LYD LEV ==>LEV-1
      DO K=1,LEV-1
         IF(H(K+1).GT.H0.AND.H(K).LE.H0)THEN
              CC=(H(K+1)-H0)/(H(K+1)-H(K))
              T0=T(K+1)+CC*(T(K)-T(K+1))             !!DT~DH
              TD0=T0-(T(K+1)-TD(K+1))
                  C1=(T(K)+T(K+1)+2*273.15)/(T0+T(K+1)+2*273.15)
                  C2=C1*CC
              U0=U(K+1)+CC*(U(K)-U(K+1))               !!DU~DH
              V0=V(K+1)+CC*(V(K)-V(K+1))               !!DV~DH
              FF0=SQRT(U0**2+V0**2)
C----------
              IF(ABS(V0).LT.0.00001)THEN
                IF(U0.GT.0) DD0=2*PI-PI/2
                IF(U0.LT.0) DD0=PI/2
              ELSE IF(ABS(U0).LT.0.00001.AND.V0.NE.0)THEN
                IF(V0.GT.0) DD0=PI
                IF(V0.LT.0) DD0=0
              ELSE
                DD0=ATAN(U0/V0)+PI
              END IF
C--------
         END IF
      END DO
      RETURN
      END 

C*********************************************

      SUBROUTINE CAL_WIN(P,H,T,TD,LEV,WINDEX)
C     THIS SUBROUTINE CALCULATE WINDEX 
C     OUTPUT  WINDEX (M/S,1KNOT=1.853KM/H=0.5147M/S)

      DIMENSION P(LEV),H(LEV),T(LEV),TD(LEV)

C-------------------------------------------
      WINDEX=0

      TM=0
C----------------
      CALL CAL_ZH(P,T,H,LEV,PM,HM)

      DO K=1,LEV
        IF(T(K).EQ.0)THEN
           TDM=TD(K+1)-T(K+1)
        ELSE
          IF((T(K).GT.0).AND.(T(K+1).LT.0))THEN        !!DT~DH
              TDM=TD(K+1)-T(K+1)
          END IF
        END IF
      END DO
C--------------------------
      HM=(HM-H(1))/1000
c    add by LYD Oct,27,2006
      if(hm.eq.0)return
c
      CALL CAL_R0(TDM,PM,QM)

      H0=H(1)+1000

      CALL CAL_H0_P0T0(H0,H,T,TD,P,LEV,P0,T0,TD0) 
      CALL CAL_R0(TD0,P0,Q0)

      SR=0
      N=0
      DO I=1,LEV
        IF((H(I)-H(1)).LT.1000)THEN
          N=N+1
          CALL CAL_R0(TD(I),P(I),RR)
            SR=SR+RR 
        END IF
      END DO
      QL=SR/N
      RQ=QL/12
      RR=(T(1)-TM)/HM
      QQ=RR**2-30+QL-2*QM
      IF(RR.LT.5.5) THEN
         WINDEX=0
      ELSEIF(QQ.LT.0)THEN
          WINDEX=0.
      ELSE
         WINDEX=5*(HM*RQ*QQ)**0.5
      END IF
      WINDEX=WINDEX*0.5147
C      WRITE(*,*)'GMA=',RR,' QL =',QL,' RQ=',RQ,' WINDEX=',WINDEX
C      PAUSE
      RETURN
      END 

      SUBROUTINE CAL_E0(TD,E)
C     THIS SUBROUTINE CALCULATE VAPOR PRESURE
C     NOTICE THE FORMULA ARE DIFFERENT WHEN DEW POINT TEMPERATURE 
C      IS ABOVE OR BELOW ZERO DEGREE
C    INPUT   TD =====> DEW POINT      (CELSIUS DEGREE)
C            LEV ====> LEVEL
C    OUTPUT  E  ====> VAPOR PRESSURE (HPA)

          IF(TD.GT.-130.AND.TD.LT.60.)THEN
              IF(TD.GE.0) E=6.112*10**(7.5*TD/(237.3+TD))
              IF(TD.LT.0) E=6.112*10**(9.5*TD/(265.5+TD))
          ELSE
              E=99999.9
          ENDIF
      RETURN
      END


      SUBROUTINE CAL_R0(TD,P,R)              
C     THIS SUBROUTINE CALCULATE MIXING RATION
C     NOTICE THIS SUBROUTINE IS SUPPORTED BY SUBROUTINE CAL_E0
C    INPUT   TD =====> DEW POINT      (CELSIUS DEGREE)
C            P  =====> PRESSURE   (HPA)
      CALL CAL_E0(TD,E)
            IF(E.LT.100.)THEN
C                R(K)=0.622*E(K)/(P(K)-E(K))
                R=622*E/(P-E)
          ELSE
              R=99999.9
          ENDIF
      RETURN
      END

      SUBROUTINE CAL_H0T0_P0(H0,T0,H,T,P,LEV,P0)
C     THIS SUBROUTINE CALCULATE PRESSURE ON GIVEN HEIGHT 
      DIMENSION H(LEV),T(LEV),P(LEV)
      INTEGER LEV 
      DO K=1,LEV
         IF(H(K+1).GT.H0.AND.H(K).LE.H0)THEN
              CC=(H(K+1)-H0)/(H(K+1)-H(K))
                  C1=(T(K)+T(K+1)+273.15*2)/(T0+T(K+1)+273.15*2)
                  C2=C1*CC
              P0=P(K+1)*(P(K)/P(K+1))**C2                !! DLNP~H/TP
         END IF
      END DO
      RETURN
      END 

      SUBROUTINE CAL_H0_P0T0(H0,H,T,TD,P,LEV,P0,T0,TD0)
C     THIS SUBROUTINE CALCULATE TEMPERATURE AND PRESSURE ON GIVEN HEIGHT 
      DIMENSION H(LEV),T(LEV),P(LEV),TD(LEV)
      INTEGER LEV 
      DO K=1,LEV
         IF(H(K+1).GT.H0.AND.H(K).LE.H0)THEN
              CC=(H(K+1)-H0)/(H(K+1)-H(K))
              T0=T(K+1)+CC*(T(K)-T(K+1))                !!DT~DH
              TD0=T0-(T(K+1)-TD(K+1))
                  C1=(T(K)+T(K+1)+2*273.15)/(T0+T(K+1)+2*273.15)
                  C2=C1*CC
              P0=P(K+1)*(P(K)/P(K+1))**C2                !! DLNP~H/TP
         END IF
      END DO
      RETURN
      END 

      SUBROUTINE CAL_ZH(P,T,H,LEV,PL,ZH)
C     THIS SUBROUTINE CALCULATE THE HEIGHT OF THE LEVEL OF TEMPERATURE IS ZERO.
      DIMENSION P(LEV),T(LEV),H(LEV)
      INTEGER LEV
      REAL ZH
      IF(T(1).LE.0.)THEN
          IF(T(1).EQ.0.)ZH=H(1)
          IF(T(1).LT.0.)ZH=99999.9
          RETURN
      ELSE
          DO K=1,LEV
              IF(T(K+1).LT.0)THEN
                  PL=EXP((T(K)*ALOG(P(K+1))-T(K+1)*ALOG(P(K)))   !!T~DLNP
     1               /(T(K)-T(K+1)))
                  ZH=H(K)+287.03*(2*273.15+T(K))/2./9.8*ALOG(P(K)/PL)
                  RETURN
                ENDIF
          ENDDO
      ENDIF
      RETURN
      END

C------------------------------------------------------------------------
      SUBROUTINE CAL_EHI(P,H,WD,WS,LEV,CAPE,EHI)
C     THIS SUBROUTINE CALCULATE ENERGETIC HELICITY INDEX
C     INPUT  P ========> PRESSURE
C            H ========> GEOPOTENTIAL HEIGHT
C            WS =======> WIND SPEED
C            WD =======> WIND DIRECTION
C            LEV ======> LEVEL
C            CAPE =====> CONVECTIVE AVAILABLE POTENTIAL ENERGY
C     OUTPUT EHI ======> ENERGETIC HELICITY INDEX
      DIMENSION P(LEV),H(LEV),WD(LEV),WS(LEV)
      INTEGER LEV
      REAL CAPE,EHI,SRH_EHI,H0
      H0=2000.
      CALL CAL_SRH(H0,P,H,WD,WS,LEV,SRH_EHI)
      IF(CAPE.GE.0.AND.SRH_EHI.GT.0.0)EHI=SRH_EHI*CAPE/160000.
      RETURN
      END
       SUBROUTINE FORM_TRAN_UV_WDS(U,V,LEV,WD,WS)
C    THIS SUBROUTINE DESIGNED BY LYD AT JUN,15,2003 
C    TO GET WIND SPEED AND DIRECTION
      DIMENSION U(LEV),V(LEV),WD(LEV),WS(LEV)
      INTEGER LEV
      DO K=1,LEV
          WS(K)=SQRT(U(K)**2+V(K)**2)
          IF(U(K).GT.0)THEN
                WD(K)=270-ATAN(V(K)/U(K))*57.2957796
          ELSEIF(U(K).LT.0)THEN
              WD(K)=90-ATAN(V(K)/U(K))*57.2957796
          ELSEIF(U(K).EQ.0)THEN
              IF(V(K).GT.0)WD(K)=180.
                IF(V(K).LE.0)WD(K)=0.
          ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE FORM_TRAN_WDS_UV(WD,WS,LEV,U,V)
C    THIS SUBROUTINE DESIGNED BY LYD AT JUN,15,2003 
C    TO GET WIND SPEED IN EAST-WEST AND NORTH-SOUTH DIRECTION
      DIMENSION U(LEV),V(LEV),WD(LEV),WS(LEV)
      INTEGER LEV
      DO K=1,LEV
          U(K)=WS(K)*SIN((360.-WD(K))/57.2957796)
          V(K)=WS(K)*COS((180.-WD(K))/57.2957796)
      ENDDO
      RETURN
      END

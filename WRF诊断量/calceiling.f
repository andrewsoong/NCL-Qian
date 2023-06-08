C NCLFORTSTART
      subroutine calceiling(qvp,qcw,qra,qci,qsn,tmk,prs,pf,
     &   clg,mkzh,mjj,mii)
       real qcw(mii,mjj,mkzh),qra(mii,mjj,mkzh),qvp(mii,mjj,mkzh),
     &   qci(mii,mjj,mkzh),qsn(mii,mjj,mkzh),tmk(mii,mjj,mkzh),
     &   prs(mii,mjj,mkzh),pf(mii,mjj,mkzh),clg(mii,mjj)
C NCLEND
c     This routine computes cloud ceiling.  It uses the same expressions
c     for beta as in routine viscalc.  However, rather than using a
c     constant beta, the routine integrates upward using variable beta,
c     until an upward beam of light would be extinguished to .02 times
c     its original intensity. qvp is passed in in kg/kg, whereas
c     the hydrometeor mixing ratios are passed in as g/kg.
C
C; set up constants
      celkel=273.15
      tice=celkel-10.
      coeflc=144.7
      coeflp=2.24
      coeffc=327.8
      coeffp=10.36
      exponlc=0.8800
      exponlp=0.7500
      exponfc=1.0000
      exponfp=0.7776
C;      CONST1= 3.912023 ;-LOG(0.02)
      thcon=0.02
      rhoice=917.
      rhowat=1000.
      rgas=287.04
      grav=9.81

      do 200 j=1,mjj-1
      do 200 i=1,mii-1
      exting=1.
      clg(i,j)=0.
      pup=pf(i,j,mkzh)
      do 150 k=mkzh,1,-1
         wmix=qvp(i,j,k)
            qprc=qra(i,j,k)+qsn(i,j,k)
            qcld=qcw(i,j,k)+qci(i,j,k)
            qrain=qra(i,j,k)
            qsnow=qsn(i,j,k)
            qclw=qcw(i,j,k)
            qclice=qci(i,j,k)
         tv=tmk(i,j,k)*(1.+0.608*qvp(i,j,k))
         rhoair=100.*prs(i,j,k)/(rgas*tv)
            vovermd=(1.+wmix)/rhoair+0.001*((qclw+qrain)/rhowat+
     &         (qclice+qsnow)/rhoice)
            conclc = qclw/vovermd
            conclp = qrain/vovermd
            concfc = qclice/vovermd
            concfp = qsnow/vovermd
         beta = coeffc * concfc ** exponfc + coeffp * concfp ** exponfp
     &        + coeflc * conclc ** exponlc + coeflp * conclp ** exponlp
     &        + 1.e-10
         pdn=pup
         if (k.eq.1) then
            pup=2.*prs(i,j,k)-pf(i,j,k)
         else
            pup=pf(i,j,k-1)
         endif
C m to km
         dz=.001*rgas*tv/grav*log(pdn/pup)   
         extingt=exting*exp(-beta*dz)
         if (extingt.le.thcon) then
            clg(i,j)=clg(i,j)-log(thcon/exting)/beta
            goto 160
         endif
         clg(i,j)=clg(i,j)+dz
         exting=extingt
  150 continue
  160 continue
C  km to m, max of 9,000m
      clg(i,j)=min(clg(i,j)*1000.,9000.)  
c
  200 continue
      return
      end subroutine calceiling
 

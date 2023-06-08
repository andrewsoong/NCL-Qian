  ;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
     load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
     load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
     load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
     load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
begin
;
; The WRF ARW input file.
; This needs to have a ".nc" appended, so just do it.
  fn_ctl  = "/public/home/wangbin/finalversion/WRFV3/run/wrfout_d03_2015-08-07_12:00:00.nc" ; define filename
  fn_NCR  = "/public/home/wangbin/finalversion/WRFV3-ncr/run/wrfout_d03_2015-08-07_12:00:00.nc" ; define filename
  fn_NIR  = "/public/home/wangbin/finalversion/WRFV3-nir/run/wrfout_d03_2015-08-07_12:00:00.nc" ; define filename
  fn_NSR  = "/public/home/wangbin/finalversion/WRFV3-nsr/run/wrfout_d03_2015-08-07_12:00:00.nc" ; define filename
  a_ctl = addfile(fn_ctl,"r")
  a_NCR = addfile(fn_NCR,"r")
  a_NIR = addfile(fn_NIR,"r")
  a_NSR = addfile(fn_NSR,"r")
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  times = wrf_user_getvar(a_ctl,"times",-1)  ; get all times in the file
  ntimes = dimsizes(times)         ; number of times in the file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     it=0
   print("Working on time: " + times(it) )
   print("There are " + ntimes + "times")
   ;; res@TimeLabel = times(it)   ; Set Valid time to use on plots   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  bottom_top=34 
    south_north=315
    west_east=294
    sn_s=110
    sn_e=309
    we_s=7
    we_e=245
    sn=200
    we=239


    grid_num=sn*we  
	step=60
   ; stimes=(ntimes-1)*step/3600 +1 ;number of hours,because number of steps*step=number of hours*3600
    stepnum=60*60/step
    grid_num=sn*we
    nspecies=4    ;number of sensitivity experiments
	;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ntends=6   ;last three term are 10-theta'*pian U'/pian x,11-U'*pian theta'/pian x,13 is sum of 11+12
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;number of theta tends:0mp(latent heat),1advection,2pbl(sensible heat)
	;3:localchange,4radiation,5:longwavw,6shortwavw,7diffusion
   p1000mb      = 100000.
   rcp=2/7 
   z   = wrf_user_getvar(a_ctl, "z",-1)         ; grid point height
   pi_phy=new((/ntimes,bottom_top,south_north,west_east/),float)
; First get the variables we will need   
theta=new((/ntends,nspecies,ntimes,bottom_top,south_north,west_east/),float)  
   theta_latlon=new((/ntends,nspecies,ntimes,bottom_top,sn,we/),float) 
   theta_areave=new((/ntends,nspecies,ntimes,bottom_top/),float) 
   theta_timeave=new((/ntends,nspecies,bottom_top/),float)    ;average over lat/lon
   theta=0
   theta_latlon=0
   theta_areave=0
   theta_timeave=0
pupx=new((/ntimes,bottom_top,south_north,west_east/),float)   
pvpy=new((/ntimes,bottom_top,south_north,west_east/),float) 
divv=new((/ntimes,bottom_top,south_north,west_east/),float)  ;pupx+pvpy
pmupt=new((/ntimes,south_north,west_east/),float) 
pmupt=0

    p   = wrf_user_getvar(a_ctl, "pressure",-1)         ; grid point height;hPa=mb
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
    ua   = wrf_user_getvar(a_ctl, "ua",-1)
    va   = wrf_user_getvar(a_ctl, "va",-1)
    th   = wrf_user_getvar(a_ctl, "th",-1)   ;theta,unit is  K
	dn  = a_ctl->DN(:,:)   ;;time,bottom_top
	mu=a_ctl->MU(:,:,:)   ;perturbation mu
    mub=a_ctl->MUB(:,:,:)
    mut=new((/ntimes,south_north,west_east/),float)
    mut=mu+mub
	
    muu=new((/ntimes,south_north,west_east/),float)
	muv=new((/ntimes,south_north,west_east/),float)
	msfuyw=a_ctl->MAPFAC_UY(:,:,:)  ; time,south_north,west_east
	msfvxw=a_ctl->MAPFAC_VX(:,:,:)  ;
	msftx=a_ctl->MAPFAC_MX(:,:,:)
    msfty=a_ctl->MAPFAC_MY(:,:,:)
	
	msfuy=new((/ntimes,south_north,west_east/),float);MAPFAC_UY(Time, south_north, west_east_stag) ;
	msfvx=new((/ntimes,south_north,west_east/),float); MAPFAC_VX(Time, south_north_stag, west_east) ;
	do it=0,west_east-1
	 msfuy(:,:,it)=(msfuyw(:,:,it+1)+msfuyw(:,:,it))/2
	end do
    do it=0,south_north-1
	msfvx(:,it,:)=(msfvxw(:,it+1,:)+msfvxw(:,it,:))/2
	end do
    muu(:,:,:) = mut(:,:,:)/msfuy(:,:,:)
    muv(:,:,:) = mut(:,:,:)/msfvx(:,:,:)
	rdx=1/3000
	rdy=1/3000
  do ix=0,west_east-2
   do it=0,bottom_top-1
	pupx(:,it,:,ix)=(muu(:,:,ix+1)*ua(:,it,:,ix+1)-muu(:,:,ix)*ua(:,it,:,ix))/3000   ;dx=3000 m
   end do
  end do
   pupx(:,:,:,west_east-1)=pupx(:,:,:,west_east-2)
  do iy=0,south_north-2
   do it=0,bottom_top-1
	pvpy(:,it,iy,:)=(muv(:,iy+1,:)*va(:,it,iy+1,:)-muv(:,iy,:)*va(:,it,iy,:))/3000   ;dx=3000 m
   end do 
 end do
   pvpy(:,:,south_north-1,:)=pvpy(:,:,south_north-2,:)
   
;Integrating with respect to nu over ALL levels, with dnu/dt=0 at top
  do it=0,ntimes-1
   do iz=0,bottom_top-1
   divv(it,iz,:,:)=(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz)
   pmupt(it,:,:)=pmupt(it,:,:)+(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz)
  end do
 end do
   pmupt=pmupt*msftx*msfty   ;coupled map factor  
   
  ;  theta(0,0,:,:,:,:)= a_ctl->H_DIABATIC(:,:,:,:)
	;theta(1,0,:,:,:,:)= a_ctl->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s
;	theta(2,0,:,:,:,:) = a_ctl->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s 	
;	theta(4,0,:,:,:,:)= a_ctl->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s
;	theta(5,0,:,:,:,:)= a_ctl->RTHRATLW(:,:,:,:)
;	theta(6,0,:,:,:,:) = a_ctl->RTHRATSW(:,:,:,:)
;	theta(7,0,:,:,:,:) = a_ctl->TDIF_TEND(:,:,:,:)
	theta(0,0,:,:,:,:) = a_ctl->ZTEND(:,:,:,:)   ; - p omegatheta p eta
 ww=new((/ntimes,bottom_top,south_north,west_east/),float)
  meanw=new((/ntimes,bottom_top/),float)
   meanth=new((/ntimes,bottom_top/),float)       ;mean theta
   pwbpeta=new((/ntimes,bottom_top/),float) ;pian(omega bar)/ pian eta
   wbthb=new((/ntimes,bottom_top/),float) 
   pwbthbpeta=new((/ntimes,bottom_top/),float)
   pwpthppeta=new((/ntimes,bottom_top,sn,we/),float)
   thww=new((/ntimes,bottom_top,south_north,west_east/),float)
   thw=new((/ntimes,bottom_top,sn,we/),float)
   pthbpeta=new((/ntimes,bottom_top/),float)
   thww(:,0,:,:)=0    ;vertical wind at the first layer and last layer is 0;
 do it=0,ntimes-1
  do iz=1,bottom_top-1
   thww(it,iz,:,:)=thww(it,iz-1,:,:)-theta(0,0,it,iz,:,:)*dn(it,iz)  ; - p omegatheta p eta
  end do
 end do
  ww(:,:,:,:)=thww(:,:,:,:)/th(:,:,:,:)
  thw(:,:,:,:)=thww(:,:,sn_s:sn_e,we_s:we_e)
  meanw(:,:)=dim_sum_n_Wrap(ww(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num
  meanth(:,:)=dim_sum_n_Wrap(th(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num   ;
 ;asciiwrite("meanw.dat",meanw(:,:))

; thw(:,:,:,:)=th(:,:,sn_s:sn_e,we_s:we_e)*ww(:,:,sn_s:sn_e,we_s:we_e)
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pthbpeta(it,iz+1)=(meanth(it,iz+1)-meanth(it,iz))/dn(it,iz+1)
  pwbpeta(it,iz+1)=(meanw(it,iz+1)-meanw(it,iz))/dn(it,iz+1)
 end do 
 end do
 pwbpeta(:,0)=pwbpeta(:,1)     ;pian omega bar /pian eta
 pthbpeta(:,0)=pthbpeta(:,1)
   wpie=new((/ntimes,bottom_top,sn,we/),float)  ;U'
   thpie=new((/ntimes,bottom_top,sn,we/),float)  ;theta'
   wpthp=new((/ntimes,bottom_top,sn,we/),float)  ;omega'*theta'
 do it=0,ntimes-1
  do iz=0,bottom_top-1
   wpie(it,iz,:,:)=ww(it,iz,sn_s:sn_e,we_s:we_e)-meanw(it,iz)  ;omega'=ww'=ww-meanww
   thpie(it,iz,:,:)=th(it,iz,sn_s:sn_e,we_s:we_e)-meanth(it,iz)
  end do
 end do 
 wpthp(:,:,:,:)=wpie(:,:,:,:)*thpie(:,:,:,:)
 wbthb=meanw*meanth
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pwbthbpeta(it,iz+1)=(wbthb(it,iz+1)-wbthb(it,iz))/dn(it,iz+1)
  pwpthppeta(it,iz+1,:,:)=(wpthp(it,iz+1,:,:)-wpthp(it,iz,:,:))/dn(it,iz+1)
 end do 
 end do
  pwbthbpeta(:,0)=pwbthbpeta(:,1)
  pwpthppeta(:,0,:,:)=pwpthppeta(:,1,:,:)
  
meanpwpthppeta=new((/ntimes,bottom_top/),float)
meanpwpthppeta=dim_sum_n_Wrap(pwpthppeta(:,:,:,:),(/2,3/))/grid_num
  
  do it=0,ntimes-1
  do iz=0,bottom_top-1
	theta(1,0,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)*pwbpeta(it,iz)    ;thetabar *pian omegabar/pian eta
	;the value for every lat/lon point is the same=constant averaged in lat/lon
     theta(2,0,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)
     theta(3,0,it,iz,sn_s:sn_e,we_s:we_e)=pwbpeta(it,iz);(/pthbpeta(it,iz)/)
	 theta(4,0,it,iz,sn_s:sn_e,we_s:we_e)=meanw(it,iz)*pthbpeta(it,iz)
	 theta(5,0,it,iz,sn_s:sn_e,we_s:we_e)=meanpwpthppeta(it,iz)  ;meanw(it,iz)
	; theta(6,0,it,iz,sn_s:sn_e,we_s:we_e)=pthbpeta(it,iz)
	end do
  end do
	
 ; do it=0,bottom_top-1
;	theta(8,0,:,it,:,:) = -pmupt(:,:,:)*th(:,it,:,:) ;;coupled by mu,;-p mu/p t is one of the right term in heat budget, 
;  end do

	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,0,:,ib,:,:)=theta(it,0,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
         do ib=0,bottom_top-1
 ;         theta(1,0,:,ib,:,:)=theta(1,0,:,ib,:,:)/mut(:,:,:)
 ;         theta(2,0,:,ib,:,:)=theta(2,0,:,ib,:,:)/mut(:,:,:)
 ;         theta(4,0,:,ib,:,:)=theta(4,0,:,ib,:,:)/mut(:,:,:)
;		  theta(8,0,:,ib,:,:)=theta(8,0,:,ib,:,:)/mut(:,:,:)
		  theta(0,0,:,ib,:,:)=theta(0,0,:,ib,:,:)/mut(:,:,:)
		  theta(1,0,:,ib,:,:)=theta(1,0,:,ib,:,:)/mut(:,:,:)
		 ; theta(11,0,:,ib,:,:)=theta(11,0,:,ib,:,:)/mut(:,:,:)
		 theta(3,0,:,ib,:,:)=theta(3,0,:,ib,:,:)/mut(:,:,:)
		 theta(4,0,:,ib,:,:)=theta(4,0,:,ib,:,:)/mut(:,:,:)
		 theta(5,0,:,ib,:,:)=theta(5,0,:,ib,:,:)/mut(:,:,:)
		 ;theta(15,0,:,ib,:,:)=theta(15,0,:,ib,:,:)/mut(:,:,:)
         end do
	delete(mut)
delete(pupx)
delete(pvpy)
delete(divv) ;pupx+pvpy 
delete(pmupt)
delete(msfuy)
delete(msfvx)
delete(muu)
delete(muv)
delete(meanth)
delete(meanw)
delete(pwbpeta)
delete(ww)
delete(wbthb)
delete(pwbthbpeta)
delete(pwpthppeta)
delete(thww)
delete(thw)
delete(pthbpeta)
delete(wpie)
delete(thpie)
delete(wpthp)
delete(meanpwpthppeta)
	
pupx=new((/ntimes,bottom_top,south_north,west_east/),float)   
pvpy=new((/ntimes,bottom_top,south_north,west_east/),float) 
divv=new((/ntimes,bottom_top,south_north,west_east/),float)   ;pupx+pvpy
pmupt=new((/ntimes,south_north,west_east/),float) 
pmupt=0

    p   = wrf_user_getvar(a_NCR, "pressure",-1)         ; grid point height;hPa=mb
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
    ua   = wrf_user_getvar(a_NCR, "ua",-1)
    va   = wrf_user_getvar(a_NCR, "va",-1)
    th   = wrf_user_getvar(a_NCR, "th",-1)   ;theta,unit is  K
	dn  = a_NCR->DN(:,:)   ;;time,bottom_top
	mu=a_NCR->MU(:,:,:)   ;perturbation mu
    mub=a_NCR->MUB(:,:,:)
    mut=new((/ntimes,south_north,west_east/),float)
    mut=mu+mub
	
    muu=new((/ntimes,south_north,west_east/),float)
	muv=new((/ntimes,south_north,west_east/),float)
	msfuyw=a_NCR->MAPFAC_UY(:,:,:)  ; 'time,south_north,west_east
	msfvxw=a_NCR->MAPFAC_VX(:,:,:)  ;
	msftx=a_NCR->MAPFAC_MX(:,:,:)
    msfty=a_NCR->MAPFAC_MY(:,:,:)
	msfuy=new((/ntimes,south_north,west_east/),float);MAPFAC_UY(Time, south_north, west_east_stag) ;
	msfvx=new((/ntimes,south_north,west_east/),float); MAPFAC_VX(Time, south_north_stag, west_east) ;
	do it=0,west_east-1
	 msfuy(:,:,it)=(msfuyw(:,:,it+1)+msfuyw(:,:,it))/2
	end do
    do it=0,south_north-1
	msfvx(:,it,:)=(msfvxw(:,it+1,:)+msfvxw(:,it,:))/2
	end do
    muu = mut/msfuy
    muv = mut/msfvx
	rdx=1/3000
	rdy=1/3000
  do ix=0,west_east-2
   do it=0,bottom_top-1
	pupx(:,it,:,ix)=(muu(:,:,ix+1)*ua(:,it,:,ix+1)-muu(:,:,ix)*ua(:,it,:,ix))/3000   ;dx=3000 m
   end do
  end do
   pupx(:,:,:,west_east-1)=pupx(:,:,:,west_east-2)
  do iy=0,south_north-2
   do it=0,bottom_top-1
	pvpy(:,it,iy,:)=(muv(:,iy+1,:)*va(:,it,iy+1,:)-muv(:,iy,:)*va(:,it,iy,:))/3000   ;dx=3000 m
   end do 
 end do
   pvpy(:,:,south_north-1,:)=pvpy(:,:,south_north-2,:)
;Integrating with respect to nu over ALL levels, with dnu/dt=0 at top

 do it=0,ntimes-1
  do iz=0,bottom_top-1
        divv(it,iz,:,:)=(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz)
   pmupt(it,:,:)=pmupt(it,:,:)+(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz) 
  end do
 end do
   pmupt=pmupt*msftx*msfty   ;coupled map factor 
  ;  theta(0,1,:,:,:,:)= a_NCR->H_DIABATIC(:,:,:,:)
;	theta(1,1,:,:,:,:)= a_NCR->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s
;	theta(2,1,:,:,:,:) = a_NCR->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s 	
;	theta(4,1,:,:,:,:)= a_NCR->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s
;	theta(5,1,:,:,:,:)= a_NCR->RTHRATLW(:,:,:,:)
;	theta(6,1,:,:,:,:) = a_NCR->RTHRATSW(:,:,:,:)
;	theta(7,1,:,:,:,:) = a_NCR->TDIF_TEND(:,:,:,:)
	theta(0,1,:,:,:,:) = a_NCR->ZTEND(:,:,:,:)
 ww=new((/ntimes,bottom_top,south_north,west_east/),float)
  meanw=new((/ntimes,bottom_top/),float)
   meanth=new((/ntimes,bottom_top/),float)       ;mean theta
   pwbpeta=new((/ntimes,bottom_top/),float) ;pian(omega bar)/ pian eta
   wbthb=new((/ntimes,bottom_top/),float) 
   pwbthbpeta=new((/ntimes,bottom_top/),float)
   pwpthppeta=new((/ntimes,bottom_top,sn,we/),float)
   thww=new((/ntimes,bottom_top,south_north,west_east/),float)
   thw=new((/ntimes,bottom_top,sn,we/),float)
   pthbpeta=new((/ntimes,bottom_top/),float)
   thww(:,0,:,:)=0    ;vertical wind at the first layer and last layer is 0;
 do it=0,ntimes-1
  do iz=1,bottom_top-1
   thww(it,iz,:,:)=thww(it,iz-1,:,:)-theta(0,1,it,iz,:,:)*dn(it,iz)   ; - p omegatheta p eta
  end do
 end do
  ww(:,:,:,:)=thww(:,:,:,:)/th(:,:,:,:)
  thw(:,:,:,:)=thww(:,:,sn_s:sn_e,we_s:we_e)
    meanw(:,:)=dim_sum_n_Wrap(ww(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num
  meanth(:,:)=dim_sum_n_Wrap(th(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num   ;
 ;thw(:,:,:,:)=th(:,:,sn_s:sn_e,we_s:we_e)*ww(:,:,sn_s:sn_e,we_s:we_e)
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pwbpeta(it,iz+1)=(meanw(it,iz+1)-meanw(it,iz))/dn(it,iz+1)
  pthbpeta(it,iz+1)=(meanth(it,iz+1)-meanth(it,iz))/dn(it,iz+1)
 end do 
 end do
 pwbpeta(:,0)=pwbpeta(:,1)     ;pian omega bar /pian eta
 pthbpeta(:,0)=pthbpeta(:,1)
   wpie=new((/ntimes,bottom_top,sn,we/),float)  ;U'
   thpie=new((/ntimes,bottom_top,sn,we/),float)  ;theta'
   wpthp=new((/ntimes,bottom_top,sn,we/),float)  ;omega'*theta'
 do it=0,ntimes-1
  do iz=0,bottom_top-1
   wpie(it,iz,:,:)=ww(it,iz,sn_s:sn_e,we_s:we_e)-meanw(it,iz)  ;omega'=ww'=ww-meanww
   thpie(it,iz,:,:)=th(it,iz,sn_s:sn_e,we_s:we_e)-meanth(it,iz)
  end do
 end do 
 wpthp(:,:,:,:)=wpie(:,:,:,:)*thpie(:,:,:,:)
 wbthb=meanw*meanth
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pwbthbpeta(it,iz+1)=(wbthb(it,iz+1)-wbthb(it,iz))/dn(it,iz+1)
  pwpthppeta(it,iz+1,:,:)=(wpthp(it,iz+1,:,:)-wpthp(it,iz,:,:))/dn(it,iz+1)
 end do 
 end do
  pwpthppeta(:,0,:,:)=pwpthppeta(:,1,:,:)
  pwbthbpeta(:,0)=pwbthbpeta(:,1)
meanpwpthppeta=new((/ntimes,bottom_top/),float)
meanpwpthppeta=dim_sum_n_Wrap(pwpthppeta(:,:,:,:),(/2,3/))/grid_num
  do it=0,ntimes-1
  do iz=0,bottom_top-1
	theta(1,1,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)*pwbpeta(it,iz)    ;thetabar *pian omegabar/pian eta
	;the value for every lat/lon point is the same=constant averaged in lat/lon
	theta(2,1,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)
    theta(3,1,it,iz,sn_s:sn_e,we_s:we_e)=pwbpeta(it,iz);(/pthbpeta(it,iz)/)
	theta(4,1,it,iz,sn_s:sn_e,we_s:we_e)=meanw(it,iz)*pthbpeta(it,iz)
	theta(5,1,it,iz,sn_s:sn_e,we_s:we_e)=meanpwpthppeta(it,iz)  ;meanw(it,iz)
	;theta(15,1,it,iz,sn_s:sn_e,we_s:we_e)=pthbpeta(it,iz)
  end do
  end do
 ; do it=0,bottom_top-1
;	theta(8,1,:,it,:,:) = -pmupt(:,:,:)*th(:,it,:,:);;coupled by mu,;-p mu/p t is one of the right term in heat budget, 
 ; end do
	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,1,:,ib,:,:)=theta(it,1,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
         do ib=0,bottom_top-1
 ;         theta(1,1,:,ib,:,:)=theta(1,1,:,ib,:,:)/mut(:,:,:)
 ;         theta(2,1,:,ib,:,:)=theta(2,1,:,ib,:,:)/mut(:,:,:)
 ;         theta(4,1,:,ib,:,:)=theta(4,1,:,ib,:,:)/mut(:,:,:)
;		  theta(8,1,:,ib,:,:)=theta(8,1,:,ib,:,:)/mut(:,:,:)
		  theta(0,1,:,ib,:,:)=theta(0,1,:,ib,:,:)/mut(:,:,:)
		  theta(1,1,:,ib,:,:)=theta(1,1,:,ib,:,:)/mut(:,:,:)
		 ; theta(11,1,:,ib,:,:)=theta(11,1,:,ib,:,:)/mut(:,:,:)
		 theta(3,1,:,ib,:,:)=theta(3,1,:,ib,:,:)/mut(:,:,:)
		 theta(4,1,:,ib,:,:)=theta(4,1,:,ib,:,:)/mut(:,:,:)
		 theta(5,1,:,ib,:,:)=theta(5,1,:,ib,:,:)/mut(:,:,:)
		 ;theta(15,1,:,ib,:,:)=theta(15,1,:,ib,:,:)/mut(:,:,:)
         end do
	delete(mut)
delete(pupx)
delete(pvpy)
delete(divv) ;pupx+pvpy 
delete(pmupt)
delete(msfuy)
delete(msfvx)
delete(muu)
delete(muv)
delete(meanth)
delete(meanw)
delete(pwbpeta)
delete(ww)
delete(wbthb)
delete(pwbthbpeta)
delete(pwpthppeta)
delete(thww)
delete(thw)
delete(pthbpeta)
delete(wpie)
delete(thpie)
delete(wpthp)
delete(meanpwpthppeta)
	
pupx=new((/ntimes,bottom_top,south_north,west_east/),float)   
pvpy=new((/ntimes,bottom_top,south_north,west_east/),float) 
divv=new((/ntimes,bottom_top,south_north,west_east/),float)   ;pupx+pvpy
pmupt=new((/ntimes,south_north,west_east/),float) 
pmupt=0
;ww=new((/ntimes,bottom_top,south_north,west_east/),float)
    p   = wrf_user_getvar(a_NIR, "pressure",-1)         ; grid point height;hPa=mb
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
    ua   = wrf_user_getvar(a_NIR, "ua",-1)
    va   = wrf_user_getvar(a_NIR, "va",-1)
    th   = wrf_user_getvar(a_NIR, "th",-1)   ;theta,unit is  K
	dn  = a_NIR->DN(:,:)   ;;time,bottom_top
	mu=a_NIR->MU(:,:,:)   ;perturbation mu
    mub=a_NIR->MUB(:,:,:)
    mut=new((/ntimes,south_north,west_east/),float)
    mut=mu+mub
	
    muu=new((/ntimes,south_north,west_east/),float)
	muv=new((/ntimes,south_north,west_east/),float)
	msfuyw=a_NIR->MAPFAC_UY(:,:,:)  ; 'time,south_north,west_east
	msfvxw=a_NIR->MAPFAC_VX(:,:,:)  ;
	msftx=a_NIR->MAPFAC_MX(:,:,:)
    msfty=a_NIR->MAPFAC_MY(:,:,:)
   msfuy=new((/ntimes,south_north,west_east/),float);MAPFAC_UY(Time, south_north, west_east_stag) ;
	msfvx=new((/ntimes,south_north,west_east/),float); MAPFAC_VX(Time, south_north_stag, west_east) ;
	do it=0,west_east-1
	 msfuy(:,:,it)=(msfuyw(:,:,it+1)+msfuyw(:,:,it))/2
	end do
    do it=0,south_north-1
	msfvx(:,it,:)=(msfvxw(:,it+1,:)+msfvxw(:,it,:))/2
	end do
    muu = mut/msfuy
    muv = mut/msfvx
	rdx=1/3000
	rdy=1/3000
  do ix=0,west_east-2
   do it=0,bottom_top-1
	pupx(:,it,:,ix)=(muu(:,:,ix+1)*ua(:,it,:,ix+1)-muu(:,:,ix)*ua(:,it,:,ix))/3000   ;dx=3000 m
   end do
  end do
   pupx(:,:,:,west_east-1)=pupx(:,:,:,west_east-2)
  do iy=0,south_north-2
   do it=0,bottom_top-1
	pvpy(:,it,iy,:)=(muv(:,iy+1,:)*va(:,it,iy+1,:)-muv(:,iy,:)*va(:,it,iy,:))/3000   ;dx=3000 m
   end do 
 end do
   pvpy(:,:,south_north-1,:)=pvpy(:,:,south_north-2,:)
;Integrating with respect to nu over ALL levels, with dnu/dt=0 at top

  do it=0,ntimes-1
    do iz=0,bottom_top-1
	      divv(it,iz,:,:)=(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz)
   pmupt(it,:,:)=pmupt(it,:,:)+(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz)
  end do
 end do
   pmupt=pmupt*msftx*msfty   ;coupled map factor 
 ;   theta(0,2,:,:,:,:)= a_NIR->H_DIABATIC(:,:,:,:)
;	theta(1,2,:,:,:,:)= a_NIR->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s
;	theta(2,2,:,:,:,:) = a_NIR->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s 	
;	theta(4,2,:,:,:,:)= a_NIR->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s
;	theta(5,2,:,:,:,:)= a_NIR->RTHRATLW(:,:,:,:)
;	theta(6,2,:,:,:,:) = a_NIR->RTHRATSW(:,:,:,:)
;	theta(7,2,:,:,:,:) = a_NIR->TDIF_TEND(:,:,:,:)
	theta(0,2,:,:,:,:) = a_NIR->ZTEND(:,:,:,:)
 ww=new((/ntimes,bottom_top,south_north,west_east/),float)
  meanw=new((/ntimes,bottom_top/),float)
   meanth=new((/ntimes,bottom_top/),float)       ;mean theta
   pwbpeta=new((/ntimes,bottom_top/),float) ;pian(omega bar)/ pian eta
      wbthb=new((/ntimes,bottom_top/),float) 
   pwbthbpeta=new((/ntimes,bottom_top/),float)
   pwpthppeta=new((/ntimes,bottom_top,sn,we/),float)
   thww=new((/ntimes,bottom_top,south_north,west_east/),float)
   thw=new((/ntimes,bottom_top,sn,we/),float)
   pthbpeta=new((/ntimes,bottom_top/),float)
   thww(:,0,:,:)=0    ;vertical wind at the first layer and last layer is 0;
 do it=0,ntimes-1
  do iz=1,bottom_top-1
   thww(it,iz,:,:)=thww(it,iz-1,:,:)-theta(0,2,it,iz,:,:)*dn(it,iz)   ; - p omegatheta p eta
  end do
 end do
  ww(:,:,:,:)=thww(:,:,:,:)/th(:,:,:,:)
  thw(:,:,:,:)=thww(:,:,sn_s:sn_e,we_s:we_e)
    meanw(:,:)=dim_sum_n_Wrap(ww(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num
  meanth(:,:)=dim_sum_n_Wrap(th(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num   ;
 ;thw(:,:,:,:)=th(:,:,sn_s:sn_e,we_s:we_e)*ww(:,:,sn_s:sn_e,we_s:we_e)
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pwbpeta(it,iz+1)=(meanw(it,iz+1)-meanw(it,iz))/dn(it,iz+1)
  pthbpeta(it,iz+1)=(meanth(it,iz+1)-meanth(it,iz))/dn(it,iz+1)
 end do 
 end do
 pwbpeta(:,0)=pwbpeta(:,1)     ;pian omega bar /pian eta
 pthbpeta(:,0)=pthbpeta(:,1)
   wpie=new((/ntimes,bottom_top,sn,we/),float)  ;U'
   thpie=new((/ntimes,bottom_top,sn,we/),float)  ;theta'
   wpthp=new((/ntimes,bottom_top,sn,we/),float)  ;omega'*theta'
 do it=0,ntimes-1
  do iz=0,bottom_top-1
   wpie(it,iz,:,:)=ww(it,iz,sn_s:sn_e,we_s:we_e)-meanw(it,iz)  ;omega'=ww'=ww-meanww
   thpie(it,iz,:,:)=th(it,iz,sn_s:sn_e,we_s:we_e)-meanth(it,iz)
  end do
 end do 
 wpthp(:,:,:,:)=wpie(:,:,:,:)*thpie(:,:,:,:)
 wbthb=meanw*meanth
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pwbthbpeta(it,iz+1)=(wbthb(it,iz+1)-wbthb(it,iz))/dn(it,iz+1)
  pwpthppeta(it,iz+1,:,:)=(wpthp(it,iz+1,:,:)-wpthp(it,iz,:,:))/dn(it,iz+1)
 end do 
 end do
  pwpthppeta(:,0,:,:)=pwpthppeta(:,1,:,:)
  pwbthbpeta(:,0)=pwbthbpeta(:,1)
  meanpwpthppeta=new((/ntimes,bottom_top/),float)
  meanpwpthppeta=dim_sum_n_Wrap(pwpthppeta(:,:,:,:),(/2,3/))/grid_num
  do it=0,ntimes-1
  do iz=0,bottom_top-1
	theta(1,2,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)*pwbpeta(it,iz)    ;thetabar *pian omegabar/pian eta
	;the value for every lat/lon point is the same=constant averaged in lat/lon
	theta(2,2,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)
    theta(3,2,it,iz,sn_s:sn_e,we_s:we_e)=pwbpeta(it,iz);(/pthbpeta(it,iz)/)
	theta(4,2,it,iz,sn_s:sn_e,we_s:we_e)=meanw(it,iz)*pthbpeta(it,iz)
	theta(5,2,it,iz,sn_s:sn_e,we_s:we_e)=meanpwpthppeta(it,iz)  ;meanw(it,iz)
	;theta(15,2,it,iz,sn_s:sn_e,we_s:we_e)=pthbpeta(it,iz)
  end do
  end do
 ; do it=0,bottom_top-1
;	theta(8,2,:,it,:,:) = -pmupt(:,:,:)*th(:,it,:,:);;coupled by mu,;-p mu/p t is one of the right term in heat budget, 
 ; end do
	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,2,:,ib,:,:)=theta(it,2,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
         do ib=0,bottom_top-1
  ;        theta(1,2,:,ib,:,:)=theta(1,2,:,ib,:,:)/mut(:,:,:)
 ;         theta(2,2,:,ib,:,:)=theta(2,2,:,ib,:,:)/mut(:,:,:)
 ;         theta(4,2,:,ib,:,:)=theta(4,2,:,ib,:,:)/mut(:,:,:)
;		  theta(8,2,:,ib,:,:)=theta(8,2,:,ib,:,:)/mut(:,:,:)
		  theta(0,2,:,ib,:,:)=theta(0,2,:,ib,:,:)/mut(:,:,:)
		  theta(1,2,:,ib,:,:)=theta(1,2,:,ib,:,:)/mut(:,:,:)
		 ; theta(11,2,:,ib,:,:)=theta(11,2,:,ib,:,:)/mut(:,:,:)
		 theta(3,2,:,ib,:,:)=theta(3,2,:,ib,:,:)/mut(:,:,:)
		 theta(4,2,:,ib,:,:)=theta(4,2,:,ib,:,:)/mut(:,:,:)
		 theta(5,2,:,ib,:,:)=theta(5,2,:,ib,:,:)/mut(:,:,:)
		 ;theta(6,2,:,ib,:,:)=theta(6,2,:,ib,:,:)/mut(:,:,:)
         end do
   	delete(mut)
delete(pupx)
delete(pvpy)
delete(divv) ;pupx+pvpy 
delete(pmupt)
delete(msfuy)
delete(msfvx)
delete(muu)
delete(muv)
delete(meanth)
delete(meanw)
delete(pwbpeta)
delete(ww)
delete(wbthb)
delete(pwbthbpeta)
delete(pwpthppeta)
delete(thww)
delete(thw)
delete(pthbpeta)
delete(wpie)
delete(thpie)
delete(wpthp)
delete(meanpwpthppeta)
	
pupx=new((/ntimes,bottom_top,south_north,west_east/),float)   
pvpy=new((/ntimes,bottom_top,south_north,west_east/),float)
divv=new((/ntimes,bottom_top,south_north,west_east/),float)   ;pupx+pvpy 
pmupt=new((/ntimes,south_north,west_east/),float) 
pmupt=0
;ww=new((/ntimes,bottom_top,south_north,west_east/),float)
    p   = wrf_user_getvar(a_NSR, "pressure",-1)         ; grid point height;hPa=mb
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
    ua   = wrf_user_getvar(a_NSR, "ua",-1)
    va   = wrf_user_getvar(a_NSR, "va",-1)
    th   = wrf_user_getvar(a_NSR, "th",-1)   ;theta,unit is  K
	dn  = a_NSR->DN(:,:)   ;;time,bottom_top
	mu=a_NSR->MU(:,:,:)   ;perturbation mu
    mub=a_NSR->MUB(:,:,:)
    mut=new((/ntimes,south_north,west_east/),float)
    mut=mu+mub
	
    muu=new((/ntimes,south_north,west_east/),float)
	muv=new((/ntimes,south_north,west_east/),float)
	msfuyw=a_NSR->MAPFAC_UY(:,:,:)  ; 'time,south_north,west_east
	msfvxw=a_NSR->MAPFAC_VX(:,:,:)  ;
	msftx=a_NSR->MAPFAC_MX(:,:,:)
    msfty=a_NSR->MAPFAC_MY(:,:,:)
	msfuy=new((/ntimes,south_north,west_east/),float);MAPFAC_UY(Time, south_north, west_east_stag) ;
	msfvx=new((/ntimes,south_north,west_east/),float); MAPFAC_VX(Time, south_north_stag, west_east) ;
	do it=0,west_east-1
	 msfuy(:,:,it)=(msfuyw(:,:,it+1)+msfuyw(:,:,it))/2
	end do
    do it=0,south_north-1
	msfvx(:,it,:)=(msfvxw(:,it+1,:)+msfvxw(:,it,:))/2
	end do
    muu = mut/msfuy
    muv = mut/msfvx
	rdx=1/3000
	rdy=1/3000
  do ix=0,west_east-2
   do it=0,bottom_top-1
	pupx(:,it,:,ix)=(muu(:,:,ix+1)*ua(:,it,:,ix+1)-muu(:,:,ix)*ua(:,it,:,ix))/3000   ;dx=3000 m
   end do
  end do
   pupx(:,:,:,west_east-1)=pupx(:,:,:,west_east-2)
  do iy=0,south_north-2
   do it=0,bottom_top-1
	pvpy(:,it,iy,:)=(muv(:,iy+1,:)*va(:,it,iy+1,:)-muv(:,iy,:)*va(:,it,iy,:))/3000   ;dx=3000 m
   end do 
 end do
   pvpy(:,:,south_north-1,:)=pvpy(:,:,south_north-2,:)
       asciiwrite("pvpy.dat",pvpy(:,:,:,:))
;Integrating with respect to nu over ALL levels, with dnu/dt=0 at top

 do it=0,ntimes-1
  do iz=0,bottom_top-1
  divv(it,iz,:,:)=(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz)
   pmupt(it,:,:)=pmupt(it,:,:)+(pupx(it,iz,:,:)+pvpy(it,iz,:,:))*dn(it,iz)
  end do
 end do
   pmupt=pmupt*msftx*msfty   ;coupled map factor 
 ;   theta(0,3,:,:,:,:)= a_NSR->H_DIABATIC(:,:,:,:)
;	theta(1,3,:,:,:,:)= a_NSR->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s
;;	theta(2,3,:,:,:,:) = a_NSR->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s 	
;	theta(4,3,:,:,:,:)= a_NSR->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s
;	theta(5,3,:,:,:,:)= a_NSR->RTHRATLW(:,:,:,:)
;	theta(6,3,:,:,:,:) = a_NSR->RTHRATSW(:,:,:,:)
;	theta(7,3,:,:,:,:) = a_NSR->TDIF_TEND(:,:,:,:)
	theta(0,3,:,:,:,:) = a_NSR->ZTEND(:,:,:,:)   ; - p omegatheta p eta
 ww=new((/ntimes,bottom_top,south_north,west_east/),float)
  meanw=new((/ntimes,bottom_top/),float)
   meanth=new((/ntimes,bottom_top/),float)       ;mean theta
      wbthb=new((/ntimes,bottom_top/),float) 
   pwbthbpeta=new((/ntimes,bottom_top/),float)
   pwbpeta=new((/ntimes,bottom_top/),float) ;pian(omega bar)/ pian eta
   pwpthppeta=new((/ntimes,bottom_top,sn,we/),float)
   thww=new((/ntimes,bottom_top,south_north,west_east/),float)
   thw=new((/ntimes,bottom_top,sn,we/),float)
   pthbpeta=new((/ntimes,bottom_top/),float)
   thww(:,0,:,:)=0    ;vertical wind at the first layer and last layer is 0;
 do it=0,ntimes-1
  do iz=1,bottom_top-1
   thww(it,iz,:,:)=thww(it,iz-1,:,:)-theta(0,3,it,iz,:,:)*dn(it,iz)   ; - p omegatheta p eta
  end do
 end do
  ww(:,:,:,:)=thww(:,:,:,:)/th(:,:,:,:)
  thw(:,:,:,:)=thww(:,:,sn_s:sn_e,we_s:we_e)
    meanw(:,:)=dim_sum_n_Wrap(ww(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num
  meanth(:,:)=dim_sum_n_Wrap(th(:,:,sn_s:sn_e,we_s:we_e),(/2,3/))/grid_num   ;
 ;thw(:,:,:,:)=th(:,:,sn_s:sn_e,we_s:we_e)*ww(:,:,sn_s:sn_e,we_s:we_e)
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pwbpeta(it,iz+1)=(meanw(it,iz+1)-meanw(it,iz))/dn(it,iz+1)
  pthbpeta(it,iz+1)=(meanth(it,iz+1)-meanth(it,iz))/dn(it,iz+1)
 end do 
 end do
 pwbpeta(:,0)=pwbpeta(:,1)     ;pian omega bar /pian eta
 pthbpeta(:,0)=pthbpeta(:,1)
   wpie=new((/ntimes,bottom_top,sn,we/),float)  ;U'
   thpie=new((/ntimes,bottom_top,sn,we/),float)  ;theta'
   wpthp=new((/ntimes,bottom_top,sn,we/),float)  ;omega'*theta'
 do it=0,ntimes-1
  do iz=0,bottom_top-1
   wpie(it,iz,:,:)=ww(it,iz,sn_s:sn_e,we_s:we_e)-meanw(it,iz)  ;omega'=ww'=ww-meanww
   thpie(it,iz,:,:)=th(it,iz,sn_s:sn_e,we_s:we_e)-meanth(it,iz)
  end do
 end do 
 wpthp(:,:,:,:)=wpie(:,:,:,:)*thpie(:,:,:,:)
 wbthb=meanw*meanth
 do it=0,ntimes-1
 do iz=0,bottom_top-2
  pwbthbpeta(it,iz+1)=(wbthb(it,iz+1)-wbthb(it,iz))/dn(it,iz+1)
  pwpthppeta(it,iz+1,:,:)=(wpthp(it,iz+1,:,:)-wpthp(it,iz,:,:))/dn(it,iz+1)
 end do 
 end do
  pwpthppeta(:,0,:,:)=pwpthppeta(:,1,:,:)
  pwbthbpeta(:,0)=pwbthbpeta(:,1)
meanpwpthppeta=new((/ntimes,bottom_top/),float)
meanpwpthppeta=dim_sum_n_Wrap(pwpthppeta(:,:,:,:),(/2,3/))/grid_num
  do it=0,ntimes-1
  do iz=0,bottom_top-1
	theta(1,3,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)*pwbpeta(it,iz)    ;thetabar *pian omegabar/pian eta
	;the value for every lat/lon point is the same=constant averaged in lat/lon
    theta(2,3,it,iz,sn_s:sn_e,we_s:we_e)=meanth(it,iz)
    theta(3,3,it,iz,sn_s:sn_e,we_s:we_e)=pwbpeta(it,iz);(/pthbpeta(it,iz)/)
	theta(4,3,it,iz,sn_s:sn_e,we_s:we_e)=meanw(it,iz)*pthbpeta(it,iz)
	theta(5,3,it,iz,sn_s:sn_e,we_s:we_e)=meanpwpthppeta(it,iz)  ;meanw(it,iz)
	;theta(15,3,it,iz,sn_s:sn_e,we_s:we_e)=pthbpeta(it,iz)
  end do
  end do
 ; do it=0,bottom_top-1
;	theta(8,3,:,it,:,:) = -pmupt(:,:,:)*th(:,it,:,:);;coupled by mu,;-p mu/p t is one of the right term in heat budget, 
 ; end do

	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,3,:,ib,:,:)=theta(it,3,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
         do ib=0,bottom_top-1
  ;        theta(1,3,:,ib,:,:)=theta(1,3,:,ib,:,:)/mut(:,:,:)
 ;         theta(2,3,:,ib,:,:)=theta(2,3,:,ib,:,:)/mut(:,:,:)
 ;         theta(4,3,:,ib,:,:)=theta(4,3,:,ib,:,:)/mut(:,:,:)
;		  theta(8,3,:,ib,:,:)=theta(8,3,:,ib,:,:)/mut(:,:,:)
		  theta(0,3,:,ib,:,:)=theta(0,3,:,ib,:,:)/mut(:,:,:)
		  theta(1,3,:,ib,:,:)=theta(1,3,:,ib,:,:)/mut(:,:,:)
		 ; theta(11,3,:,ib,:,:)=theta(11,3,:,ib,:,:)/mut(:,:,:)
		 theta(3,3,:,ib,:,:)=theta(3,3,:,ib,:,:)/mut(:,:,:)
		 theta(4,3,:,ib,:,:)=theta(4,3,:,ib,:,:)/mut(:,:,:)
		 theta(5,3,:,ib,:,:)=theta(5,3,:,ib,:,:)/mut(:,:,:)
		 ;theta(15,3,:,ib,:,:)=theta(15,3,:,ib,:,:)/mut(:,:,:)
         end do

;theta(3,:,:,:,:,:)= theta(0,:,:,:,:,:)+theta(1,:,:,:,:,:)+theta(2,:,:,:,:,:)+theta(4,:,:,:,:,:)+theta(8,:,:,:,:,:)  ;local change
	
 ; theta_latlon=new((/ntends,nspecies,ntimes,bottom_top,sn,we/),float) 
  theta_latlon(:,:,:,:,:,:)=theta(:,:,:,:,sn_s:sn_e,we_s:we_e)
  
 ;  theta_areave=new((/ntends,nspecies,ntimes,bottom_top/),float) 
   theta_areave=dim_sum_n_Wrap(theta_latlon(:,:,:,:,:,:),(/4,5/))  ;sum over lat/lon
   theta_areave = theta_areave/grid_num     ;per grid
;asciiwrite("test0.dat",theta_areave)   
 ; theta_timeave=new((/ntends,nspecies,bottom_top/),float)    ;average over lat/lon
  starttime=13
  endtime=36
  timespan=endtime-starttime+1
  do nt=0,ntends-1
   do i=0,nspecies-1
    do j=0,bottom_top-1
     theta_timeave(nt,i,j)=sum(theta_areave(nt,i,starttime:endtime,j))/timespan  ;average theta from startime to endtime
    end do
   end do
  end do


;ntends=10:0-mp(latent heat),1-advection,2-pbl(sensible heat),3-:localchange,4-radiation,5-longwavw,6-shortwavw,7-diffussion,
;;8-pmu/pt,9-xtend  ;total 10 terms
;tc_bias=new((/ntends,nspecies+1,bottom_top/),float) 
;tc_bias(:,0,:)=theta_timeave(:,3,:)-theta_timeave(:,0,:)   ;NSR-CTL
;tc_bias(:,1,:)=theta_timeave(:,1,:)-theta_timeave(:,2,:)   ;NCR-NIR
;tc_bias(:,2,:)=theta_timeave(:,2,:)-theta_timeave(:,0,:)   ;NIR-CTL
;tc_bias(:,3,:)=theta_timeave(:,1,:)-theta_timeave(:,3,:)   ;NCR-NSR
;tc_bias(:,4,:)=0     ;line:x=0
;tc_bias=tc_bias*3600    ;convert from celsius per second to celsius per hour
;;9-"ztend",10-"thb*p omegab/p eta",11-"thetabar",12-"p omegab/p eta",13-"wb*p thb/p eta",14-"wbar",15-"p thb/p eta"
pt_bias=new((/6,nspecies-1,bottom_top/),float);7 kinds  ,ntends=16,here present 9,10,11,12,last term 0 line  
                                            ;nspecies=4,sensitivity experiments
do i=0,5  ;0,3 corresponds 9-15
 pt_bias(i,0,:)=theta_timeave(i,1,:)-theta_timeave(i,0,:)   ;NSR-CTL
 pt_bias(i,1,:)=theta_timeave(i,2,:)-theta_timeave(i,0,:)   ;NCR-NIR
 pt_bias(i,2,:)=theta_timeave(i,3,:)-theta_timeave(i,0,:)   ;NIR-CTL
; pt_bias(i,3,:)=theta_timeave(i+9,1,:)-theta_timeave(i+9,3,:)   ;NCR-NSR
; pt_bias(i,4,:)=theta_timeave(i+9,1,:)-theta_timeave(i+9,0,:)   ;NCR-CTL 
end do
 ;pt_bias(6,:,:)=0  
;pt_bias(1,:,:)=-pt_bias(1,:,:)
;pt_bias(2,:,:)=pt_bias(2,:,:)/10000   ;11-thetabar the magtitude of meantheta is about 2-3 larger than other terms;
;pt_bias(6,:,:)=pt_bias(6,:,:)/1000/100 ;15-"p thb/p eta"
;pt_bias(3,:,:)=-100*pt_bias(3,:,:)  ;12-"p omegab/p eta"
;pt_bias(4,:,:)=-10*pt_bias(4,:,:)  ;13-"wb*p thb/p eta"
;pt_bias(5,:,:)=-100*pt_bias(5,:,:)  ;14-"wbar"
pt_bias(0,:,:)=pt_bias(0,:,:)*3600*24 ;ztend ;convert from celsius per second to celsius per day
pt_bias(1,:,:)=-pt_bias(1,:,:)*3600*24    ;"thb*p omegab/p eta",do not need convert from celsius per second to celsius per day
pt_bias(2,:,:)=10*pt_bias(2,:,:)    ;"thetabar" do not need convert from celsius per second to celsius per day
pt_bias(3,:,:)=-100*pt_bias(3,:,:)*3600*24    ;"p omegab/p eta" do not need convert from celsius per second to celsius per day
pt_bias(4,:,:)=-pt_bias(4,:,:)*3600*24    ;"wb*p thb/p eta" do not need convert from celsius per second to celsius per day
pt_bias(5,:,:)=-pt_bias(5,:,:)*3600*24    ;"p omega'th' p eta" convert from celsius per second to celsius per day
;pt_bias(6,:,:)=pt_bias    ;do not need convert from celsius per second to celsius per day

pt_bias=-pt_bias ;;;reverse to 0:CTL-NSR,1:NIR-NCR,2:CTL-NIR,3:NSR-NCR

delete(theta)
delete(theta_areave)
delete(theta_latlon)
delete(theta_timeave)
delete(mu)
delete(mut)
delete(mub)
delete(pupx)
delete(pvpy)
delete(divv) ;pupx+pvpy 
delete(pmupt)
delete(msfuy)
delete(msfvx)
delete(muu)
delete(muv)
delete(meanth)
delete(meanw)
delete(pwbpeta)

type="ps"  
wks = gsn_open_wks(type,"profile_CTL-NCR-NIR-NSR")  ;w bar pian theta pian eta
gsn_define_colormap(wks,"WhViBlGrYeOrReWh")
  ;************************************************
; create plot
;************************************************
; title = "Profiles"
 res                   = True                       ; plot mods desired
 ;res@tiMainString      = title                     ; add title
; res@gsnCenterString   = "Time = average" 
;res@trYReverse        = True                       ; reverse Y-axis
; res@tiXAxisString     = "difference in heat budget(~S~o~N~C day~S~-1~N~)"
 res@tiYAxisString     = "Height (km)"
 res@gsnFrame           = False                      ; don't advance frame yet
 res@trYMinF = 1
 res@trYMaxF = 20
 res@trXMinF = -20.
 res@trXMaxF = 20.
; res@tmXBTickSpacingF=2.
; add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend
 
 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .85                   ; move units right
 res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down
 
 res@pmLegendWidthF         = 0.1                  ; Change width and
 res@pmLegendHeightF        = 0.14                  ; height of legend.
 res@lgLabelFontHeightF     = 15                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend

  res@xyLineLabelFontColors  = (/"black","black","black","black","black","black"/)
 res@xyExplicitLegendLabels = (/"ztend","ztend2","~F14~F~B~~F33~q~","~F14~F~B~d~F33~w~","ztend1","ztend3"/)
  res@xyLineColors      = (/"black","red","orange","magenta","blue","green"/)  ; line colors"deepskyblue"
  res@xyLineThicknesses = (/2.5,2.5,2.5,2.5,2.5,2.5/)        ; line thicknesses
  res@xyDashPatterns    = (/0,0,0,0,0,0/)        ; line patterns
 ;asciiwrite("testme.txt",z(0,:,0,0))
z(7,0:33,7,7)=z(7,0:33,7,7)/1000
res=True
res@gsnDraw = False
res@gsnFrame =False
plots = new(3,graphic)

plots(0)=gsn_csm_xy(wks,pt_bias(0:5,0,2:33),z(7,2:33,7,7),res)
plots(1)=gsn_csm_xy(wks,pt_bias(0:5,1,2:33),z(7,2:33,7,7),res)
plots(2)=gsn_csm_xy(wks,pt_bias(0:5,2,2:33),z(7,2:33,7,7),res)
pres=True

pres@gsnPanelFigureStrings =(/"(a)","(b)","(c)"/)
pres@gsnPanelFigureStringsPerimOn=False

gsn_panel(wks,plots,(/2,2/),pres)
 ;plot  = gsn_csm_xy (wks,pt_bias(0:5,4,2:33),z(7,2:33,7,7),res)    ;-(NCR-CTL) 
 	   
; create plot

end
                                       



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
  fn  = "/public/home/wangbin/finalversion/WRFV3/run/wrfout_d03_2015-08-07_12:00:00.nc" ; define filename
  a = addfile(fn,"r")

; We generate plots, but what kind do we prefer?
; type = "x11"
; type = "pdf"
 type = "ps"
; type = "ncgm"
 ; wks = gsn_open_wks(type,"Qvaporadvection")
;  gsn_define_colormap(wks,"BlAqGrYeOrReVi200") ; choose colormap
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  times = wrf_user_getvar(a,"times",-1)  ; get all times in the file
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
  ;  sn_s=1
  ;  sn_e=315
  ;  we_s=1
  ;  we_e=294
    sn_s=110
    sn_e=309
    we_s=7
    we_e=245
    sn=200
    we=239
    mpnum=43 
    grid_num=sn*we  
	step=60
   ; ntimes=(ntimes-1)*step/3600 +1 ;number of hours,because number of steps*step=number of hours*3600
    stepnum=1
    grid_num=sn*we
    nspecies=6     


;fn2  = "/public/home/wangbin/WRFV3-new/run/wrfout_d03.nc" ; define filename
;	a2 = addfile(fn2,"r")
    rain_grid = wrf_user_getvar(a,"RAINNC",-1)
 
    rain_s=new((/ntimes,sn,we/),float)

  do it=0,ntimes-1
    if(it.eq.0)then
      rain_s(it,:,:)=rain_grid(it,sn_s:sn_e,we_s:we_e)
else
     rain_s(it,:,:)=rain_grid(it,sn_s:sn_e,we_s:we_e)-rain_grid(it-1,sn_s:sn_e,we_s:we_e)
end if
end do
 
                     
 
  delete([/rain_grid/]) 
  ;;;;;;;;;;;;;;;;calculate microphysical terms(source and sink of qv,qc,qr,qi,qs,qg)
  ;qv
   revpm = wrf_user_getvar(a,"REVP",-1)
    smltevp = wrf_user_getvar(a,"SMLTEVP",-1)
    gmltevp = wrf_user_getvar(a,"GMLTEVP",-1)
    ssub=wrf_user_getvar(a,"SSUB",-1)
    sdep = wrf_user_getvar(a,"SDEP",-1)
    gsub = wrf_user_getvar(a,"GSUB",-1)
    gdep = wrf_user_getvar(a,"GDEP",-1)
	ladj = wrf_user_getvar(a,"LADJ",-1)
	iadj = wrf_user_getvar(a,"IADJ",-1)
	;qc
	raut = wrf_user_getvar(a,"RAUT",-1)
    racw = wrf_user_getvar(a,"RACW",-1)
    sacw = wrf_user_getvar(a,"SACW",-1)
    sfw=wrf_user_getvar(a,"SFW",-1)
    gacw = wrf_user_getvar(a,"GACW",-1)
    ihom= wrf_user_getvar(a,"IHOM",-1)
	imlt= wrf_user_getvar(a,"IMLT",-1)
	idw= wrf_user_getvar(a,"IDW",-1)
	;qr
	iacr = wrf_user_getvar(a,"IACR",-1)
   sacr = wrf_user_getvar(a,"SACR",-1)
    gfrg = wrf_user_getvar(a,"GFRG",-1)
    gacrgdelt4=wrf_user_getvar(a,"GACRGDELT4",-1)
    gacrpg1sdelt4 = wrf_user_getvar(a,"GACRPG1SDELT4",-1)
    smlt = wrf_user_getvar(a,"SMLT",-1)
    gmltg = wrf_user_getvar(a,"GMLTG",-1)
	;qi
	 saut = wrf_user_getvar(a,"SAUT",-1)
    sfi = wrf_user_getvar(a,"SFI",-1)
    saci = wrf_user_getvar(a,"SACI",-1)
    raci=wrf_user_getvar(a,"RACI",-1)
    gacigdelt4 = wrf_user_getvar(a,"GACIGDELT4",-1)
    gacipg1sdelt4 = wrf_user_getvar(a,"GACIPG1SDELT4",-1)
	;qs
    racidelt3 = wrf_user_getvar(a,"RACIDELT3",-1)
    iacrdelt3 = wrf_user_getvar(a,"IACRDELT3",-1)
    gautg=wrf_user_getvar(a,"GAUTG",-1)
    gacsgdelt4 = wrf_user_getvar(a,"GACSGDELT4",-1)
    gacspg1sdelt4 = wrf_user_getvar(a,"GACSPG1SDELT4",-1)
    racs1sdelt2 = wrf_user_getvar(a,"RACS1SDELT2",-1)
	sacrdelt2 = wrf_user_getvar(a,"SACRDELT2",-1)
	gacsg = wrf_user_getvar(a,"GACSG",-1)
	;qg
    gacwdelt4 = wrf_user_getvar(a,"GACWDELT4",-1)
    gwet1sdelt4 = wrf_user_getvar(a,"GWET1SDELT4",-1)
     gacwgezero=wrf_user_getvar(a,"GACWGE0",-1) ;ppgacwge0
	 sacwgezero=wrf_user_getvar(a,"SACWGE0",-1) ;ppsacwge0
	 sacwltzero=wrf_user_getvar(a,"SACWLT0",-1) ;ppsacwlt0
	
     revp= new((/mpnum,ntimes,south_north,west_east/),float)

    
    do t=0,ntimes-1
	  revp(0,t,:,:)=revpm(t,:,:)
	  revp(1,t,:,:)=smltevp(t,:,:)
	  revp(2,t,:,:)=gmltevp(t,:,:)
	  revp(3,t,:,:)=ssub(t,:,:)
	  revp(4,t,:,:)=sdep(t,:,:)
	  revp(5,t,:,:)=gsub(t,:,:)
      revp(6,t,:,:)=gdep(t,:,:)
	  revp(7,t,:,:)=ladj(t,:,:)
	  revp(8,t,:,:)=iadj(t,:,:)
	  revp(9,t,:,:)=raut(t,:,:)
	  revp(10,t,:,:)=racw(t,:,:) 
	  revp(11,t,:,:)=sacw(t,:,:)
	  revp(12,t,:,:)=sfw(t,:,:)
	  revp(13,t,:,:)=gacw(t,:,:)
	  revp(14,t,:,:)=ihom(t,:,:)
	  revp(15,t,:,:)=imlt(t,:,:)
	  revp(16,t,:,:)=idw(t,:,:)
	  revp(17,t,:,:)=iacr(t,:,:)
	  revp(18,t,:,:)=sacr(t,:,:) 
	  revp(19,t,:,:)= gfrg(t,:,:) 
	  revp(20,t,:,:)= gacrgdelt4(t,:,:)
	  revp(21,t,:,:)=gacrpg1sdelt4(t,:,:)
	  revp(22,t,:,:)=smlt(t,:,:)
	  revp(23,t,:,:)=gmltg(t,:,:)
	  revp(24,t,:,:)=saut(t,:,:)
	  revp(25,t,:,:)=sfi(t,:,:)
	  revp(26,t,:,:)=saci(t,:,:)
	  revp(27,t,:,:)=raci(t,:,:)
	  revp(28,t,:,:)=gacigdelt4(t,:,:)
	  revp(29,t,:,:)=gacipg1sdelt4(t,:,:)
	  revp(30,t,:,:)=racidelt3(t,:,:)
	  revp(31,t,:,:)= iacrdelt3(t,:,:) 
	  revp(32,t,:,:)=gautg(t,:,:)
	  revp(33,t,:,:)=gacsgdelt4(t,:,:)
	  revp(34,t,:,:)=gacspg1sdelt4(t,:,:)
	  revp(35,t,:,:)=racs1sdelt2(t,:,:)
	  revp(36,t,:,:)=sacrdelt2(t,:,:)
	  revp(37,t,:,:)=gacsg(t,:,:)
	  revp(38,t,:,:)=gacwdelt4(t,:,:)
	  revp(39,t,:,:)=gwet1sdelt4(t,:,:)
	  revp(40,t,:,:)=gacwgezero(t,:,:)  ;=wrf_user_getvar(a,"GACWGE0",-1) ;ppgacwge0
	  revp(41,t,:,:)=sacwgezero(t,:,:)   ;=wrf_user_getvar(a,"SACWGE0",-1) ;ppsacwge0
	  revp(42,t,:,:)=sacwltzero(t,:,:)   ;=wrf_user_getvar(a,"SACWLT0",-1) ;ppsacwlt0
    end do
	delete([/smltevp,ssub,sdep,gsub,gdep,ladj,iadj/])	
	delete([/raut,racw,sacw,sfw,gacw,ihom,imlt,idw/])
    delete([/iacr,sacr,gfrg,gacrgdelt4,gacrpg1sdelt4,smlt,gmltg/])
	delete([/saut,sfi,raci,gacigdelt4,gacipg1sdelt4/])
	delete([/racidelt3,iacrdelt3,gautg,gacsgdelt4,gacspg1sdelt4,racs1sdelt2,sacrdelt2,gacsg/])
	delete([/gacwdelt4,gwet1sdelt4,gacwgezero,sacwltzero,sacwgezero/])	
     revps=new((/mpnum,ntimes,sn,we/),float)    
     revps(:,:,:,:)=revp(:,:,sn_s:sn_e,we_s:we_e)
    revps=3600*revps    ;convert s to hour
	
     Sqx=new((/nspecies,ntimes,sn,we/),float)   ;nspecies=6
  

 
  types=8
    Qnc=new((/ntimes,sn,we/),float)
       test=new((/ntimes,sn,we/),float)

	Qcm=new((/ntimes,sn,we/),float)
	Qcmw=new((/ntimes,sn,we/),float)   ;Qcmc+Qcmr
	Qcmi=new((/ntimes,sn,we/),float)   ;Qcmi+Qcms+Qcmg
	Qcmc=new((/ntimes,sn,we/),float)
	Qcmr=new((/ntimes,sn,we/),float)
	Qcmice=new((/ntimes,sn,we/),float)
	Qcms=new((/ntimes,sn,we/),float)
	Qcmg=new((/ntimes,sn,we/),float)
  

	Qcm_total=new((/types/),float)
	Qcmw_total=new((/types/),float)
	Qcmi_total=new((/types/),float)
    
	Ps=new((/types/),float)
	Fc=new((/types/),float)
    qssub=new((/ntimes/),float)
    qsdep=new((/ntimes/),float)
	qrevp=new((/ntimes/),float)
	qgsub=new((/ntimes/),float)
	qgdep=new((/ntimes/),float)
	qsmltevp=new((/ntimes/),float)
	qgmltevp=new((/ntimes/),float)
	qladj=new((/ntimes/),float)
	qiadj=new((/ntimes/),float)
Qcm=0



Qcm_total=0

Ps=0
Fc=0
Qcmw_total=0
Qcmi_total=0



 ;
 qssub(:)=dim_sum_n_Wrap(revps(22,:,:,:),(/1,2/))
    
	
;	print(ps + "/"+cm + "/"+cmc + "/"+cmr + "/"+cmi + "/"+cms + "/"+cmg)
	print(qssub)
end
                                       



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

   stimes=ntimes
    stepnum=60*60/step
    grid_num=sn*we
    nspecies=4    ;number of sensitivity experiments
    ntends=8  
	;number of theta tends:0mp(latent heat),1advection,2pbl(sensible heat)
	;3:localchange,4radiation,5:longwavw,6shortwavw,7diffusion
   p1000mb      = 100000.
   rcp=2/7
   pi_phy=new((/stimes,bottom_top,south_north,west_east/),float)
; First get the variables we will need   
theta=new((/ntends,nspecies,ntimes,bottom_top,south_north,west_east/),float)   
theta=0   
    z   = wrf_user_getvar(a_ctl, "z",-1)         ; grid point height
    p   = wrf_user_getvar(a_ctl, "pressure",-1)         ; grid point height
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
	mu=a_ctl->MU(:,:,:)
    mub=a_ctl->MUB(:,:,:)
    mut=new((/stimes,south_north,west_east/),float)
    mut=mu+mub
    theta(0,0,:,:,:,:)= a_ctl->H_DIABATIC(:,:,:,:)
	theta(1,0,:,:,:,:)= a_ctl->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s
	theta(2,0,:,:,:,:) = a_ctl->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s 	
	theta(4,0,:,:,:,:)= a_ctl->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s
	theta(5,0,:,:,:,:)= a_ctl->RTHRATLW(:,:,:,:)
	theta(6,0,:,:,:,:) = a_ctl->RTHRATSW(:,:,:,:)
	theta(7,0,:,:,:,:) = a_ctl->TDIF_TEND(:,:,:,:)
	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,0,:,ib,:,:)=theta(it,0,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
         do ib=0,bottom_top-1
          theta(1,0,:,ib,:,:)=theta(1,0,:,ib,:,:)/mut(:,:,:)
          theta(2,0,:,ib,:,:)=theta(2,0,:,ib,:,:)/mut(:,:,:)
          theta(4,0,:,ib,:,:)=theta(4,0,:,ib,:,:)/mut(:,:,:)
         end do
	p   = wrf_user_getvar(a_NCR, "pressure",-1)         ; grid point height
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
	mu=a_NCR->MU(:,:,:)
    mub=a_NCR->MUB(:,:,:)
    mut=new((/stimes,south_north,west_east/),float)
    mut=mu+mub
    theta(0,1,:,:,:,:)= a_NCR->H_DIABATIC(:,:,:,:);MICROPHYSICS LATENT HEATING 
	theta(1,1,:,:,:,:)= a_NCR->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s ADVECTIVE POTENTIAL TEMPERATURE TENDENCY 
	theta(2,1,:,:,:,:) = a_NCR->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s COUPLED THETA TENDENCY DUE TO PBL PARAMETERIZATION 	
	theta(4,1,:,:,:,:)= a_NCR->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s THETA TENDENCY DUE TO RADIATION 
	theta(5,1,:,:,:,:)= a_NCR->RTHRATLW(:,:,:,:);UNCOUPLED THETA TENDENCY DUE TO LONG WAVE RADIATION 
	theta(6,1,:,:,:,:) = a_NCR->RTHRATSW(:,:,:,:);UNCOUPLED THETA TENDENCY DUE TO SHORT WAVE RADIATION 
	theta(7,1,:,:,:,:) = a_NCR->TDIF_TEND(:,:,:,:) ;DIFFUSION POTENTIAL TEMPERATURE TENDENCY 
	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,1,:,ib,:,:)=theta(it,1,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
        do ib=0,bottom_top-1
          theta(1,1,:,ib,:,:)=theta(1,1,:,ib,:,:)/mut(:,:,:)
          theta(2,1,:,ib,:,:)=theta(2,1,:,ib,:,:)/mut(:,:,:)
          theta(4,1,:,ib,:,:)=theta(4,1,:,ib,:,:)/mut(:,:,:)
         end do
	
	p   = wrf_user_getvar(a_NIR, "pressure",-1)         ; grid point height
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
	mu=a_NIR->MU(:,:,:)
    mub=a_NIR->MUB(:,:,:)
    mut=new((/stimes,south_north,west_east/),float)
    mut=mu+mub
    theta(0,2,:,:,:,:)= a_NIR->H_DIABATIC(:,:,:,:)
	theta(1,2,:,:,:,:)= a_NIR->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s
	theta(2,2,:,:,:,:) = a_NIR->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s 	
	theta(4,2,:,:,:,:)= a_NIR->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s
	theta(5,2,:,:,:,:)= a_NIR->RTHRATLW(:,:,:,:)
	theta(6,2,:,:,:,:) = a_NIR->RTHRATSW(:,:,:,:)
	theta(7,2,:,:,:,:) = a_NIR->TDIF_TEND(:,:,:,:)
	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,2,:,ib,:,:)=theta(it,2,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
        do ib=0,bottom_top-1
          theta(1,2,:,ib,:,:)=theta(1,2,:,ib,:,:)/mut(:,:,:)
          theta(2,2,:,ib,:,:)=theta(2,2,:,ib,:,:)/mut(:,:,:)
          theta(4,2,:,ib,:,:)=theta(4,2,:,ib,:,:)/mut(:,:,:)
         end do
   	
    p   = wrf_user_getvar(a_NSR, "pressure",-1)         ; grid point height
	p_phy=100*p
    pi_phy(:,:,:,:) = (p_phy(:,:,:,:)/p1000mb)^rcp   ;t=t0thz,temperature/theta
	mu=a_NSR->MU(:,:,:)
    mub=a_NSR->MUB(:,:,:)
    mut=new((/stimes,south_north,west_east/),float)
    mut=mu+mub
    theta(0,3,:,:,:,:)= a_NSR->H_DIABATIC(:,:,:,:)
	theta(1,3,:,:,:,:)= a_NSR->RTHFTEN(:,:,:,:)    ;coupled by mut,unit is Pa K/s
	theta(2,3,:,:,:,:) = a_NSR->RTHBLTEN(:,:,:,:) ;coupled by mut,unit is Pa K/s 	
	theta(4,3,:,:,:,:)= a_NSR->RTHRATEN(:,:,:,:);coupled by mut,unit is Pa K/s
	theta(5,3,:,:,:,:)= a_NSR->RTHRATLW(:,:,:,:)
	theta(6,3,:,:,:,:) = a_NSR->RTHRATSW(:,:,:,:)	
	theta(7,3,:,:,:,:) = a_NSR->TDIF_TEND(:,:,:,:)
	do it=0,ntends-1
	 do ib=0,bottom_top-1
	  theta(it,3,:,ib,:,:)=theta(it,3,:,ib,:,:)*pi_phy(:,ib,:,:)   ;convert from theta to temperature
	 end do
	end do
         do ib=0,bottom_top-1
          theta(1,3,:,ib,:,:)=theta(1,3,:,ib,:,:)/mut(:,:,:)
          theta(2,3,:,ib,:,:)=theta(2,3,:,ib,:,:)/mut(:,:,:)
          theta(4,3,:,ib,:,:)=theta(4,3,:,ib,:,:)/mut(:,:,:)
         end do

	theta(3,:,:,:,:,:)= theta(0,:,:,:,:,:)+theta(1,:,:,:,:,:)+theta(2,:,:,:,:,:)+theta(4,:,:,:,:,:)  ;local change
	
  theta_latlon=new((/ntends,nspecies,ntimes,bottom_top,sn,we/),float) 
  theta_latlon(:,:,:,:,:,:)=theta(:,:,:,:,sn_s:sn_e,we_s:we_e)
  
   theta_areave=new((/ntends,nspecies,ntimes,bottom_top/),float) 
   theta_areave=dim_sum_n_Wrap(theta_latlon(:,:,:,:,:,:),(/4,5/))  ;sum over lat/lon
   theta_areave = theta_areave/grid_num     ;per grid
   
  theta_timeave=new((/ntends,nspecies,bottom_top/),float)    ;average over lat/lon
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
   theta_timeave=theta_timeave
;ntends=8:0-mp(latent heat),1-advection,2-pbl(sensible heat),3-:localchange,4-radiation,5-longwavw,6-shortwavw,7-diffussion
;tc_bias=new((/ntends,nspecies-1,bottom_top/),float) 
;tc_bias(:,0,:)=theta_timeave(:,0,:)-theta_timeave(:,1,:)   ;CTL-NCR
;tc_bias(:,1,:)=theta_timeave(:,0,:)-theta_timeave(:,2,:)   ;CTL-NIR
;tc_bias(:,2,:)=theta_timeave(:,0,:)-theta_timeave(:,3,:)   ;CTL-NSR
;tc_bias(:,3,:)=theta_timeave(:,1,:)-theta_timeave(:,3,:)   ;NCR-NWR
;tc_bias(:,4,:)=0     ;line:x=0
pt_bias=new((/ntends-2,nspecies-1,bottom_top/),float)
do i=0,ntends-4
 pt_bias(i,0,:)=theta_timeave(i,1,:)-theta_timeave(i,0,:)   ;NCR-CTL
 pt_bias(i,1,:)=theta_timeave(i,2,:)-theta_timeave(i,0,:)   ;NIR-CTL
 pt_bias(i,2,:)=theta_timeave(i,3,:)-theta_timeave(i,0,:)   ;NSR-CTL
 ;pt_bias(i,3,:)=theta_timeave(i,1,:)-theta_timeave(i,3,:)   ;NCR-NWR
 ;pt_bias(i,4,:)=theta_timeave(i,1,:)-theta_timeave(i,0,:)   ;NCR-CTL 
end do
 pt_bias(ntends-3,:,:)=0
asciiwrite("hdctl.txt",theta_timeave(1,0,:))
asciiwrite("hd2.txt",theta_timeave(1,2,:))
asciiwrite("hd.txt",theta_timeave(1,1,:))
asciiwrite("hd3.txt",theta_timeave(1,3,:))
delete(theta)
delete(theta_areave)
delete(theta_latlon)
delete(theta_timeave)
delete(mu)
delete(mut)
delete(mub)

;tc_bias=tc_bias*3600    ;convert from celsius per minute to celsius per hour
pt_bias=pt_bias*3600*24    ;convert from celsius per second to celsius per day

pt_bias = -pt_bias  ;;reverse to 0:CTL-NWR,1:NIR-NCR,2:CTL-NIR,3:NWR-NCR  4:CTL-NCR

type="ps"  
wks = gsn_open_wks(type,"profile_CTL-NCR")
gsn_define_colormap(wks,"WhViBlGrYeOrReWh")
  ;************************************************
; create plot
;************************************************
; title = "Profiles"
 res                   = True                       ; plot mods desired
; res@tiMainString      = title                     ; add title
; res@gsnCenterString   = "Time = average" 
;res@trYReverse        = True                       ; reverse Y-axis
; res@tiXAxisString     = "difference in heat budget(~S~o~N~C day~S~-1~N~)"
 res@tiYAxisString     = "Height (km)"
 res@gsnFrame           = False                      ; don't advance frame yet
   res@trYMinF = 1
 res@trYMaxF = 20
  res@trXMinF = -6.
 res@trXMaxF = 6.
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
txres=True
txres@txFontHeightF=0.015
gsn_text_ndc(wks,"(a)",0.22,0.78,txres)
;gsn_text_ndc(wks,"d02",0.4,0.42,txres)
;gsn_text_ndc(wks,"d03",0.5,0.48,txres)
  res@xyLineLabelFontColors  = (/"black","black","black","black","black"/)
 res@xyExplicitLegendLabels = (/"Fmp","Fhd","Fpbl","Floc","Frad"/)
 ; res@xyExplicitLegendLabels = (/"CTL-NCR","CTL-NIR","CTL-NSR","NCR-NWR"/)
  res@xyLineColors      = (/"black","red","green","blue","orange"/)  ; line colors
  res@xyLineThicknesses = (/2.5,2.5,2.5,2.5,2.5/)        ; line thicknesses
  res@xyDashPatterns    = (/0,0,0,0.0,0.0/)        ; line patterns
 ;asciiwrite("testme.txt",z(0,:,0,0))
z(7,0:33,7,7)=z(7,0:33,7,7)/1000
 plot  = gsn_csm_xy (wks,pt_bias(0:4,0,2:33),z(7,2:33,7,7),res)            
; create plot

end
                                       


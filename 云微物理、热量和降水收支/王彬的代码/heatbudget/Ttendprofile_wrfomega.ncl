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
    stimes=ntimes
    stepnum=60*60/step
    grid_num=sn*we
    nspecies=4    ;number of sensitivity experiments:0-ctl,1-NCR,2-NIR,3-NSR
    ntends=3+2+1 ;;0:temperature,1:vertical velocity,2:relative humidity,3-qvapor,4-pressure,5-omega 
	
   z   = wrf_user_getvar(a_ctl, "z",-1)         ; grid point height
	;number of theta tends:0mp(latent heat),1advection,2pbl(sensible heat)
	;3:localchange,4radiation,5:longwavw,6shortwavw,7diffusion
;;;theta used to represent all the varaibles,for I did not find the right name for all the variables.
; First get the variables we will need   
theta=new((/ntends,nspecies,ntimes,bottom_top,south_north,west_east/),float)   
theta=0   

    theta(0,0,:,:,:,:)= wrf_user_getvar(a_ctl,"tc",-1)+273.15   ; Temperature [C]
	theta(1,0,:,:,:,:)= wrf_user_getvar(a_ctl,"wa",-1)   ;vertical velocity,W component of wind on mass points
	theta(2,0,:,:,:,:) = wrf_user_getvar(a_ctl,"rh",-1) ;Relative Humidity [%]
	theta(3,0,:,:,:,:) = wrf_user_getvar(a_ctl,"QVAPOR",-1)  ;qvapor
	theta(4,0,:,:,:,:) = wrf_user_getvar(a_ctl,"pressure",-1)  ;pressure	
	
    theta(0,1,:,:,:,:)= wrf_user_getvar(a_NCR,"tc",-1)+273.15   ; Temperature [C]
	theta(1,1,:,:,:,:)= wrf_user_getvar(a_NCR,"wa",-1)   ;vertical velocity,W component of wind on mass points
	theta(2,1,:,:,:,:) = wrf_user_getvar(a_NCR,"rh",-1) ;Relative Humidity [%] 	
	theta(3,1,:,:,:,:) = wrf_user_getvar(a_NCR,"QVAPOR",-1)  ;qvapor
	theta(4,1,:,:,:,:) = wrf_user_getvar(a_NCR,"pressure",-1)  ;pressure
	
	theta(0,2,:,:,:,:)= wrf_user_getvar(a_NIR,"tc",-1)+273.15   ; Temperature [C]
	theta(1,2,:,:,:,:)= wrf_user_getvar(a_NIR,"wa",-1)   ;vertical velocity,W component of wind on mass points
	theta(2,2,:,:,:,:) = wrf_user_getvar(a_NIR,"rh",-1) ;Relative Humidity [%] 	
	theta(3,2,:,:,:,:) = wrf_user_getvar(a_NIR,"QVAPOR",-1)  ;qvapor
	theta(4,2,:,:,:,:) = wrf_user_getvar(a_NIR,"pressure",-1)  ;pressure
	
	theta(0,3,:,:,:,:)= wrf_user_getvar(a_NSR,"tc",-1)+273.15   ; Temperature [C]
	theta(1,3,:,:,:,:)= wrf_user_getvar(a_NSR,"wa",-1)   ;vertical velocity,W component of wind on mass points
	theta(2,3,:,:,:,:) = wrf_user_getvar(a_NSR,"rh",-1) ;Relative Humidity [%] 	
	theta(3,3,:,:,:,:) = wrf_user_getvar(a_NSR,"QVAPOR",-1)  ;qvapor
	theta(4,3,:,:,:,:) = wrf_user_getvar(a_NSR,"pressure",-1)  ;pressure hPa	
	
	theta(5,:,:,:,:,:) = wrf_omega(theta(3,:,:,:,:,:),theta(0,:,:,:,:,:),theta(1,:,:,:,:,:),theta(4,:,:,:,:,:))	;;qv,tk,w,p
	;return_val@units = "Pa/s"
  theta_latlon=new((/ntends,nspecies,ntimes,bottom_top,sn,we/),float) ;nspecies=4:number of sensitivity tests 0:ctl,1:NCR,2:NIR,3:NSR
  theta_latlon(:,:,:,:,:,:)=theta(:,:,:,:,sn_s:sn_e,we_s:we_e)
  
   theta_areave=new((/ntends,nspecies,ntimes,bottom_top/),float) 
   theta_areave=dim_sum_n_Wrap(theta_latlon(:,:,:,:,:,:),(/4,5/))  ;sum over lat/lon
   theta_areave = theta_areave/grid_num                           ;per grid

 ;ntends=3;;0:temperature,1:vertical velocity,2:relative humidity  
;;nspecies=4 ;number of sensitivity experiments:0-ctl,1-NCR,2-NIR,3-NSR 
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
need=new((/ntends,nspecies+1,bottom_top/),float) 
need(:,0:3,:)=theta_timeave(:,:,:)
need(:,4,:)=0

tc_bias=new((/ntends,nspecies+1+1,bottom_top/),float) ;ntends=3;;0:temperature,1:vertical velocity,2:relative humidity  
tc_bias(:,0,:)=theta_timeave(:,1,:)-theta_timeave(:,0,:)   ;NSR-CTL
tc_bias(:,1,:)=theta_timeave(:,2,:)-theta_timeave(:,0,:)   ;NCR-NIR
tc_bias(:,2,:)=theta_timeave(:,3,:)-theta_timeave(:,0,:)   ;NIR-CTL
tc_bias(:,3,:)=theta_timeave(:,1,:)-theta_timeave(:,3,:)   ;NCR-NSR
tc_bias(:,4,:)=theta_timeave(:,1,:)-theta_timeave(:,0,:)   ;NCR-CTL
tc_bias(:,5,:)=0     ;line:x=0

tc_bias=-(-tc_bias)   ;;convert omega to -omega (similar to W)
;;convert to from  "NSR-CTL","NCR-NIR","NIR-CTL","NCR-NSR" to "CTL-NSR","NIR-NCR","CTL-NIR","NSR-NCR"
asciiwrite("omegactl.txt",theta_timeave(5,0,:))
delete(theta)
delete(theta_areave)
delete(theta_latlon)
delete(theta_timeave)


;;tc_bias=tc_bias*3600*24   ;convert from celsius per second to celsius per day

type="ps"  
wks = gsn_open_wks(type,"profile-1000omega_ctl_NCR_NIR_NSR+CTL-NCR") ;;get the needed w profile
gsn_define_colormap(wks,"WhViBlGrYeOrReWh")
  ;************************************************
; create plot
;************************************************
; title = "Profiles_W"
 res                   = True                       ; plot mods desired
; res@tiMainString      = title                     ; add title
; res@gsnCenterString   = "Time = average" 
;res@trYReverse        = True                       ; reverse Y-axis
; res@tiXAxisString     = "omega difference(ms~S~-1~N~)"  ;;(~S~o~N~C s~S~-1~N~)"
 res@tiYAxisString     = "Height (km)"
 res@gsnFrame           = False                      ; don't advance frame yet
   res@trYMinF = 1
  res@trYMaxF = 20
 ;; res@trXMinF = -4.
;; res@trXMaxF = 5.
; add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend
 txres=True
txres@txFontHeightF=0.015
gsn_text_ndc(wks,"(b)",0.22,0.78,txres)
 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .85                   ; move units right
 res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down
 res@gsnMaximize = True
 res@pmLegendWidthF         = 0.1                  ; Change width and
 res@pmLegendHeightF        = 0.14                  ; height of legend.
 res@lgLabelFontHeightF     = 15                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend
res@gsnXRefLine           = 0.0             ; draw a reference line
  res@xyLineLabelFontColors  = (/"black","black","black"/)
  ;;res@xyExplicitLegendLabels = (/"CTL","NCR","NIR","NSR"," "/)
  res@xyExplicitLegendLabels = (/"CTL-NCR","CTL-NIR","CTL-NSR"/)
  res@xyLineColors      = (/"black","red","blue"/)  ; line colors
  res@xyLineThicknesses = (/2.5,2.5,2.5/)        ; line thicknesses
  res@xyDashPatterns    = (/0,0,0/)        ; line patterns
 ;asciiwrite("testme.txt",z(0,:,0,0))
z(7,0:33,7,7)=z(7,0:33,7,7)/1000
tc_bias(5,0:2,2:33)=1000*tc_bias(5,0:2,2:33)  ;;convert unit Pa/s to 1000 Pa/s
 plot  = gsn_csm_xy (wks,tc_bias(5,0:2,2:33),z(7,2:33,7,7),res)            
; create plot

end
                                       


load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


f1 = addfile("ctl.nc", "r")
f2 = addfile("ice.nc", "r")

lat = f1->lat
lon = f1->lon
level = f1->lev

uu_ctl = f1->uwnd(0:2639,:,:,:)
vv_ctl = f1->vwnd(0:2639,:,:,:)
q_ctl = f1->q(0:2639,:,:,:)

uu_snt = f2->uwnd(0:2639,:,:,:)
vv_snt = f2->vwnd(0:2639,:,:,:)
q_snt = f2->q(0:2639,:,:,:)

uq_ctl = uu_ctl * q_ctl * 1000.0
vq_ctl = vv_ctl * q_ctl * 1000.0

dpp = new(dimsizes(level)-1, "float")
dpp = level(0:dimsizes(level)-2) - level(1:dimsizes(level)-1)

uq_new_ctl = new((/2640,dimsizes(level)-1,dimsizes(lat),dimsizes(lon)/), "float")
uq_new_ctl = (uq_ctl(:,0:dimsizes(level)-2,:,:)+uq_ctl(:,1:dimsizes(level)-1,:,:))/2.0

dp = conform(uq_new_ctl,dpp,1)

uq_int_ctl = dim_sum_n_Wrap(uq_new_ctl*dp,1) / 9.8
copy_VarMeta(uu_ctl(:,0,:,:),uq_int_ctl)

vq_new_ctl = new((/2640,dimsizes(level)-1,dimsizes(lat),dimsizes(lon)/), "float")
vq_new_ctl = (vq_ctl(:,0:dimsizes(level)-2,:,:)+vq_ctl(:,1:dimsizes(level)-1,:,:))/2.0

vq_int_ctl = dim_sum_n_Wrap(vq_new_ctl*dp,1) / 9.8
copy_VarMeta(uu_ctl(:,0,:,:),vq_int_ctl)

qVdiv_int_ctl = uv2dv_cfd(uq_int_ctl,vq_int_ctl,lat,lon,3)
copy_VarMeta(uu_ctl(:,0,:,:),qVdiv_int_ctl)

; snt

uq_snt = uu_snt * q_snt * 1000.0
vq_snt = vv_snt * q_snt * 1000.0

uq_new_snt = new((/2640,dimsizes(level)-1,dimsizes(lat),dimsizes(lon)/), "float")
uq_new_snt = (uq_snt(:,0:dimsizes(level)-2,:,:)+uq_snt(:,1:dimsizes(level)-1,:,:))/2.0

uq_int_snt = dim_sum_n_Wrap(uq_new_snt*dp,1) / 9.8
copy_VarMeta(uu_ctl(:,0,:,:),uq_int_snt)

vq_new_snt = new((/2640,dimsizes(level)-1,dimsizes(lat),dimsizes(lon)/), "float")
vq_new_snt = (vq_snt(:,0:dimsizes(level)-2,:,:)+vq_snt(:,1:dimsizes(level)-1,:,:))/2.0

vq_int_snt = dim_sum_n_Wrap(vq_new_snt*dp,1) / 9.8
copy_VarMeta(uu_ctl(:,0,:,:),vq_int_snt)

qVdiv_int_snt = uv2dv_cfd(uq_int_snt,vq_int_snt,lat,lon,3)
copy_VarMeta(uu_ctl(:,0,:,:),qVdiv_int_snt)


system("rm -rf watflux.nc")
setfileoption("nc","Format","LargeFile")
out = addfile("watflux.nc","c")
out->uq_int_ctl = uq_int_ctl
out->vq_int_ctl = vq_int_ctl
out->qVdiv_int_ctl = qVdiv_int_ctl

out->uq_int_snt = uq_int_snt
out->vq_int_snt = vq_int_snt
out->qVdiv_int_snt = qVdiv_int_snt

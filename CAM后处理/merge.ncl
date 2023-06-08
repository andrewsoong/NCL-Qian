load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  

case="ice"

filename = systemfunc("ls "+case+"/"+"*.nc")

HGT=new((/2640,17,96,144/),"float")
UWND=new((/2640,17,96,144/),"float")
VWND=new((/2640,17,96,144/),"float")
TK=new((/2640,17,96,144/),"float")
OMEGA=new((/2640,17,96,144/),"float")
PRE=new((/2640,96,144/),"float")
Q=new((/2640,17,96,144/),"float")

do idx=0, dimsizes(filename)-1
fall = addfile(filename(idx), "r")

hgt = fall->hgt
time = fall->time
lev = fall->lev
lat = fall->lat
lon = fall->lon
uwnd = fall->uwnd
vwnd = fall->vwnd
tk = fall->tk
omega = fall->omega
pre = fall->pre
q = fall->q

HGT(idx,:,:,:) = hgt(0,:,:,:)
UWND(idx,:,:,:) = uwnd(0,:,:,:)
VWND(idx,:,:,:) = vwnd(0,:,:,:)
TK(idx,:,:,:) = tk(0,:,:,:)
OMEGA(idx,:,:,:) = omega(0,:,:,:)
PRE(idx,:,:) = pre(0,:,:)
Q(idx,:,:,:) = q(0,:,:,:)
end do


setfileoption("nc","Format","LargeFile")
out = addfile(case+".nc", "c")
out->hgt = HGT
out->uwnd = UWND
out->vwnd = VWND
out->tk = TK
out->omega = OMEGA
out->pre = PRE
out->q = Q

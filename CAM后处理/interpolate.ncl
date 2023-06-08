load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  

case="ice"

filename = systemfunc("ls "+case+"/"+case+".long.cam.h0.*.nc")

do idx=0, dimsizes(filename)-1
myname = str_get_cols(filename(idx),20,30)

fall = addfile(filename(idx), "r")

Z3 = fall->Z3
T = fall->T
U = fall->U
V = fall->V
OMG = fall->OMEGA
PRE = fall->PRECT
SPH = fall->Q

PS = fall->PS

hbcofa1=fall->hyai
hbcofb1=fall->hybi

plevo=(/1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10/)

intyp=3
p0=1000.0
kxtrp=False
ii=1

HGT = vinth2p(Z3,hbcofa1,hbcofb1,plevo,PS,intyp,p0,ii,kxtrp)
TK = vinth2p(T,hbcofa1,hbcofb1,plevo,PS,intyp,p0,ii,kxtrp)
UWND = vinth2p(U,hbcofa1,hbcofb1,plevo,PS,intyp,p0,ii,kxtrp)
VWND = vinth2p(V,hbcofa1,hbcofb1,plevo,PS,intyp,p0,ii,kxtrp)
OMEGA = vinth2p(OMG,hbcofa1,hbcofb1,plevo,PS,intyp,p0,ii,kxtrp)
Q = vinth2p(SPH,hbcofa1,hbcofb1,plevo,PS,intyp,p0,ii,kxtrp)

HGT!1="lev"
TK!1="lev"
UWND!1="lev"
VWND!1="lev"
OMEGA!1="lev"
Q!1="lev"

setfileoption("nc","Format","LargeFile")
out = addfile("./extract/"+case+"/"+myname+".nc", "c")
out->hgt = HGT
out->uwnd = UWND
out->vwnd = VWND
out->tk = TK
out->omega = OMEGA
out->pre = PRE
out->q = Q

delete(HGT)
delete(TK)
delete(UWND)
delete(VWND)
delete(OMEGA)
delete(Q)
delete(out)
delete(plevo)
delete(intyp)
delete(p0)
delete(kxtrp)
delete(ii)
delete(hbcofa1)
delete(hbcofb1)
delete(Z3)
delete(T)
delete(U)
delete(V)
delete(OMG)
delete(PRE)
delete(SPH)

print(filename(idx))

end do

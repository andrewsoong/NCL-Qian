C NCLFORTSTART
	SUBROUTINE latlon2xy(grdlat,grdlon,grdj,grdi,np,
     +  xlat_cnt,xlon_cnt,iy,ix,trulat1,trulat2,xlon_std,
     +  dx,dy)
        real grdlat(np),grdlon(np),grdi(np),grdj(np)
C NCLEND
  	real pi, pi2, pi4, d2r, r2d, radius
	real gcon,ogcon,ahem,deg,cn1,cn2,cn3,
     +  cn4,rih,xih,yih,rrih,check
  	real alnfix,alon,x,y
  	real vals(0:8)
        vals(0) = xlat_cnt
        vals(1) = xlon_cnt
        vals(2) = float(ix+1)/2
        vals(3) = float(iy+1)/2
        vals(4) = trulat1      
        vals(5) = trulat2      
        vals(6) = xlon_std
        vals(7) = dx*1000.
        vals(8) = dy*1000.
  	pi = 4.0*atan(1.0)
  	pi2 = pi/2.0
  	pi4 = pi/4.0
  	d2r = pi/180.0
  	r2d = 180.0/pi
  	radius = 6371000.0

C*case where standard lats are the same */

  	if(vals(4) .eq.vals(5)) then
    	gcon = sin(vals(4)*d2r)
	 else
    	gcon = (log(sin((90.0-vals(4))*d2r))
     +  -log(sin((90.0-vals(5))*d2r)))
     +  /(log(tan((90.0-vals(4))*0.5*d2r))
     +  -log(tan((90.0-vals(5))*0.5*d2r)))
	endif
  	ogcon = 1.0/gcon
 	ahem = abs(vals(4))/(vals(4))
  	deg = (90.0-abs(vals(4)))*d2r
  	cn1 = sin(deg)
  	cn2 = radius*cn1*ogcon
  	deg = deg*0.5
  	cn3 = tan(deg)
        deg = (90.0-abs(vals(0)))*0.5*d2r
  	cn4 = tan(deg)
  	rih = cn2*(cn4/cn3)**gcon
  	deg = (vals(1)-vals(6))*d2r*gcon
  	xih = rih*sin(deg)
  	yih = -rih*cos(deg)*ahem
        do i = 1, np
  	deg = (90.0-grdlat(i)*ahem)*0.5*d2r
  	cn4 = tan(deg)
  	rrih = cn2*((cn4/cn3)**gcon)
  	check = 180.0-vals(6)
  	alnfix = vals(6)+check
  	alon = grdlon(i)+check
  	if (alon.lt.0.0) alon = alon+360.0
  	if (alon.gt.360.0) alon = alon-360.0
  	deg = (alon-alnfix)*gcon*d2r
  	x = rrih*sin(deg)
  	y = -rrih*cos(deg)*ahem
  	grdi(i) = vals(2)+(x-xih)/(vals(7))
  	grdj(i) = vals(3)+(y-yih)/(vals(8))
        end do
        RETURN
	END

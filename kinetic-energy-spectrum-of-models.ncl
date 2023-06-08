; ------------------------------------------------------------------------------
; Description:
;
;   This script is used to calculate and plot the kinetic energy spectrum of the
;   model data (contain U and V variables).
;
;   Original scripts and codes are from Haiyang Yu.
;
; Authors:
;
;   Li Dong <dongli@lasg.iap.ac.cn>
;   Haiyang Yu <yuhaiyang1001student@126.com>
; ------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function get_energy_spectrum_dims(num_time, num_lev, num_lon)
begin

    trunc_num = (num_lon-1)/3

    return (/num_time,3,num_lev,trunc_num+1/)

end

procedure calc_energy_spectrum(rad, U, V, spectrum, verbose)
local num_lon, num_lat, num_lev, lon, lat, num_glat, glat, \
    k, u, v, div, vrt, im, in
begin

    dims = dimsizes(U)
    num_lev = dims(0)
    num_lat = dims(1)
    num_lon = dims(2)

    lon = U&lon
    lat = U&lat

    ; --------------------------------------------------------------------------
    if (verbose) then
        system("echo [Notice]: Interpolate data onto Gaussian grid.")
    end if
    num_glat = num_lon/2
    gaus_grid = gaus(num_glat/2)
    glat = gaus_grid(:,0)
    glat@units = "degree_north"
    gwgt = gaus_grid(:,1)
    u = new((/num_lev,num_glat,num_lon/), "float")
    v = new((/num_lev,num_glat,num_lon/), "float")
    do k = 0, num_lev-1
        u(k,:,:) = linint2(lon, lat, U(k,:,:), True, lon, glat, 0)
        v(k,:,:) = linint2(lon, lat, V(k,:,:), True, lon, glat, 0)
    end do
    u!0 = "lev"
    u!1 = "lat"
    u!2 = "lon"
    u&lat = glat
    u&lon = lon
    v!0 = "lev"
    v!1 = "lat"
    v!2 = "lon"
    v&lat = glat
    v&lon = lon
    ; --------------------------------------------------------------------------
    if (verbose) then
        system("echo [Notice]: Calculate divergence.")
    end if
    div = uv2dvG_Wrap(u, v)
    ; --------------------------------------------------------------------------
    if (verbose) then
        system("echo [Notice]: Calculate vorticity.")
    end if
    vrt = uv2vrG_Wrap(u, v)
    ; --------------------------------------------------------------------------
    if (verbose) then
        system("echo [Notice]: Calculate spherical harmonic of both components")
    end if
    coefd = shagC(div)
    coefr = shagC(vrt)
    trunc_num = (num_lon-1)/3
    coefd = tri_trunC(coefd, trunc_num)
    coefr = tri_trunC(coefr, trunc_num)
    ; --------------------------------------------------------------------------
    if (verbose) then
        system("echo [Notice]: Calculate kinetic energy spectrum")
    end if
    ; kinetic energy spectrum
    ; 0 - divergent part; 1 - rotational part; 2 - total
    spectrum(0,:,0) = 0.0
    spectrum(1,:,0) = 0.0
    do in = 1, trunc_num
        spectrum(0,:,in) = coefd(0,:,in,0)^2+coefd(1,:,in,0)^2
        spectrum(1,:,in) = coefr(0,:,in,0)^2+coefr(1,:,in,0)^2
        do im = 1, in
            spectrum(0,:,in) = spectrum(0,:,in)+2* \
                (coefd(0,:,in,im)^2+coefd(1,:,in,im)^2)
            spectrum(1,:,in) = spectrum(1,:,in)+2* \
                (coefr(0,:,in,im)^2+coefr(1,:,in,im)^2)
        end do
        spectrum(0,:,in) = spectrum(0,:,in)*rad*rad/(4.0*in*(in+1))
        spectrum(1,:,in) = spectrum(1,:,in)*rad*rad/(4.0*in*(in+1))
    end do
    spectrum(2,:,:) = spectrum(0,:,:)+spectrum(1,:,:)

end

procedure plot_energy_spectrum(fig_name, spectrum)
begin

    wks = gsn_open_wks("pdf", fig_name)
    ; ==========================================================================
    res = True
    res@gsnFrame = False
    res@trXLog = True
    res@trYLog = True
    res@trXMinF = 1.0
    res@trXMaxF = 400.0
    res@trYMinF = 1.0e-4
    res@trYMaxF = 1.0e5
    res@xyLineThicknessF = 3.0
    res@tiXAxisString = "Wave number"
    res@tiYAxisString = "Kinetic energy (m2 s-2)"
    res@tiMainString = "Kinetic Energy Spectrum"
    ; ==========================================================================
    ; reference lines (power -3 and power -5/3)
    x1 = 2.7
    x2 = 27.0
    x3 = 400.0
    y1 = 1.0e5
    y2 = exp(-3.0*log(x2/x1))*y1
    y3 = exp(-5.0/3.0*log(x3/x2))*y2
    res@xyLineColors = (/"grey"/)
    plot = gsn_csm_xy(wks, (/x1,x2/), (/y1,y2/), res)
    plot = gsn_csm_xy(wks, (/x2,x3/), (/y2,y3/), res)
    ; ==========================================================================
    ; kinetic energy spectrum
    trunc_num = dimsizes(spectrum(0,0,0,:))-1
    res@xyLineColors = (/"blue"/)
    res@xyDashPatterns = (/2/)
    plot = gsn_csm_xy(wks, ispan(1, trunc_num, 1), \
        dim_avg_n_Wrap(spectrum(:,0,:,1:trunc_num), (/0,1/)), res)
    res@xyDashPatterns = (/1/)
    plot = gsn_csm_xy(wks, ispan(1, trunc_num, 1), \
        dim_avg_n_Wrap(spectrum(:,1,:,1:trunc_num), (/0,1/)), res)
    res@xyDashPatterns = (/0/)
    res@gsnDraw = False
    plot = gsn_csm_xy(wks, ispan(1, trunc_num, 1), \
        dim_avg_n_Wrap(spectrum(:,2,:,1:trunc_num), (/0,1/)), res)
    ; ==========================================================================
    ; add legend
    lgres = True
    lgres@lgDashIndexes = (/2,1,0/)
    lgres@lgLabelFontHeightF = 0.1
    lgres@vpWidthF = 0.15
    lgres@vpHeightF = 0.1
    lgres@lgPerimOn = False
    labels = (/"Divergent","Rotational","Total"/)
    legend = gsn_create_legend(wks, 3, labels, lgres)
    amres = True
    amres@amParallelPosF = 0.27
    amres@amOrthogonalPosF = -0.35
    annota = gsn_add_annotation(plot, legend, amres)

    draw(plot)
    frame(wks)

end

function calc_for_hybrid_coordinate(file_names, rad, pres)
local f, num_file, num_time, i, U, V, spectrum
begin

    ; --------------------------------------------------------------------------
    f = addfiles(file_names, "r")
    num_file = dimsizes(file_names)
    num_time = 0
    do i = 0, num_file-1
        num_time = num_time+dimsizes(f[i]->time)
    end do
    ; --------------------------------------------------------------------------
    spectrum = new(get_energy_spectrum_dims(num_time, dimsizes(pres), \
        dimsizes(f[0]->lon)), "float")
    k = 0
    do i = 0, num_file-1
        system("echo [Notice]: Use file "+file_names(i)+".")
        do j = 0, dimsizes(f[i]->time)-1
            ; ==================================================================
            ; interpolate model data onto specific levels
            if (isfilevar(f[i], "P0")) then
                U = vinth2p(f[i]->U(j,:,:,:), f[i]->hyam, f[i]->hybm, pres, \
                    f[i]->PS(j,:,:), 1, f[i]->P0, 1, False)
                V = vinth2p(f[i]->V(j,:,:,:), f[i]->hyam, f[i]->hybm, pres, \
                    f[i]->PS(j,:,:), 1, f[i]->P0, 1, False)
            else
                U = vinth2p(f[i]->U(j,:,:,:), f[i]->hyam, f[i]->hybm, pres, \
                    f[i]->PS(j,:,:), 1, 1000.0, 1, False)
                V = vinth2p(f[i]->V(j,:,:,:), f[i]->hyam, f[i]->hybm, pres, \
                    f[i]->PS(j,:,:), 1, 1000.0, 1, False)
            end if
            ; ==================================================================
            calc_energy_spectrum(rad, U, V, spectrum(k,:,:,:), False)
            k = k+1
        end do
    end do

    return spectrum

end

function calc_for_sigma_coordinate(file_names, rad, pres)
begin

    ; --------------------------------------------------------------------------
    f = addfiles(file_names, "r")
    num_file = dimsizes(file_names)
    num_time = 0
    do i = 0, num_file-1
        num_time = num_time+dimsizes(f[i]->time)
    end do
    ; --------------------------------------------------------------------------
    spectrum = new(get_energy_spectrum_dims(num_time, dimsizes(pres), \
        dimsizes(f[0]->lon)), "float")
    k = 0
    do i = 0, num_file-1
        system("echo [Notice]: Use file "+file_names(i)+".")
        do j = 0, dimsizes(f[i]->time)-1
            ; ==================================================================
            ; interpolate model data onto specific levels
            fake = fspan(0.0, 0.0, dimsizes(f[i]->lev))
            U = vinth2p(f[i]->ua(j,:,:,:), fake, f[i]->lev, pres, \
                f[i]->ps(j,:,:), 1, 1000.0, 1, False)
            V = vinth2p(f[i]->va(j,:,:,:), fake, f[i]->lev, pres, \
                f[i]->ps(j,:,:), 1, 1000.0, 1, False)
            delete(fake)
            ; ==================================================================
            calc_energy_spectrum(rad, U, V, spectrum(k,:,:,:), False)
            k = k+1
        end do
    end do

    return spectrum

end

begin

    ; --------------------------------------------------------------------------
    if (isdefined("file_pattern")) then
        file_names = systemfunc("ls "+file_pattern)
    else
        system("echo [Notice]: Input the patterns of the data file names:")
        file_names = systemfunc("read -p '> ' -e ans && echo $ans")
    end if
    ; --------------------------------------------------------------------------
    rad = 6371e3
    pres = (/200.0/)
    spectrum = calc_for_hybrid_coordinate(file_names, rad, pres)
    ; --------------------------------------------------------------------------
    plot_energy_spectrum("energy_spectrum", spectrum)

end

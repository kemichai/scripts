#
# GMT code to plot seismicity and cross sections
# KM Sep 2017
#

out=quadtree_2.eps

gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set PS_MEDIA A0
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 12
gmt set FONT_LABEL Helvetica
gmt set LABEL_FONT_SIZE 7

#datadir="/Users/home/michaiko/Desktop/map/Data"
#eqdir="/Volumes/GeoPhysics_05/users-data/michaiko/NLL_HDD_2014_2017/input_files"
#topodir="/Users/home/michaiko/Desktop/map/Topo_faults_etc"
#DEMdir="/Users/home/michaiko/Desktop/map/DEMs/south_island/"

#
set -o nounset
set -o errexit


# Define map characteristics
# Define your area
north=-42.5
south=-44.5
east=171.8
west=168.65

# tick='-B0.5/0.5WSen'
# proj='-JM22'
proj='-JM6i'


echo Make basemap ...
# make a basemap
gmt psbasemap -R$west/$east/$south/$north $proj -B0.5wSEn -P -Y12 -K > $out
#gmt makecpt -Ctopo -T-1000/5000/110 -Z -I -A100 > z.cpt
gmt makecpt -Cgray -Z -T0/6000/200 -I > topo.cpt

# grdimage -R -J $topodir/SI_100m_dem_wgs84.grd -Cz.cpt -O -K >> $out

# gmt pscoast -R -Swhite -J  -P -W0.05 -Df  -L169.3/-44.13/-42./50+l -O -K --FONT_ANNOT_PRIMARY=9p >> $out

#pscoast -R -J -Df -W1p -Saliceblue -Gbeige -Lf169.7/-43.9/-43.9/40+l -O -K >> $out

echo Using this clipped grid ....
# gmt grdimage -R -J $DEMdir/clipped_topo.grd -CFrance2.cpt -I$DEMdir/SAMBA_relief.grd  -O -K >> $out

# psxy -R -J $topodir/map.faults_NZ_WGS84.gmt.xy -Wgray0 -W.5p -O -K >> $out

# ************************************************************************************
#Plot events depth by colour 
# ************************************************************************************
# This time, let's make it rainbow colored and call it seis.cpt (0 to 100 with an increament of 20)
# makecpt -Cjet -T0/25/1 -Z -I > seis.cpt

# makecpt -Chot -T10/14/0.5 -Z -A70 > seis.cpt
gmt makecpt -Cviridis -T10/14/0.5 -Z > seis.cpt
#makecpt -Ccool -T10/14/0.5 -Z -I -A50 > seis.cpt


gmt pscoast -W1/0.05 -Df -Swhite -J -R -K -O -L169.3/-44.13/-42./50+l+u >> $out




echo Plotting lakes...
#gmt psxy -R -J $topodir/nz.gmt -W0.05,black -Gwhite -O -K >> $out
gmt psxy boxes_gmt_2area.dat -R -J -W0.25p -O -K -L -Cseis.cpt -t0  >> $out 
echo Plotting faults and eqz ...
#active faults
#gmt psxy -R -J $topodir/activefaults.xy -Wgray5 -W.8p -O -K >> $out

# plots the values of the moments
# awk '{print $3 , $4, $2}' test/values.dat |
# gmt pstext -R -J -O -K  -F+f5p,Helvetica,gray5+jBL+a0 -Gwhite >> $out

# echo Plotting Toponyms as squares...
# gmt psxy -R -J -Ss.1 -W1p -Gblack -O -K  >> $out << END
# 170.56 -43.15
# 170.96 -42.71
# 170.017778 -43.464444
# 170.181944  -43.389167
# 170.814167 -42.895833 #ross
# 170.36 -43.262   # whataroa
# 169.042222 -43.881111 #haast
# END

# 97.7 122.05
awk '{print $1 , $2, $3, $4}' gmt_psxy_bowties.dat |
gmt psxy -R -J -SW0.2i -Gblack -W0.5p -O -K >> $out 



gmt psscale -Dx7/11+o0/0.5i+w1.5i/0.08i+h+e -R -J -Cseis.cpt  -Bx1f1 -By+l"Log@-10@-(@~\123@~@~\115@~@-o@-/@~\101@~)" -O -K --FONT_ANNOT_PRIMARY=7p >> $out
# -t50
# # Carolin
# psxy -R -J -SW0.2i -Gblue -W0.5p -O -K >> $out << END
# 170.02 -43.52 348.8 315.2
# 170.02 -43.52 168.8 135.2
# 170.34 -43.52 345.2 324.8
# 170.34 -43.52 165.2 144.8
# 170.544 -43.396 352.6 319.4
# 170.544 -43.396 172.6 139.39
# 170.583 -43.592 341.1 314.9
# 170.583 -43.592 161.1 134.89
# 170.753 -43.435 346.5 317.5
# 170.753 -43.435 166.5 137.5
# 170.888 -43.218 340 316
# 170.888 -43.218 160 136
# END

# # # Emily
# psxy -R -J -SW0.2i -Gpurple -W0.5p -O -K >> $out << END
# 168.72 -44.92 357.9 325.9
# 168.72 -44.92 177.89 145.89
# 169.49 -44.71 354.2 311.2
# 169.49 -44.71 174.2 131.2
# 168.46 -44.21 351.9 313.9
# 168.46 -44.21 171.8 133.8
# 168.77 -44.04 0.600 322.6
# 168.77 -44.04 180.6 142.6
# 169.52 -44.22 343.3 314.3
# 169.52 -44.22 163.3 134.3
# 168.77 -44.51 340.6 303.6
# 168.77 -44.51 160.60 123.6
# 168.21 -44.51 347.0 320.2
# 168.21 -44.51 167.0 140.2
# 167.98 -44.71 356.1 315.1
# 167.98 -44.71 176.1 135.1
# END


# # John
# psxy -R -J -SW0.2i -Gdarkgreen -W0.5p -O -K >> $out << END
# 169.3	-44.3 350 311
# 169.3	-44.3 170 131
# 170.8	-43.0 0   304
# 170.8	-43.0 180 124
# 171.9	-43.5 327 304
# 171.9	-43.5 147 124
# END

# rose_radius=5
# psrose  -W0.5p -R0/20/0/360 -Gred -S0.4n -A50 -T -Xa8.1 -Ya4 -O -K >> $out << END
# 1 110
# END

# rose_radius=5
# #after
# psrose summfiles/rose_file/SVR.???.after_rose.tmp -R0/$rose_radius/0/360 -Gblue -S0.8n -T -A10   -Xa9.03176 -Ya12.4126 -Wthinnest,0/0/0 -O -K >>$out



# psxy -SV0.03/0.2/0.2 -Wthin -GBlack -O -K -R -J >> $out << END
# 177 -44.2 -105 1.5
# END

#seismicity
# gmt psxy $datadir/locations_4stat.dat -R -J -O -K -h1 -Sci -i6,7,9s0.03  -Wthin,black >> $out
# gmt psxy $datadir/locations_4stat.dat -R -J -O -K -h1 -Sci -i6,7,9s0.03 -Gdarkblue -W0.4p >> $out
# gmt psxy $datadir/nll.dat -R -J -O -K -h1 -Sci -i0,1,3s0.03 -Gdarkblue -W0.4p >> $out
# gmt psxy locations.dat -R -J -O -K -h1 -Sci -i0,1,3s0.04  -Wthin,green >> $out
# gmt psxy $datadir/hypoDD.reloc -R -J -O -K -h1 -Sci -i2,1,16s0.05  -Wthin,darkblue >> $out






echo Plotting stuff...
gmt pstext -R$west/$east/$south/$north -J -O -K  -F+f10p,Helvetica,gray10+jBL+a32  >> $out << END
169.076 -43.876 Alpine Fault
# 170.857 -43.098 Alpine Fault
END
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a0 -Gwhite >> $out << END
# 171.45 -42.7 Hope Fault
171.65 -42.7 HF
END
#Mount Cook
gmt psxy -R -J -Sx.3 -W1p -Gwhite -O -K  >> $out << END
170.1410417 -43.5957472
END
gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jB  >> $out << END
# 170.17 -43.609 Aoraki/Mt Cook
170.065 -43.65 Aoraki 
170.045 -43.69 Mt Cook
170.05 -43.73
END

awk '{print $1, $2, $3}' box_names.dat | gmt pstext -R -J -O -K -F+f6p,Helvetica,gray10+jB -Gwhite >> $out


# gmt pstext -R -J -O -K -F+f4p,Helvetica,gray10+jB -Gwhite >> $out << END
# 169.043785839	-42.9358244297	1
# 170.310091064	-42.3110149502	2
# 169.733970462	-43.4439663678	4
# 170.054487521	-43.2870307202	5
# 170.372363865	-43.1296543917	13
# 170.689641986	-42.9708428898	14
# 170.455690583	-43.3440065372	15
# 170.656247763	-43.3716010368	18
# 170.710654857	-43.4295221738	20
# 170.707206275	-43.1759024544	22
# 170.786278796	-43.1361778096	23
# 170.870381485	-43.3495177534	29
# 170.949581637	-43.3096837743	30
# 171.270828976	-42.8465107036	32
# 170.685611809	-43.5271162003	34
# 170.660481349	-43.6247013613	35
# 170.740181214	-43.5850256111	36
# 170.782359907	-43.6917168216	38
# 170.941682601	-43.612138456	39
# 171.027646051	-43.8261115734	41
# 171.346987473	-43.6655304389	42
# 171.711485694	-43.3083625186	43
# END
# gmt pstext -R -J -O -K -F+f4p,Helvetica,gray10+jB >> $out << END
# 169.357300035	-43.7943940364	3
# 169.949518644	-43.6779153691	6
# 170.136867747	-43.5015839792	7
# 170.296483147	-43.4229124547	8
# 170.244986801	-43.6178745744	9
# 170.404783374	-43.5390575507	10
# 169.789392576	-44.2633148653	11
# 170.43751261	-43.9484361257	12
# 170.614490683	-43.2648680188	16
# 170.761496142	-43.2337983648	24
# 170.840611292	-43.194037392	25
# 170.931882915	-43.1053893556	26
# 170.815887749	-43.2916701836	27
# 170.895045443	-43.2518727661	28
# 171.040877495	-43.220950701	31
# 170.832206721	-43.4964790874	37
# 170.735611111	-43.331919338	19
# 170.790061448	-43.3898039301	21
# 171.125985283	-43.4341192541	40
# 170.605955806	-43.5667551875	33
# 170.564168823	-43.4600060435	17
# END

###############################################################################
# R stuff
gmt makecpt -Cpolar -T0/1/0.05 -Z > pol.cpt
gmt psscale -Dx10.5/0.5+o0/0.5i+w1.5i/0.08i+h+e -R -J -Cpol.cpt  -Bx1f1 -By+l"R values" -O -K --FONT_ANNOT_PRIMARY=7p >> $out
# strike-slip
#awk '{print $1 - 0.05, $2-0.05, $3, $4}' Rvalues_ss.dat |
#gmt psxy -R -J -Wthinner -Cpol.cpt -Ss -i0,1,2,3s0.09 -O -K -V  >> $out
## reverse
#gmt psxy -R -J -Wthinner -Cpol.cpt -Si -i0,1,2,3s0.09 -O -K -V  >> $out << END
#169.95  -43.68 0.8 1
#170.14  -43.5 0.6 1
#170.24  -43.62 0.8 1
#END
###############################################################################



gmt set FONT_ANNOT_PRIMARY 9

# pscoast -R  -J  -P -W0.05 -Df  -L169.3/-44.13/-42./50+l -O -K --FONT_ANNOT_PRIMARY=9p >> $out

gmt pscoast -W1/0.05 -Df  -J -R -K -O -L169.3/-44.13/-42./50+l+u >> $out




echo Plotting velocity of Pacific Plate relative to Australia...
gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jBL+a26 >> $out << END
170.95 -44.15 39.5 mm/yr
END
gmt psxy -N -SV0.15i+e -Gblack -W2p -O -K -R -J >> $out << END
171.2 -44.1 244 1.5
END

#--------------------------------------------------------   
# Inset map of New Zealand showing study area
#--------------------------------------------------------
echo Plotting inset ...
echo ...
region2="-R166/175/-47/-40."
projection2="-JM4.9"
boundaries2="-B2nsew"

gmt psbasemap $region2 $projection2 $boundaries2 -X.01 -Y8.02 -O -K >> $out

# grdimage -R -J $topodir/100m_dem_wgs84.grd -Ctopo.cpt -O -K >> $out
# grdimage -R -J $topodir/SI_100m_dem_wgs84.grd -Ctopo.cpt -O -K >> $out


# pscoast -R -J -Df -W1/0.05 -Swhite  -L176.1/-47/-47/400+l -O -K >> $out
gmt pscoast -R -J -Df -W1/0.05 -Swhite -Gwhite -L174/-46/-42./100+l+u -O -K >> $out

#gmt psxy -R -J $topodir/PB_UTIG_Transform.xy -Sf0.5c/0.03i+l+t -Gblack -W -O -K  >> $out
# psxy -R -J $topodir/Alpine_fault.xy -Sf0.5c/0.03i+l+t -Gred -W -O -K  >> $out

#gmt psxy -R -J $topodir/PB_UTIG_Transform.xy -Sf2c/0.1i+r+s+o1 -Gblack -W -O -K  >> $out
# psxy $topodir/qmap_faults.gmt -R -J -V -O -K -W0.5,black -m >> $out

echo Plot study area in the inset...
# #study area
gmt psxy -R -J -Wthinner,red -O -K  >> $out << END
# 169.5 -44.9
# 172   -43.5
# 171   -42.5
# 168.5 -43.8
# 169.5 -44.9
168.65 -44.2
168.65 -42.5
171.8 -42.5
171.8  -44.2
168.65 -44.2
END

echo Creating legend...
# # construct legend
gmt pslegend <<END -R -J -Dx0.05i/3i+w1.9i/1.1i/TC -C0.1i/0.1i -F+gwhite+pthin -P -O -K >> $out
G -.06i
S .1i t .11i blue 0.01p
G -.09i
S .018i i .11i blue 0.01p 0.2i Boese et al. (2012)
G .02i
S .1i t .11i darkgreen 0.01p
G -.09i
S .018i i .11i darkgreen 0.01p 0.2i Townend et al. (2012)
G .02i
S .1i t .11i purple 0.01p
G -.09i
S .02i i .11i purple 0.01p 0.2i Warren-Smith et al. (2017)
END



gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jBL >> $out << END
172 -46 PAC
168 -42 AUS
END


gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a32  >> $out << END
169.0 -43.5 AF 
END
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a65  >> $out << END
178 -40 HT
END

gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL  >> $out << END
166.5 -46.5 PT
END
# gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jBL+a15 >> $out << END
# 175 -44.2 39.8 mm/yr
# END
# psxy -Sv0.15i+ea -Wthin -GBlack -O -K -R -J >> $out << END
# 177 -44.2 -165 1.5
# END

# Carolin
awk '{print $1 , $2, $3, $4}' CB_gmt_psxy_bowties.dat |
gmt psxy -R -J -SW0.15i -Gred -W0.5p -O -K >> $out


# # Emily
awk '{print $1 , $2, $3, $4}' EWS_gmt_psxy_bowties.dat |
gmt psxy -R -J -SW0.15i -Gpurple -W0.5p -O -K >> $out


# John
awk '{print $1, $2, $3, $4}' Townend_2012_gmt_psxy_bowties.dat |
gmt psxy -R -J -SW0.15i -Gdarkgreen -W0.5p -O -K >> $out





gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
#ps2raster -Tf -A map.ps
evince ${out%.*}.pdf

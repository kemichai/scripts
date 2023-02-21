#
out=station_map.ps
gmtset FORMAT_GEO_MAP D
gmtset FORMAT_GEO_MAP D
gmtset PS_MEDIA A3
gmtset FONT_ANNOT_PRIMARY Helvetica
gmtset FONT_ANNOT_PRIMARY 12
gmtset FONT_LABEL Helvetica
gmtset LABEL_FONT_SIZE 7
#

R=-R168.8/172/-44.2/-42.5r

set -o nounset
set -o errexit

datadir="/data"
topodir="/topo"
# Define map characteristics
# Define your area

north=-42.5
south=-44.2
east=171.8
west=168.65

proj='-JM6i'

echo Creating basemap with or withouth grid...
psbasemap -R$west/$east/$south/$north $proj -B0.5wSEn -P -K > $out

echo Make cpts for topography...
makecpt -Cgray -Z -T0/6000/200 -I > topo.cpt

# grdclip $DEMdir/25m_dem_wgs84.grd -Gclipped_topo.grd -Sb0/NaN
# echo Clipped
# grdgradient clipped_topo.grd -GSAMBA_relief.grd -A270 -Nt0.3
# echo Calculated gradient
# grdview clipped_topo.grd -Qc -JM -R -O -K -Ctopo.cpt -ISAMBA_relief.grd >>tmp.ps
# grdclip SAMBA_relief.grd -Gclipped_bath.grd -Sa0/NaN
# grdgradient clipped_bath.grd -GSAMBA_bath.grd -A270 -Nt0.3
# grdview clipped_bath.grd -Qc -JM -R -O -K -Ctopo.cpt -ISAMBA_bath.grd >> $out
# grdview clipped_topo.grd -Qc -JM -R -O -K -Ctopo.cpt -ISAMBA_relief.grd >> $out
echo Plotted grid

# psbasemap  -Y14 -X2 -P -JOa170/-43/147/4.4i -R$R -K -Ba0 > tmp.ps
# grdview clipped_topo.grd -Qs -J -R -O -K -Ctopo.cpt -ISAMBA_relief.grd 
#  >>tmp.ps
# pscoast -W1/0.05 -Df -J -R -K -O -L169.3/-44.13/-42./50+l >> $out
# psxy faults_NZ_WGS84.gmt -J -R$R -m -O -K -Wthick >> $out
# psxy nz-lakes.gmt -J -R -m -O -K -Wblack -Gwhite >> tmp.ps
# pscoast -W1/thin -Ba0.25swEN -Df -J -R$R -K -O -T$T -L$L >> $out

# grdimage -R -J $topodir/SI_100m_dem_wgs84.grd -Ctopo.cpt -O -K >> $out

# The one to finally use.....!!!!
echo Using this clipped grid ....
grdimage -R -J clipped_topo.grd -Ctopo.cpt -ISAMBA_relief.grd  -O -K >> $out

gmt psscale -Dx4.5/9.5+o0/0.5i+w1.2i/0.08i+h+e -R -J -Ctopo.cpt -Bx2000f1000 -By+l"Topography (m)" -O -K --FONT_ANNOT_PRIMARY=7p >> $out

echo Plotting faults...
psxy -R -J $topodir/activefaults.xy -Wgray5 -W0.8p -O -K  >> $out

echo Plotting lakes...
psxy -R -J $topodir/nz.gmt -W0.05,black -Gwhite -O -K >> $out

pscoast -W1/0.05 -Df -J -R -K -O -L169.3/-44.13/-42./50+l+u >> $out

echo Plotting Alpine and Hope Fault...
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a32  >> $out << END
169.076 -43.876 Alpine Fault
# 170.857 -43.098 Alpine Fault
END
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a0 -Gwhite >> $out << END
# 171.45 -42.7 Hope Fault
171.65 -42.7 HF
END

echo Plotting Lake names...
gmt pstext -R -J -O -K  -F+f8p,Helvetica,navy+jBL+a0 -Gwhite >> $out << END
170.5 -43.9 Lake 
170.46 -43.95 Tekapo
170.12 -44.05 Lake 
170.08 -44.1 Pukaki
END


# gmt pstext -R -J -O -K  -F+f9p,Helvetica,gray10+jBL+a-40  >> $out << END
# 170.9 -43.57 Transect 1
# 170.4 -43.93 Transect 2
# END

echo Plotting seismic stations...
# Plot stations
awk '{print $3, $2}' $datadir/sta_SAMBA.txt |
    psxy -R -J -Si.25 -W.1p -Gred -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_SAMBA_new.txt |
    psxy -R -J -St.25 -W.1p -Gred -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_GEONET.txt |
    psxy -R -J -Si.25 -W.1p -Gdarkorange2 -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_WIZARD.txt |
    psxy -R -J -Si.25 -W.1p -Gmediumpurple -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_DFDP10.txt |
    psxy -R -J -St.25 -W.1p -Gyellow -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_DFDP13.txt |
    psxy -R -J -Si.25 -W.1p -Gdodgerblue -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_ALFA_08.txt |
    psxy -R -J -Si.25 -W.1p -Gcyan -O -K  >> $out
#  awk '{print $3, $2}' $datadir/GPS_sta_GEONET.txt |
#     psxy -R -J -Ss.25 -W.1p -Gdarkorange2 -O -K  >> $out   

# pstext -Wwhite -J -R -O -K -Dj0.25c <<EOF >> $out
# 170.05951 -43.44833 12 0 13 RB COSA
# 170.16933 -43.42643 12 0 13 LT EORO
# 170.00347 -43.51237 12 0 13 RB MTFO
# 170.37158 -43.44123 12 0 13 LT WHYM
# 170.24512 -43.54639 12 0 13 LT LABE
# 170.50255 -43.63881 12 0 13 LT GOVA
# 170.16041 -43.385 12 0 13 RB FRAN
# 170.22333 -43.35161 12 0 13 RB POCR2
# 170.36037 -43.27918 12 0 13 RB WHAT2
# 169.606 -43.908 12 0 13 RB SOLU
# 169.641 -43.688 12 0 13 RB MTBA
# 169.968 -43.613 12 0 13 RB COVA
# 169.936 -43.749 12 0 13 RB LARB
# EOF

#Rivers
gmt pstext -R -J -O -K -F+f8p,Helvetica,gray10+jB -Gwhite >> $out << END
170.16 -43.1 WH
170.37  -43.01 WA
169.68  -43.45 KA
END
# gmt pstext -R -J -O -K  -F+f9p,Helvetica,gray10+jBL -Gwhite >> $out << END
# 171.25 -43.8 T1
# 170.6 -44.1 T2
# END

#Mount Cook
gmt psxy -R -J -Sx.3 -W1p -Gwhite -O -K  >> $out << END
170.1410417 -43.5957472
END
gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jB  >> $out << END
# 170.17 -43.609 Aoraki/Mt Cook
170.065 -43.65 Aoraki/ 
169.985 -43.69 Mount Cook
170.05 -43.73
END

echo Plotting GeoNet seimic site names ...
#GeoNet stations
gmt pstext -R -J -O -K -F+f7p,Helvetica,gray10+jB -Gwhite >> $out << END
170.685 -43.055 WVZ
169.775 -43.51 FOZ
171.014 -43.695 RPZ
# 171.412  -42.764  INZ
168.745  -44.05  JCZ
END

##############################################################################
#Toponyms
##############################################################################
echo Plotting Toponyms labels...
gmt pstext -R -J -O -K -F+f9p,Helvetica,gray9+jB  >> $out << END
# gmt pstext -R -J -O -K -F+f12p,Times-Italic+jLM >> $out << END
170.175 -42.980 Harihari
170.7 -42.7 Hokitika
169.595 -43.3 Fox
169.60 -43.35 Glacier
169.79 -43.19 Franz Josef
169.79  -43.24 Glacier
170.63 -42.87 Ross
170.0 -43.10 Whataroa
168.83 -43.85 Haast
# 170.166667 -44.116667 Lake Pukaki
END

####################################################
#lines
####################################################
echo Plot lines that connect Toponyms to their labels...
#Harihari
gmt psxy -R -J -Wblack -W0.5p -O -K  >> $out << END
170.56 -43.15
170.28 -43.0
END
#fox
gmt psxy -R -J -Wblack -W0.5p -O -K  >> $out << END
170.017778 -43.464444
169.685 -43.371
END
# #Franz
gmt psxy -R -J -Wblack -W0.5p -O -K  >> $out << END
170.181944 -43.389167
169.890 -43.266
END
#whataroa
gmt psxy -R -J -Wblack -W0.5p -O -K  >> $out << END
170.11 -43.115
170.36 -43.262
END

####################################################
echo Plotting Toponyms as squares...
gmt psxy -R -J -Ss.1 -W1p -Gblack -O -K  >> $out << END
170.56 -43.15
170.96 -42.71
170.017778 -43.464444
170.181944  -43.389167
170.814167 -42.895833 #ross
170.36 -43.262   # whataroa
169.042222 -43.881111 #haast
END

gmtset FONT_ANNOT_PRIMARY 9

echo Creating legend...
pslegend <<END -R -J -Dx5.02i/0.06i+w0.93i/1.7i/TC -C0.1i/0.1i -F+gwhite+pthin -P -O -K >> $out
G -.01i
S .04i i .11i red 0.2p 0.18i SAMBA
G .07i
S .04i i .11i mediumpurple 0.2p 0.18i WIZARD
G .07i
S .04i i .11i dodgerblue 0.2p 0.18i DFDP-13
G .07i
S .04i t .11i yellow 0.2p 0.18i DFDP-10
G .07i
S .04i i .11i cyan 0.2p 0.18i ALFA-08
G .07i
S .04i i .11i darkorange2 0.2p 0.18i GeoNet
G .07i
S .04i s .08i black 0.2p 0.18i Towns
G .065i
S .04i - .14i black thick 0.18i Active fault
END


#--------------------------------------------------------   
# Inset map of New Zealand showing study area
#--------------------------------------------------------
echo Plotting inset ...
echo ...
region2="-R165/180/-48/-34."
projection2="-JM4"
boundaries2="-B80nsew"

psbasemap $region2 $projection2 $boundaries2 -X0.0255 -Y6.285 -O -K >> $out

# grdimage -R -J $topodir/100m_dem_wgs84.grd -Ctopo.cpt -O -K >> $out
# grdimage -R -J $topodir/SI_100m_dem_wgs84.grd -Ctopo.cpt -O -K >> $out


# pscoast -R -J -Df -W1/0.05 -Swhite  -L176.1/-47/-47/400+l -O -K >> $out
pscoast -R -J -Df -W0.05p -Swhite -L176.1/-47/-47/400+l -O -K >> $out

psxy -R -J $topodir/PB_UTIG_Transform.xy -Sf0.5c/0.03i+l+t -Gblack -W -O -K  >> $out
# psxy -R -J $topodir/Alpine_fault.xy -Sf0.5c/0.03i+l+t -Gred -W -O -K  >> $out

psxy -R -J $topodir/PB_UTIG_Transform.xy -Sf2c/0.1i+r+s+o1 -Gblack -W -O -K  >> $out
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

gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jBL >> $out << END
176.1 -46.5 PAC
167   -38 AUS
END

gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a32  >> $out << END
169.0 -43.5 AF 
END
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a65  >> $out << END
178 -40 HT
END

gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL  >> $out << END
166 -47.5 PT
END
gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jBL+a15 >> $out << END
175 -44.2 39.8 mm/yr
END
psxy -Sv0.15i+ea -Wthin -GBlack -O -K -R -J >> $out << END
177 -44.2 -165 1.5
END


psxy -R -J -T -O >> $out
ps2raster -Tf -A $out
#ps2raster -Tf -A map.ps
evince ${out%.*}.pdf
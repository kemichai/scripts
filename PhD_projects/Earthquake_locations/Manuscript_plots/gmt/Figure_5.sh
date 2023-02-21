#
out=Relocations_map_view.ps
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set PS_MEDIA A0
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 12
gmt set FONT_LABEL Helvetica
gmt set LABEL_FONT_SIZE 7

R=-R168.8/172/-44.2/-42.5r

set -o nounset
set -o errexit

datadir="/data"
topodir="/topo"

north=-42.5
south=-44.2
east=171.8
west=168.65

proj='-JM6i'

echo Creating basemap with or withouth grid...
gmt psbasemap -R$west/$east/$south/$north $proj -B0.5wSEn -P -K > $out

echo Make cpts for topography and seismicity ...
gmt makecpt -Cgray -Z -T0/6000/200 -I > topo.cpt
gmt makecpt -Cjet -T0/25/1 -Z -I > seis.cpt

echo Using this clipped grid ....
#grdimage -R -J clipped_topo.grd -Ctopo.cpt -ISAMBA_relief.grd  -O -K >> $out

echo Plotting faults...
#gmt psxy -R -J $topodir/activefaults.xy -Wgray5 -W0.8p -O -K  >> $out

echo Plotting lakes...
#gmt psxy -R -J $topodir/nz.gmt -W0.05,black -Gwhite -O -K >> $out

gmt pscoast -W1/0.05 -Df -J -R -K -O  -L169.3/-44.13/-42./50+l+u >> $out


echo Plotting Alpine Fault...
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a32  >> $out << END
169.076 -43.876 Alpine Fault
END
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a0 -Gwhite >> $out << END
# 171.45 -42.7 Hope Fault
171.65 -42.7 HF
END

echo Plot earthquake epicenters...
awk '{if ($17<=4) print $3, $2, $4, 1+$17}' hypoDD.reloc3 | gmt psxy -i0,1,2,3s0.018 -Sc -R -J \
-O -K  -W.25 -Cseis.cpt >> $out
awk '{if ($17>4) print $3, $2, $4, 1+$17}' hypoDD.reloc3 | gmt psxy -i0,1,2,3s0.018 -Sc -R -J \
-O -K  -W.25 -Cseis.cpt >> $out

# Magnitude scale
# 1 to 4 stands for -1 to 4
gmt psxy  -i0,1,2,3s0.018 -Sc -R -J -O -K -W.25 -Cseis.cpt >> $out << END
169.57 -43.30  9.0 0.5 
169.57 -43.325  9.0 1.0 
169.57 -43.35  9.0 2.0 
169.57 -43.375  9.0 3.0
169.57 -43.40  9.0 4.0
169.57 -43.425  9.0 5.0
END
gmt pstext -R -J -O -K  -F+f4p,Helvetica,gray10+jBL+a0 -Gwhite >> $out << END
169.60 -43.308  M<0
169.60 -43.333  M=0
169.60 -43.358  M=1
169.60 -43.383  M=2
169.60 -43.408  M=3
169.60 -43.433  M=4
END

gmt psscale -Dx5/9+o0/0.5i+w1.5i/0.08i+h+e -R -J  -Cseis.cpt -Bx5f5 -By+l"Depth (km)" -O -K --FONT_ANNOT_PRIMARY=7p >> $out

echo Plotting seismic stations...
awk '{print $3, $2}' $datadir/sta_SAMBA.txt |
    gmt psxy -R -J -Si.25 -W.1p -Gred -O -K  >> $out
 awk '{print $3, $2}' $datadir/sta_SAMBA_new.txt |
    gmt psxy -R -J -St.25 -W.1p -Gred -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_GEONET.txt |
    gmt psxy -R -J -Si.25 -W.1p -Gdarkorange2 -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_WIZARD.txt |
    gmt psxy -R -J -Si.25 -W.1p -Gmediumpurple -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_DFDP10.txt |
    gmt psxy -R -J -St.25 -W.1p -Gyellow -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_DFDP13.txt |
    gmt psxy -R -J -Si.25 -W.1p -Gdodgerblue -O -K  >> $out
awk '{print $3, $2}' $datadir/sta_ALFA_08.txt |
    gmt psxy -R -J -Si.25 -W.1p -Gcyan -O -K  >> $out

#Rivers
# gmt pstext -R -J -O -K -F+f8p,Helvetica,gray10+jB  >> $out << END
# 170.16 -43.1 WH
# 170.37  -43.01 WA
# 169.68  -43.45 KA
# END
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
170.065 -43.65 Aoraki 
170.045 -43.69 Mt Cook
170.05 -43.73
END


echo Plotting GeoNet seimic site names ...
#GeoNet stations
# gmt pstext -R -J -O -K -F+f6p,Helvetica,gray10+jB -Gwhite >> $out << END
# 170.685 -43.055 WVZ
# 169.775 -43.51 FOZ
# 171.014 -43.695 RPZ
# # 171.412  -42.764  INZ
# 168.745  -44.05  JCZ
# END

echo Plotting cross section lines parallel to AF ...
start_lon_par='169.'
start_lat_par='-44'
end_lon_par='171.4'
end_lat_par='-42.8'

gmt psxy << END -R -J -O -W0.5,black,- -K>> $out
$start_lon_par $start_lat_par
$end_lon_par $end_lat_par
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$start_lon_par $start_lat_par F
END
gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$end_lon_par $end_lat_par F'
END
# #..................'

start_lon_par2='169.3'
start_lat_par2='-44.1'
end_lon_par2='171.7'
end_lat_par2='-42.9'

gmt psxy << END -R -J -O -W0.5,black,- -K>> $out
$start_lon_par2 $start_lat_par2
$end_lon_par2 $end_lat_par2
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$start_lon_par2 $start_lat_par2 G
END
gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$end_lon_par2 $end_lat_par2 G'
END

# #..................
echo Plotting cross section lines peprendicular to AF ...

start_lon_per='170.6'
start_lat_per='-42.95'
end_lon_per='171.27'
end_lat_per='-43.65'

gmt psxy << END -R -J -O -W0.5,black,- -K>> $out
$start_lon_per $start_lat_per
$end_lon_per $end_lat_per
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$start_lon_per $start_lat_per A
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$end_lon_per $end_lat_per A'
END

# #..................''
start_lon_per1='170.4'
start_lat_per1='-43.05'
end_lon_per1='171.07'
end_lat_per1='-43.75'


gmt psxy << END -R -J -O -W0.5,black,- -K>> $out
$start_lon_per1 $start_lat_per1
$end_lon_per1 $end_lat_per1
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$start_lon_per1 $start_lat_per1 B
END
gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$end_lon_per1 $end_lat_per1 B'
END

# #..................'
start_lon_per2='170.2'
start_lat_per2='-43.15'
end_lon_per2='170.87'
end_lat_per2='-43.85'

gmt psxy << END -R -J -O -W0.5,black,- -K>> $out
$start_lon_per2 $start_lat_per2
$end_lon_per2 $end_lat_per2
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$start_lon_per2 $start_lat_per2 C
END
gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$end_lon_per2 $end_lat_per2 C'
END

# #..................'
start_lon_per3='170'
start_lat_per3='-43.25'
end_lon_per3='170.67'
end_lat_per3='-43.95'

gmt psxy << END -R -J -O -W0.5,black,- -K>> $out
$start_lon_per3 $start_lat_per3
$end_lon_per3 $end_lat_per3
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$start_lon_per3 $start_lat_per3 D
END
gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
$end_lon_per3 $end_lat_per3 D'
END


# #..................'
start_lon_per4='169.8'
start_lat_per4='-43.35'
end_lon_per4='170.47'
end_lat_per4='-44.05'

gmt psxy << END -R -J -O -W0.5,black,- -K>> $out
$start_lon_per4 $start_lat_per4
$end_lon_per4 $end_lat_per4
END

gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$start_lon_per4 $start_lat_per4 E
END
gmt pstext -R -J -D0/0.23 -O -K -F+f9p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$end_lon_per4 $end_lat_per4 E'
END


gmt set FONT_ANNOT_PRIMARY 9

echo Creating legend...
gmt pslegend <<END -R -J -Dx5.07i/0.06i+w0.88i/1.7i/TC -C0.1i/0.1i -F+gwhite+p -P -O -K >> $out
G -.01i
S .05i i .12i red 0.2p .2i SAMBA
G .07i
S .05i i .12i mediumpurple 0.2p .2i WIZARD
G .07i
S .05i i .12i dodgerblue 0.2p .2i DFDP-13
G .07i
S .05i t .12i yellow 0.2p .2i DFDP-10
G .07i
S .05i i .12i cyan 0.2p .2i ALFA-08
G .07i
S .05i i .12i darkorange2 0.2p .2i GeoNet
G .07i
# S .05i i .08i darkgreen 0.2p .2i SIGHT
# G .07i
S .05i s .08i black 0.2p .2i Towns
G .07i
S .05i - .15i black thick 0.2i Active fault
END




#--------------------------------------------------------
# Inset map of New Zealand showing study area
#--------------------------------------------------------

# dir="/home/konstantinos/Desktop/gmt"
# region2="-R165/180/-48/-34."
# projection2="-JM4"
# boundaries2="-B80nsew"

# psbasemap $region2 $projection2 $boundaries2 -X0.01 -Y6.32 -O -K >> $out

# # grdimage -R -J $dir/100m_dem_wgs84.grd -Cz.cpt -O -K >> $out
# # grdimage -R -J $dir/SI_100m_dem_wgs84.grd -Cz.cpt -O -K >> $out

# pscoast -R -J -Df -W0.1p -Swhite -L176.1/-47/-47/400+l -O -K >> $out

# # psxy -R -J $dir/activefaults.xy -Wblack -W0.2p -O -K  >> $out
# # psxy -R -J $dir/PB_UTIG_Transform.xy -Sf0.5c/0.03i+l+t -Gblack -W -O -K  >> $out
# # psxy -R -J $dir/PB_UTIG_Transform.xy -Sf2c/0.1i+r+s+o1 -Gblack -W -O -K  >> $out




# gmt pstext -R -J -O -K -F+f7p,Helvetica,gray10+jBL >> $out << END
# 176.1 -46.5 PA
# 167   -38 AU
# END


# gmt pstext -R -J -O -K  -F+f7p,Helvetica,gray10+jBL+a32  >> $out << END
# 169.0 -43.5 Alpine Fault
# END
# gmt pstext -R -J -O -K  -F+f7p,Helvetica,gray10+jBL+a65  >> $out << END
# 178 -40 Hikurangi Trough
# END

# gmt pstext -R -J -O -K  -F+f7p,Helvetica,gray10+jBL  >> $out << END
# 166 -46 Puysegur Trench
# END
# gmt pstext -R -J -O -K -F+f7p,Helvetica,gray10+jBL+a15 >> $out << END
# 175.4 -44.2 40 mm/yr
# END
# psxy -SV0.03/0.2/0.2 -Wthin -GBlack -O -K -R -J >> $out << END
# 177 -44.2 -105 1.5
# END


# echo Plot study area in the inset...
# # #study area
# gmt psxy -R -J -Wthinner,red -O -K  >> $out << END
# # 169.5 -44.9
# # 172   -43.5
# # 171   -42.5
# # 168.5 -43.8
# # 169.5 -44.9
# 168.65 -44.2
# 168.65 -42.5
# 171.8 -42.5
# 171.8  -44.2
# 168.65 -44.2
# END


gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
#ps2raster -Tf -A map.ps
evince ${out%.*}.pdf

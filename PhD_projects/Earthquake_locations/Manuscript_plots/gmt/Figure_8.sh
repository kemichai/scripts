#
# GMT code to plot quadtree boxplots
# KM Sep 2017
#

out=depth_variations.eps

gmtset FORMAT_GEO_MAP D
gmtset FORMAT_GEO_MAP D
gmtset PS_MEDIA A0
gmtset FONT_ANNOT_PRIMARY Helvetica
gmtset FONT_ANNOT_PRIMARY 12
gmtset FONT_LABEL Helvetica
gmtset LABEL_FONT_SIZE 7

datadir="/data"
topodir="/topo"

set -o nounset
set -o errexit


# Define map characteristics
# Define your area
north=-42.5
south=-44.2
east=171.8
west=168.65

# tick='-B0.5/0.5WSen'
# proj='-JM22'
proj='-JM6i'


echo Make basemap ...
psbasemap -R$west/$east/$south/$north $proj -B0.5wSEn -P -Y12 -K > $out

echo Make cpts for topography and seismicity ...
makecpt -Cgray -Z -T0/6000/200 -I > topo.cpt
makecpt -Cjet -T0/25/1 -Z -I > seis.cpt


echo Using this clipped grid ....
grdimage -R -J $topo/clipped_topo.grd -Ctopo.cpt -ISAMBA_relief.grd  -O -K >> $out

echo Create scale...
psscale -DJBC+o0/0.4i+w4i/0.15i+h -R -J -Cseis.cpt -Bx5f5 -By+l"Depth (km)" -O -K >> $out

echo Plot lakes...
psxy -R -J $topodir/nz.gmt -W0.05,black -Gwhite -O -K >> $out


echo Plot earthquake epicenters and the boxes on top...
awk '{print $3, $2, $4, 2+$17}' hypoDD.reloc3 | psxy -i0,1,2,3s0.012 -Scc -R -J \
-O -K  -W.25 -Cseis.cpt -H15 >> $out
# Quadtree boxes
psxy 95per_box_.dat -R -J -W0.25p -O -K -L -Cseis.cpt -t0  >> $out 

echo Plotting various stuff...
gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a32  >> $out << END
169.076 -43.876 Alpine Fault
# 170.857 -43.098 Alpine Fault
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

gmtset FONT_ANNOT_PRIMARY 9

pscoast -W1/0.05 -Df -J -R -K -O -L169.3/-44.13/-42./50+l+u >> $out

echo Plotting faults and eqz ...
psxy -R -J $topodir/activefaults.xy -Wgray5 -W.8p -O -K >> $out


echo Plotting seismic sites...
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

echo Creating legend ...
pslegend <<END -R -J -Dx5.07i/0.06i+w0.88i/1.7i/TC -C0.1i/0.1i -F+gwhite+p -P -O -K >> $out
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
S .05i s .08i black 0.2p .2i Towns
G .07i
S .05i - .15i black thick 0.2i Active fault
END


#--------------------------------------------------------   
# Inset map of New Zealand showing study area
#--------------------------------------------------------

# region2="-R165/180/-48/-34."
# projection2="-JM4"
# boundaries2="-B80nsew"

# echo Plotting inset ...
# echo ...
# psbasemap $region2 $projection2 $boundaries2 -X0.01 -Y6.32 -O -K >> $out

# # grdimage -R -J $topodir/100m_dem_wgs84.grd -Cz.cpt -O -K >> $out
# # grdimage -R -J $topodir/SI_100m_dem_wgs84.grd -Cz.cpt -O -K >> $out

# pscoast -R -J -Df -W0.1p -Swhite -L176.1/-47/-47/400+l -O -K >> $out

# # psxy -R -J $dir/activefaults.xy -Wblack -W0.2p -O -K  >> $out
# psxy -R -J $topodir/PB_UTIG_Transform.xy -Sf0.5c/0.03i+l+t -Gblack -W -O -K  >> $out
# # psxy -R -J $topodir/Alpine_fault.xy -Sf0.5c/0.03i+l+t -Gred -W -O -K  >> $out

# psxy -R -J $topodir/PB_UTIG_Transform.xy -Sf2c/0.1i+r+s+o1 -Gblack -W -O -K  >> $out
# # psxy -R -J $topodir/Alpine_fault.xy -Sf2c/0.1i+r+s+o1 -Gred -W -O -K  >> $out




# gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jBL >> $out << END
# 176.1 -46.5 PA
# 167   -38 AU
# END


# gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a32  >> $out << END
# 169.0 -43.5 Alpine Fault
# END
# gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL+a65  >> $out << END
# 178 -40 Hikurangi Trough
# END

# gmt pstext -R -J -O -K  -F+f10p,Helvetica,gray10+jBL  >> $out << END
# 166 -46 Puysegur Trench
# END
# gmt pstext -R -J -O -K -F+f10p,Helvetica,gray10+jBL+a15 >> $out << END
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


psxy -R -J -T -O >> $out
ps2raster -Tf -A $out
#ps2raster -Tf -A map.ps
evince ${out%.*}.pdf

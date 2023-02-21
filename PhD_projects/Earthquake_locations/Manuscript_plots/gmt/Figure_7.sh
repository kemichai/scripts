out=parallel.eps




datadir="/Users/home/michaiko/Desktop/map/Data"
topodir="/Users/home/michaiko/Desktop/map/Topo_faults_etc"
growcldir="/Users/home/michaiko/Desktop/Growclust/"

#

gmtset FORMAT_GEO_MAP D
gmtset PS_MEDIA A0
gmtset FONT_ANNOT_PRIMARY Helvetica
gmtset FONT_ANNOT_PRIMARY 25
gmtset FONT_LABEL Helvetica
gmtset LABEL_FONT_SIZE 35

gmt psbasemap -R0/6.5/0/7.5 -JX45/50 -P -B -K > $out


width='10'
start_lon='169.3'
start_lat='-44.1'
end_lon='171.7'
end_lat='-42.9'


awk '{ print $3, $2, $4 ,$17}' hypoDD.reloc3 | project -C$start_lon/$start_lat -E$end_lon/$end_lat\
 -W-$width/$width -Q -Fpz > projection1.dat
# awk '{if ($14>10) print($9, $8, $10, $11)}' $growcldir/growclust_cat.dat | project -C$start_lon/$start_lat -E$end_lon/$end_lat\
#  -W-$width/$width -Q -Fpz > projection1.dat
awk '{print($9, $8, $10, $11)}' $growcldir/growclust_cat.dat | project -C$start_lon/$start_lat -E$end_lon/$end_lat\
 -W-$width/$width -Q -Fpz > projection1_.dat

awk '{print($1, $2, -2)}' Aoraki.dat | project -C$start_lon/$start_lat -E$end_lon/$end_lat \
-W-$width/$width -Q -Fpz > Aoraki_g.dat
# LFEs
awk '{print($1, $2, $3)}' LFE_LMB.txt | project -C$start_lon/$start_lat -E$end_lon/$end_lat\
 -W-$width/$width -Q -Fpz > LFE_proj.dat


# Create a text file that contains the meeting points with the other cross sections
awk '{print($1, $2, -2)}' Cross_section_A.dat | project -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-1/1 -Q -Fpz > Cross_A.dat
awk '{print($1, $2, -2)}' Cross_section_B.dat | project -C$start_lon/$start_lat\
 -E$end_lon/$end_lat -W-1/1 -Q -Fpz > Cross_B.dat
awk '{print($1, $2, -2)}' Cross_section_C.dat | project -C$start_lon/$start_lat\
 -E$end_lon/$end_lat -W-1/1 -Q -Fpz > Cross_C.dat
awk '{print($1, $2, -2)}' Cross_section_D.dat | project -C$start_lon/$start_lat\
 -E$end_lon/$end_lat -W-1/1 -Q -Fpz > Cross_D.dat
awk '{print($1, $2, -2)}' Cross_section_E.dat | project -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-1/1 -Q -Fpz > Cross_E.dat




##############################################################################
# Create a file that contains the cross section points in Lat/Lon
awk '{ print $3, $2, $4 ,$17}' hypoDD.reloc | project -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-1/1 -Q -Fxyz > Cross_section_G.dat
##############################################################################

awk '{print $3, $2, "-1" , $1,$5}' sta_ALFA_08.txt |  project -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-$width/$width -Q -Fpz > station_proj_a_1.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_old.txt |  project  -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-$width/$width -Q -Fpz > station_proj_s_1.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_new.txt |  project  -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-$width/$width -Q -Fpz > station_proj_sn_1.dat
awk '{print $3, $2, "-1" , $1}' sta_WIZARD.txt |  project  -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-$width/$width -Q -Fpz > station_proj_w_1.dat
awk '{print $3, $2, "-1" , $1}' sta_GEONET.txt |  project  -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-$width/$width -Q -Fpz > station_proj_g_1.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP10.txt |  project  -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-$width/$width -Q -Fpz > station_proj_dfdpa_1.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP13.txt |  project  -C$start_lon/$start_lat \
-E$end_lon/$end_lat -W-$width/$width -Q -Fpz > station_proj_dfdpb_1.dat



gmtset FONT_TITLE 12p,Helvetica,black
# awk '{print($1,$2,$3)}' projection1_.dat | psxy -Sci -i0,1,2s0.045 -W.25 -Gdimgrey -R0/235/-3/25 -JX40/-8 \
# -Bxafg1000+l"Distance (km)" -By5+l"Depth (km)" -Y4 -X1 -BwSnE -O -K >> $out
# awk '{print($1,$2,$3)}' projection1_.dat | psxy -Sci -i0,1,2s0.045 -W.25 -Gdimgrey -R0/235/-3/25 -JX40/-8 \
# -Bxafg1000+l"Distance (km)" -By5+l"Depth (km)" -BwSnE -O -K >> $out

awk '{print($1,$2,$3)}' projection1.dat | psxy -Sci -i0,1,2s0.045 -W.25 -Gdimgrey -R0/235/-3/25 -JX40/-8 \
-Bxafg1000+l"Distance (km)" -By5+l"Depth (km)" -X1 -Y1 -BwSnE -O -K >> $out
# awk '{print($1,$2,$3)}' projection1.dat | psxy -Sci -i0,1,2s0.045 -W.25 -Gred -R0/235/-3/25 -JX40/-8 \
# -Bxafg1000+l"Distance (km)" -By5+l"Depth (km)" -BwSnE -O -K >> $out



awk '{print($1,$2,$2)}' LFE_proj.dat | psxy -Sa0.4  -W.25 -Cseis.cpt -R -J -O -K -V >> $out


# plot seismic networks
psxy station_proj_a_1.dat -Si.4 -W.5 -Gcyan -R -J -O -K >> $out
psxy station_proj_s_1.dat -Si.4 -W.5 -Gred -R -J  -O -K >> $out
psxy station_proj_sn_1.dat -St.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_w_1.dat -Si.4 -W.5 -Gmediumpurple -R -J -O -K >> $out
psxy station_proj_g_1.dat -Si.4 -W.5 -Gorange -R -J -O -K >> $out
psxy station_proj_dfdpa_1.dat -Si.4 -W.5 -Gyellow -R -J -O -K >> $out
psxy station_proj_dfdpb_1.dat -Si.4 -W.5 -Gdodgerblue -R -J -O -K >> $out
awk '{print($1,$2,$2)}' Cross_A.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_B.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J  -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_C.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J  -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_D.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J  -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_E.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J  -O -K -V >> $out
awk '{print($1,$2,$2)}' Aoraki_g.dat| psxy -Sd0.5  -W.5 -Gblack -R \
-J  -O -K -V >> $out

gmt pstext -R -JX -O -K -F+f18p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
3.5 0 G
226 0 G'
END
#'

# gmt grdcut $topodir/clipped_topo.grd -R168.65E/171.8E/44.2S/42.5S -Gspac_33.nc

# cat << EOF > ridge.txt
# 169.25 -44.05
# 169.343 -44.15
# EOF

# gmt grdtrack ridge.txt -G@spac_33.nc -C470k/1k/5k+v -Sm+sstack.txt > table.txt
# # gmt psxy -R -J -O -K -W0.5p table.txt >> $out
# # # Show upper/lower values encountered as an envelope
# gmt convert stack.txt -o0,5 > env.txt
# gmt convert stack.txt -o0,6 -I -T >> env.txt

# gmt psxy -R0/235/0/4000 -Bxafg1000+l"Distance from ridge (km)" -Byaf+l"Depth (m)" -BwSnE \
#     -JX40/-3 -O -K -Glightgray env.txt -Y8 >> $out
# gmt psxy -R -J -O -K -W1.5p stack.txt -Byaf+l"Depth (m)"  >> $out



############################################################################################
# F to F' cross section
############################################################################################
start_lon_par='169.'
start_lat_par='-44'
end_lon_par='171.4'
end_lat_par='-42.8'
width='10'

##############################################################################
# Create a file that contains the cross section points in Lat/Lon
awk '{ print $3, $2, $4 ,$17}' hypoDD.reloc3 | project -C$start_lon_par/$start_lat_par \
-E$end_lon_par/$end_lat_par -W-1/1 -Q -Fxyz > Cross_section_F.dat
##############################################################################

awk '{ print $3, $2, $4 ,$17}' hypoDD.reloc3 | project -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par \
-W-$width/$width -Q -Fpz > projection2.dat

awk '{print($1, $2, $3)}' LFE_LMB.txt | project -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > LFE_proj2.dat
width_par='50'

awk '{print($1, $2, -2)}' Aoraki.dat | project -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par \
-W-$width_par/$width_par -Q -Fpz > Aoraki_f.dat
# Wallace's locking depths
awk '{print($1, $2, $3)}' Wallace.txt | project -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par\
 -W-50/50 -Q -Fpz > GEOD_proj.dat
awk '{print($1, $2, $3)}' Lamb.txt | project -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par\
 -W-50/50 -Q -Fpz > GEOD2_proj.dat

awk '{print $3, $2, "-1" , $1}' sta_ALFA_08.txt |  project -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > station_proj_a_1a.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_old.txt |  project  -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > station_proj_s_1a.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_new.txt |  project  -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > station_proj_sn_1a.dat
awk '{print $3, $2, "-1" , $1}' sta_WIZARD.txt |  project  -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > station_proj_w_1a.dat
awk '{print $3, $2, "-1" , $1}' sta_GEONET.txt |  project  -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > station_proj_g_1a.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP10.txt |  project  -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > station_proj_dfdpa_1a.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP13.txt |  project  -C$start_lon_par/$start_lat_par -E$end_lon_par/$end_lat_par -W-$width/$width -Q -Fpz > station_proj_dfdpb_1a.dat

# Create a text file that contains the meeting points with the other cross sections
awk '{print($1, $2, -2)}' Cross_section_A.dat | project -C$start_lon_par/$start_lat_par \
-E$end_lon_par/$end_lat_par -W-0.5/0.5 -Q -Fpz > Cross_A2.dat
awk '{print($1, $2, -2)}' Cross_section_B.dat | project -C$start_lon_par/$start_lat_par \
-E$end_lon_par/$end_lat_par -W-0.5/0.5 -Q -Fpz > Cross_B2.dat
awk '{print($1, $2, -2)}' Cross_section_C.dat | project -C$start_lon_par/$start_lat_par \
-E$end_lon_par/$end_lat_par -W-0.5/0.5 -Q -Fpz > Cross_C2.dat
awk '{print($1, $2, -2)}' Cross_section_D.dat | project -C$start_lon_par/$start_lat_par \
-E$end_lon_par/$end_lat_par -W-0.5/0.5 -Q -Fpz > Cross_D2.dat
awk '{print($1, $2, -2)}' Cross_section_E.dat | project -C$start_lon_par/$start_lat_par \
-E$end_lon_par/$end_lat_par -W-0.5/0.5 -Q -Fpz > Cross_E2.dat





gmtset FONT_TITLE 12p,Helvetica,black



awk '{print($1,$2,$3)}' projection2.dat | psxy -Sci -i0,1,2s0.045 -W.25 -Gdimgrey\
 -R0/235/-3/25 -JX40/-8 -Bxafg1000 -By5+l"Depth (km)" -Y9 -BwSnE -O -K >> $out
awk '{print($1,$2,$2)}' LFE_proj2.dat| psxy -Sa0.4  -W.25 -Cseis.cpt -R0/235/-3/25 \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' GEOD_proj.dat | psxy -Ss0.2  -W.25 -Ggreen -R -J -O -K -V >> $out
awk '{print($1,$2,$2)}' GEOD2_proj.dat | psxy -Ss0.5  -W.25 -Gmediumpurple -R -J -O -K -V >> $out

# stations
psxy station_proj_a_1a.dat -Si.4 -W.5 -Gcyan -R -J -O -K >> $out
psxy station_proj_s_1a.dat -Si.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_sn_1a.dat -St.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_w_1a.dat -Si.4 -W.5 -Gmediumpurple -R -J -O -K >> $out
psxy station_proj_g_1a.dat -Si.4 -W.5 -Gorange -R -J -O -K >> $out
psxy station_proj_dfdpa_1a.dat -Si.4 -W.5 -Gyellow -R -J -O -K >> $out
psxy station_proj_dfdpb_1a.dat -Si.4 -W.5 -Gdodgerblue -R -J -O -K >> $out
awk '{print($1,$2,$2)}' Aoraki_f.dat| psxy -Sd0.5  -W.5 -Gblack -R \
-J  -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_A2.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_B2.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_C2.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_D2.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_E2.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out

gmt pstext -R -JX -O -K -F+f18p,Helvetica,gray10+jB  -TO -Gwhite -W0.1 >> $out << END
3.5 0 F
226 0 F'
END

awk '{print($1,$2,1)}' AGES_F.dat | psxy -Ssi -i0,1,2s0.09 -W.25 -Gblack\
 -R0/235/0/80 -JX40/-4 -Bxafg1000 -By20+l"Ar/Ar K-Ar age data for hornblendes" -Y9 -BwSnE -O -K >> $out






psxy -R -J -T -O >> $out
ps2raster -Tf -A $out
#ps2raster -Tf -A map.ps
evince ${out%.*}.pdf

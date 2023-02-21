out=Mult_cross.eps

gmtset FORMAT_GEO_MAP D
gmtset PS_MEDIA A0
gmtset FONT_ANNOT_PRIMARY Helvetica
gmtset FONT_ANNOT_PRIMARY 25
gmtset FONT_LABEL Helvetica
gmtset LABEL_FONT_SIZE 35


datadir="../data"


gmt psbasemap -R0/6.5/0/7.5 -JX50/60 -P -B -K > $out

start_lona='170.6'
start_lata='-42.95'
end_lona='171.27'
end_lata='-43.65'
width='10'


awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_lona/$start_lata \
-E$end_lona/$end_lata -W-$width/$width -Q -Fpz > projection_a.dat
awk '{print($1, $2, -2)}' Aoraki.dat | project -C$start_lona/$start_lata -E$end_lona/$end_lata \
-W-$width/$width -Q -Fpz > Aoraki_a.dat
awk '{print($1, $2, $3)}' LFE_LMB.txt | project -C$start_lona/$start_lata -E$end_lona/$end_lata\
 -W-$width/$width -Q -Fpz > LFE_proj_a.dat


awk '{print($1, $2, -2)}' Cross_section_G.dat | project -C$start_lona/$start_lata\
 -E$end_lona/$end_lata -W-1/1 -Q -Fpz > Cross_G.dat
awk '{print($1, $2, -2)}' Cross_section_F.dat | project -C$start_lona/$start_lata\
-E$end_lona/$end_lata -W-1/1 -Q -Fpz > Cross_F.dat


##############################################################################
# Create a file that contains the cross section points in Lat/Lon
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_lona/$start_lata \
-E$end_lona/$end_lata -W-1/1 -Q -Fxyz > Cross_section_A.dat
##############################################################################


awk '{print $3, $2, "-1" , $1}' $datadir/sta_ALFA_08.txt |  project -C$start_lona/$start_lata \
-E$end_lona/$end_lata -W-$width/$width -Q -Fpz > station_proj_a_4.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_SAMBA.txt |  project  -C$start_lona/$start_lata\
 -E$end_lona/$end_lata -W-$width/$width -Q -Fpz > station_proj_sn_4.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_SAMBA_new.txt |  project  -C$start_lona/$start_lata \
-E$end_lona/$end_lata -W-$width/$width -Q -Fpz > station_proj_so_4.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_WIZARD.txt |  project  -C$start_lona/$start_lata \
-E$end_lona/$end_lata -W-$width/$width -Q -Fpz > station_proj_w_4.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_GEONET.txt |  project  -C$start_lona/$start_lata\
 -E$end_lona/$end_lata -W-$width/$width -Q -Fpz > station_proj_g_4.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_DFDP10.txt |  project  -C$start_lona/$start_lata \
-E$end_lona/$end_lata -W-$width/$width -Q -Fpz > station_proj_dfdpa_4.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_DFDP13.txt |  project  -C$start_lona/$start_lata \
-E$end_lona/$end_lata -W-$width/$width -Q -Fpz > station_proj_dfdpb_4.dat

gmtset FONT_TITLE 12p,Helvetica,black

awk '{print($1,$2,$3)}' projection_a.dat | psxy -Sci -i0,1,2s0.048 -W.25 -Gdimgrey \
-R0/100/-3/25 -JX17/-8 -Bx10 -By5+l"Depth (km)" -BwSnE -Y50 -X1 -O -K >> $out
awk '{print($1,$2,$2)}' LFE_proj_a.dat | psxy -Sa0.4  -W.25 -Cseis.cpt -R -J -O -K -V >> $out

# awk '{print($1,-$2/1000)}' profile.dat | psxy -W.25 -R0/100/-3/25 $new_projection $new_boundaries  -K -O >> $out
awk '{print($1,$2,$2)}' Cross_F.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_G.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out

psxy station_proj_a_4.dat -Si.4 -W.5 -Gcyan -R -J -O -K >> $out
psxy station_proj_so_4.dat -Si.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_sn_4.dat -St.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_w_4.dat -Si.4 -W.5 -Gmediumpurple -R -J -O -K >> $out
psxy station_proj_g_4.dat -Si.4 -W.5 -Gorange -R -J -O -K >> $out
psxy station_proj_dfdpa_4.dat -Si.4 -W.5 -Gyellow -R -J -O -K >> $out
psxy station_proj_dfdpb_4.dat -Si.4 -W.5 -Gdodgerblue -R -J -O -K >> $out

gmt pstext -R -J -O -K -F+f18p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
2 0 A
95.5 0 A'
END

# #..................
start_lone='170.4'
start_late='-43.05'
end_lone='171.07'
end_late='-43.75'
width='10'

##############################################################################
# Create a file that contains the cross section points in Lat/Lon
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_lone/$start_late \
-E$end_lone/$end_late -W-1/1 -Q -Fxyz > Cross_section_B.dat
##############################################################################'


# awk '{if ($4<25)print($2, $1, $3,$4+1)}' $datadir/loc_mag.dat | project -C$start_lone/$start_late -E$end_lone/$end_late -W-$width/$width -Q -Fpz > projection5.dat
# awk '{if ($4>6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat |  project -C$start_lone/$start_late -E$end_lone/$end_late -W-$width/$width -Q -Fpz > projection5.dat
# awk '{if ($4<6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat |  project -C$start_lone/$start_late -E$end_lone/$end_late -W-$width/$width -Q -Fpz > projection5_.dat
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_lone/$start_late \
-E$end_lone/$end_late -W-$width/$width -Q -Fpz > projection5.dat

awk '{print($1, $2, -2)}' Cross_section_G.dat | project -C$start_lone/$start_late \
-E$end_lone/$end_late -W-1/1 -Q -Fpz > Cross_G_.dat
awk '{print($1, $2, -2)}' Cross_section_F.dat | project -C$start_lone/$start_late \
-E$end_lone/$end_late -W-1/1 -Q -Fpz > Cross_F_.dat


# awk '{ print $3, $2, $4 ,$17}' $datadir/../hypoDD.reloc3 | project -C$start_lone/$start_late -E$end_lone/$end_late -W-$width/$width -Q -Fpz > projection5.dat



awk '{print $3, $2, "-1" , $1}' $datadir/sta_ALFA_08.txt |  project -C$start_lone/$start_late\
 -E$end_lone/$end_late -W-$width/$width -Q -Fpz > station_proj_a_5.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_SAMBA_old.txt |  project  -C$start_lone/$start_late\
 -E$end_lone/$end_late -W-$width/$width -Q -Fpz > station_proj_so_5.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_SAMBA.txt |  project  -C$start_lone/$start_late\
 -E$end_lone/$end_late -W-$width/$width -Q -Fpz > station_proj_sn_5.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_WIZARD.txt |  project  -C$start_lone/$start_late \
-E$end_lone/$end_late -W-$width/$width -Q -Fpz > station_proj_w_5.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_GEONET.txt |  project  -C$start_lone/$start_late\
 -E$end_lone/$end_late -W-$width/$width -Q -Fpz > station_proj_g_5.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_DFDP10.txt |  project  -C$start_lone/$start_late \
-E$end_lone/$end_late -W-$width/$width -Q -Fpz > station_proj_dfdpa_5.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_DFDP13.txt |  project  -C$start_lone/$start_late \
-E$end_lone/$end_late -W-$width/$width -Q -Fpz > station_proj_dfdpb_5.dat


awk '{print($1,$2,$3)}' projection5.dat | psxy -Sci -i0,1,2s0.048 -W.25 -Gdimgrey \
-R0/100/-3/25 -J -B -Y-9 -X0 -O -K >> $out

psxy station_proj_a_5.dat -Si.4 -W.5 -Gcyan -R -J -O -K >> $out
psxy station_proj_sn_5.dat -Si.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_so_5.dat -St.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_w_5.dat -Si.4 -W.5 -Gmediumpurple -R -J -O -K >> $out
psxy station_proj_g_5.dat -Si.4 -W.5 -Gorange -R -J -O -K >> $out
psxy station_proj_dfdpa_5.dat -Si.4 -W.5 -Gyellow -R -J -O -K >> $out
psxy station_proj_dfdpb_5.dat -Si.4 -W.5 -Gdodgerblue -R -J -O -K >> $out

awk '{print($1,$2,$2)}' Cross_F_.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' Cross_G_.dat| psxy -Sd0.5  -W.5 -Gwhite -R \
-J -O -K -V >> $out

gmt pstext -R -J -O -K -F+f18p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
2 0 B
95.5 0 B'
END


# #..................
start_lonf='170.2'
start_latf='-43.15'
end_lonf='170.87'
end_latf='-43.85'
width='10'


##############################################################################
# Create a file that contains the cross section points in Lat/Lon
 awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_lonf/$start_latf \
-E$end_lonf/$end_latf -W-1/1 -Q -Fxyz > Cross_section_C.dat
#############################################################################'
# awk '{if ($4<25)print($2, $1, $3,$4+1)}' $datadir/loc_mag.dat | project -C$start_lonf/$start_latf -E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > projection6.dat

# awk '{ print $3, $2, $4 ,$17}' $datadir/../hypoDD.reloc33 | project -C$start_lonf/$start_latf -E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > projection6.dat
# awk '{if ($4>6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat |  project -C$start_lonf/$start_latf -E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > projection6.dat
# awk '{if ($4<6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat |  project -C$start_lonf/$start_latf -E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > projection6_.dat
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3  |  project -C$start_lonf/$start_latf \
-E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > projection6.dat

awk '{print $3, $2, "-1" , $1}' $datadir/sta_ALFA_08.txt |  project -C$start_lonf/$start_latf\
 -E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > station_proj_a_6.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_SAMBA.txt |  project  -C$start_lonf/$start_latf\
 -E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > station_proj_so_6.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_SAMBA_new.txt |  project  -C$start_lonf/$start_latf\
 -E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > station_proj_sn_6.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_WIZARD.txt |  project  -C$start_lonf/$start_latf \
-E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > station_proj_w_6.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_GEONET.txt |  project  -C$start_lonf/$start_latf \
-E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > station_proj_g_6.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_DFDP10.txt |  project  -C$start_lonf/$start_latf \
-E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > station_proj_dfdpa_6.dat
awk '{print $3, $2, "-1" , $1}' $datadir/sta_DFDP13.txt |  project  -C$start_lonf/$start_latf \
-E$end_lonf/$end_latf -W-$width/$width -Q -Fpz > station_proj_dfdpb_6.dat


awk '{print($1,$2,$3)}' projection6.dat | psxy -i0,1,2s0.048 -Sci -W.25p -Gdimgrey \
 -R0/100/-3/25 -B -J -Y-9 -X0 -O -K >> $out
psxy station_proj_a_6.dat -Si.4 -W.5 -Gcyan -R -J -O -K >> $out
psxy station_proj_so_6.dat -Si.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_sn_6.dat -St.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_w_6.dat -Si.4 -W.5 -Gmediumpurple -R -J -O -K >> $out
psxy station_proj_g_6.dat -Si.4 -W.5 -Gorange -R -J -O -K >> $out
psxy station_proj_dfdpa_6.dat -Si.4 -W.5 -Gyellow -R -J -O -K >> $out
psxy station_proj_dfdpb_6.dat -Si.4 -W.5 -Gdodgerblue -R -J -O -K >> $out

gmt pstext -R -J -O -K -F+f18p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
2 0 C
95.5 0 C'
END

# #..................
start_long='170'
start_latg='-43.25'
end_long='170.67'
end_latg='-43.95'
width='10'

##############################################################################
# Create a file that contains the cross section points in Lat/Lon
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_long/$start_latg \
-E$end_long/$end_latg -W-1/1 -Q -Fxyz > Cross_section_D.dat
#############################################################################'

# awk '{if ($4<25)print($2, $1, $3,$4+1)}' $datadir/loc_mag.dat | project -C$start_long/$start_latg -E$end_long/$end_latg -W-$width/$width -Q -Fpz > projection7.dat

# awk '{ print $3, $2, $4 ,$17}' $datadir/../hypoDD.reloc3 | project -C$start_long/$start_latg -E$end_long/$end_latg -W-$width/$width -Q -Fpz > projection7.dat
# awk '{if ($4>6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat | project -C$start_long/$start_latg -E$end_long/$end_latg -W-$width/$width -Q -Fpz > projection7.dat
# awk '{if ($4<6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat | project -C$start_long/$start_latg -E$end_long/$end_latg -W-$width/$width -Q -Fpz > projection7_.dat
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_long/$start_latg -E$end_long/$end_latg\
 -W-$width/$width -Q -Fpz > projection7.dat
awk '{print($1, $2, -2)}' Aoraki.dat | project -C$start_long/$start_latg -E$end_long/$end_latg\
 -W-$width/$width -Q -Fpz > Aoraki_d.dat
awk '{print($1, $2, $3)}' LFE_LMB.txt | project -C$start_long/$start_latg -E$end_long/$end_latg\
 -W-$width/$width -Q -Fpz > LFE_proj_g.dat



awk '{print $3, $2, "-1" , $1}' sta_ALFA_08.txt |  project -C$start_long/$start_latg \
-E$end_long/$end_latg -W-$width/$width -Q -Fpz > station_proj_a_7.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_new.txt |  project  -C$start_long/$start_latg\
 -E$end_long/$end_latg -W-$width/$width -Q -Fpz > station_proj_sn_7.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_old.txt |  project  -C$start_long/$start_latg\
 -E$end_long/$end_latg -W-$width/$width -Q -Fpz > station_proj_so_7.dat
awk '{print $3, $2, "-1" , $1}' sta_WIZARD.txt |  project  -C$start_long/$start_latg \
-E$end_long/$end_latg -W-$width/$width -Q -Fpz > station_proj_w_7.dat
awk '{print $3, $2, "-1" , $1}' sta_GEONET.txt |  project  -C$start_long/$start_latg \
-E$end_long/$end_latg -W-$width/$width -Q -Fpz > station_proj_g_7.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP10.txt |  project  -C$start_long/$start_latg\
 -E$end_long/$end_latg -W-$width/$width -Q -Fpz > station_proj_dfdpa_7.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP13.txt |  project  -C$start_long/$start_latg \
-E$end_long/$end_latg -W-$width/$width -Q -Fpz > station_proj_dfdpb_7.dat




awk '{print($1,$2,$3)}' projection7.dat | psxy -i0,1,2s0.048 -Sci -W.25p -Gdimgrey  \
-R0/100/-3/25 -J -B -Y-9 -X0 -O -K >> $out
awk '{print($1,$2,$2)}' Aoraki_d.dat| psxy -Sd0.5  -W.5 -Gblack -R0/100/-3/25 \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' LFE_proj_g.dat | psxy -Sa0.4  -W.25 -Cseis.cpt -R -J -O -K -V >> $out

psxy station_proj_a_7.dat -Si.4 -W.5 -Gcyan -R -J -O -K >> $out
psxy station_proj_so_7.dat -Si.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_sn_7.dat -St.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_w_7.dat -Si.4 -W.5 -Gmediumpurple -R -J -O -K >> $out
psxy station_proj_g_7.dat -Si.4 -W.5 -Gorange -R -J -O -K >> $out
psxy station_proj_dfdpa_7.dat -Si.4 -W.5 -Gyellow -R -J -O -K >> $out
psxy station_proj_dfdpb_7.dat -Si.4 -W.5 -Gdodgerblue -R -J -O -K >> $out


gmt pstext -R -J -O -K -F+f18p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
2 0 D
95.5 0 D'
END



# #..................
start_lonh='169.8'
start_lath='-43.35'
end_lonh='170.47'
end_lath='-44.05'
width='10'

##############################################################################
# Create a file that contains the cross section points in Lat/Lon
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-1/1 -Q -Fxyz > Cross_section_E.dat
#############################################################################'

# awk '{if ($4<25)print($2, $1, $3,$4+1)}' $datadir/loc_mag.dat | project -C$start_lonh/$start_lath -E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > projection8.dat

# awk '{if ($4>6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat | project -C$start_lonh/$start_lath -E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > projection8.dat
# awk '{if ($4<6.0) print($2, $1, $3,$4)}' reloc_mag_nll.dat | project -C$start_lonh/$start_lath -E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > projection8_.dat
awk '{ print $3, $2, $4 ,$17}' ../hypoDD.reloc3 | project -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > projection8.dat
awk '{print($1, $2, -2)}' Aoraki.dat | project -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > Aoraki_e.dat
awk '{print($1, $2, $3)}' LFE_LMB.txt | project -C$start_lonh/$start_lath -E$end_lonh/$end_lath\
 -W-$width/$width -Q -Fpz > LFE_proj_h.dat


awk '{print $3, $2, "-1" , $1}' sta_ALFA_08.txt |  project -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > station_proj_a_8.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_new.txt |  project  -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > station_proj_sn_8.dat
awk '{print $3, $2, "-1" , $1}' sta_SAMBA_old.txt |  project  -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > station_proj_so_8.dat
awk '{print $3, $2, "-1" , $1}' sta_WIZARD.txt |  project  -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > station_proj_w_8.dat
awk '{print $3, $2, "-1" , $1}' sta_GEONET.txt |  project  -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > station_proj_g_8.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP10.txt |  project  -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > station_proj_dfdpa_8.dat
awk '{print $3, $2, "-1" , $1}' sta_DFDP13.txt |  project  -C$start_lonh/$start_lath \
-E$end_lonh/$end_lath -W-$width/$width -Q -Fpz > station_proj_dfdpb_8.dat

new_projection="-JX20/-5"
new_boundaries="-Bf10a10/f5a5:"":wEnS"

awk '{print($1,$2,$3)}' projection8.dat | psxy -i0,1,2s0.048 -Sci -W.25p -Gdimgrey\
 -R0/100/-3/25 -J -Bx10+l"Distance (km)" -By5+l"Depth (km)" -BwSnE -Y-9 -X0 -O -K >> $out
awk '{print($1,$2,$2)}' Aoraki_e.dat| psxy -Sd0.5  -W.5 -Gblack -R0/100/-3/25 \
-J -O -K -V >> $out
awk '{print($1,$2,$2)}' LFE_proj_h.dat | psxy -Sa0.4  -W.25 -Cseis.cpt -R -J -O -K -V >> $out

psxy station_proj_a_8.dat -Si.4 -W.5 -Gcyan -R -J -O -K >> $out
psxy station_proj_so_8.dat -Si.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_sn_8.dat -St.4 -W.5 -Gred -R -J -O -K >> $out
psxy station_proj_w_8.dat -Si.4 -W.5 -Gmediumpurple -R -J -O -K >> $out
psxy station_proj_g_8.dat -Si.4 -W.5 -Gorange -R -J -O -K >> $out
psxy station_proj_dfdpa_8.dat -Si.4 -W.5 -Gyellow -R -J  -O -K >> $out
psxy station_proj_dfdpb_8.dat -Si.4 -W.5 -Gdodgerblue -R -J  -O -K >> $out

gmt pstext -R -J -O -K -F+f18p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
2 0 E
95.5 0 E'
END

psxy -R -J -T -O >> $out
ps2raster -Tf -A $out
#ps2raster -Tf -A map.ps
evince ${out%.*}.pdf



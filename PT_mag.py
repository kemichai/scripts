"""
P and T plunges and etc with magnitudes....
read Dropbox/VUW/PhD/Codes/Kaverina_plots/output-file.dat
"""
import os


directory = "/Users/home/michaiko/Dropbox/VUW/PhD/Codes/Kaverina_plots/"


print '>>> Reading focal mechanism details...'
working_dir = "/Users/home/michaiko/Dropbox/VUW/PhD/Codes/Kaverina_plots/"
fm_details = os.path.join(os.path.join(working_dir,
                          "output-file.dat"))

lon = []
lat = []
dep = []
mag = []
Pplu = []
Tplu = []
with open(fm_details, "r") as open_file:
    for line in open_file:
        Lon, Lat, Dep, _, _, _, _, _, _, _, \
        _, Mw, Strike_A, Dip_A, Rake_A, Strike_B, \
        Dip_B, Rake_B, P_Trend, P_plunge, B_trend, B_plunge, T_trend, \
        T_plunge, fclvd, X_Kaverina_diagram, Y_Kaverina_diagram, ID, \
        Mech_type = line.split()
        lon.append(Lon)
        mag.append(Mw)
        Pplu.append(P_plunge)
        Tplu.append(T_plunge)



mag = [float(i) for i in mag]
Pplu = [float(i) for i in Pplu]
Tplu = [float(i) for i in Tplu]


        
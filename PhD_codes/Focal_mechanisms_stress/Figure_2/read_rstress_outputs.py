#!/usr/bin/env python
"""
Read rstress outputs and store a file with each bin's
details as well as a file to plot the bowties with GMT.

VUW
July 2019
Konstantinos Michailos
"""

import os
from obspy import read_events, Catalog
from pyproj import Proj, transform
from functions import read_quakeml2create_lists
from functions import *
import numpy as np
import csv

# Define working directory
WORK_DIR = os.getcwd()

rstress_root = WORK_DIR + '/rstress/rcodes'

# Read cluster details file
cluster_details_file = WORK_DIR + '/rstress/details/cluster_details.csv'
c_reader = csv.reader(open(cluster_details_file, 'r'), delimiter=',')
cluster_name = list(zip(*c_reader))[0]
c_reader = csv.reader(open(cluster_details_file, 'r'), delimiter=',')
lon = list(zip(*c_reader))[1]
c_reader = csv.reader(open(cluster_details_file, 'r'), delimiter=',')
lat = list(zip(*c_reader))[2]
c_reader = csv.reader(open(cluster_details_file, 'r'), delimiter=',')
dep = list(zip(*c_reader))[3]
c_reader = csv.reader(open(cluster_details_file, 'r'), delimiter=',')
num_of_obs = list(zip(*c_reader))[4]

cluster_name = [int(i) for i in cluster_name]
lon = [float(i,) for i in lon]
lat = [float(i)for i in lat]
dep = [float(i) for i in dep]
num_of_obs = [int(i) for i in num_of_obs]

# Create a new file which includes the stress inversion details
with open(WORK_DIR + '/rstress/details/'
            'final_cluster_details.csv', 'a') as of:
    writer = csv.DictWriter(of, fieldnames=["Cluster", "Lat",
                                            "Lon","Dep", "N", "S1",
                                            "S2",
                                            "S3", "nu",
                                            "Shmax"]
                            , delimiter=',')
    writer.writeheader()

for i, val in enumerate(cluster_name):
    Rstress_out = '{}/outdata/{}'.format(rstress_root,
                                         'cluster_' + str(val) +
                                         '.1dparameters.dat')
    Sigma_out = '{}/outdata/{}'.format(rstress_root,
                                       'cluster_' + str(val) +
                                       '.2dparameters.dat')
    # Only if we have at least 15 observations
    if num_of_obs[i] >= 15:
        try:
            c_reader = csv.reader(open(Rstress_out, 'r'), delimiter=',')
            par = list(zip(*c_reader))[0]
            c_reader = csv.reader(open(Rstress_out, 'r'), delimiter=',')
            mean = list(zip(*c_reader))[1]
            c_reader = csv.reader(open(Rstress_out, 'r'), delimiter=',')
            X10 = list(zip(*c_reader))[4]
            c_reader = csv.reader(open(Rstress_out, 'r'), delimiter=',')
            X90 = list(zip(*c_reader))[5]

            Shmax_mean = float(mean[4])
            Shmax_10 = float(X10[4])
            Shmax_90 = float(X90[4])
            Shmax_m = Shmax_mean - Shmax_10
            Shmax_p = Shmax_90 - Shmax_mean

            nu_mean = float(mean[5])
            nu_10 = float(X10[5])
            nu_90 = float(X90[5])
            nu_m = nu_mean - nu_10
            nu_p = nu_90 - nu_mean

            c_reader = csv.reader(open(Sigma_out, 'r'), delimiter=',')
            par = list(zip(*c_reader))[0]
            c_reader = csv.reader(open(Sigma_out, 'r'), delimiter=',')
            mean = list(zip(*c_reader))[1]

            S1_phi = float(mean[5])
            S1_theta = float(mean[6])
            S2_phi = float(mean[7])
            S2_theta = float(mean[8])
            S3_phi = float(mean[9])
            S3_theta = float(mean[10])

            def phi_theta2trend_plunge(phi, theta):
                # Sort out upwards vectors
                if theta > 90:
                    plunge = 180. - theta
                    if phi < 0:
                        trend = 180. + phi
                    else:
                        trend = phi
                else:
                    plunge = theta
                    if phi < 0:
                        trend = 180. + phi
                    else:
                        trend = phi

                return trend, plunge


            S1_trend = phi_theta2trend_plunge(S1_phi, S1_theta)[0]
            S1_plunge = phi_theta2trend_plunge(S1_phi, S1_theta)[1]
            S2_trend = phi_theta2trend_plunge(S2_phi, S2_theta)[0]
            S2_plunge = phi_theta2trend_plunge(S2_phi, S2_theta)[1]
            S3_trend = phi_theta2trend_plunge(S3_phi, S3_theta)[0]
            S3_plunge = phi_theta2trend_plunge(S3_phi, S3_theta)[1]


            def gmt_azim(value):
                if value <= 180:
                    out = value
                    out_2 = out + 180
                else:
                    out_2 = out - 180
                    out = value
                return out, out_2


            az_gmt10, az_gmt10_ = gmt_azim(Shmax_10)
            az_gmt90, az_gmt90_ = gmt_azim(Shmax_90)

            with open(WORK_DIR + '/gmt_files/' +
                      'gmt_psxy_bowties.dat', 'a') as f:
                f.write('{} {} {} {}\n {} {} {} {}\n'.format(lon[i], lat[i],
                                                             az_gmt10, az_gmt90, lon[i], lat[i], az_gmt10_, az_gmt90_))

            with open(WORK_DIR + '/rstress/details/' +
                                  'final_cluster_details.csv', 'a') as of:
                of.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(cluster_name[i], round(lat[i], 2),
                                                                           round(lon[i], 2), round(dep[i], 1),
                                                                           num_of_obs[i],
                                                                           str(round(S1_trend, 1)) + '/' + str(
                                                                               round(S1_plunge)),
                                                                           str(round(S2_trend, 1)) + '/' + str(
                                                                               round(S2_plunge)),
                                                                           str(round(S3_trend, 1)) + '/' + str(
                                                                               round(S3_plunge)),
                                                                           str(round(nu_mean, 1)) + ' (-' + str(
                                                                               round(nu_m, 1)) + '/+' + str(
                                                                               round(nu_p, 1)) + ')',
                                                                           str(round(Shmax_mean, 1)) + ' (-' + str(
                                                                               round(Shmax_m, 1)) + '/+' + str(
                                                                               round(Shmax_p, 1)) + ')'))
        except Exception as e:
            print(e)
    else:
        continue








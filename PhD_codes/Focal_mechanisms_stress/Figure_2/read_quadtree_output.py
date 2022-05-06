#!/usr/bin/env python
"""
Reads quadtree output and creates input for rstress codes.

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

# Read matlab outputs
cat_file = (WORK_DIR + '/matlab/output/'
            'cat_box_num.csv')

# Read csv files
cat_csv_file = (cat_file)
c_reader = csv.reader(open(cat_csv_file, 'r'), delimiter=',')
cluster_names_str = list(zip(*c_reader))[-1]
c_reader = csv.reader(open(cat_csv_file, 'r'), delimiter=',')
error_str = list(zip(*c_reader))[-3]
c_reader = csv.reader(open(cat_csv_file, 'r'), delimiter=',')
rake_str = list(zip(*c_reader))[-4]
c_reader = csv.reader(open(cat_csv_file, 'r'), delimiter=',')
dip_str = list(zip(*c_reader))[-5]
c_reader = csv.reader(open(cat_csv_file, 'r'), delimiter=',')
strike_str = list(zip(*c_reader))[-6]

strike = convert_str2float(strike_str)
dip = convert_str2float(dip_str)
rake = convert_str2float(rake_str)
error = convert_str2float(error_str)
cluster_names = convert_str2int(cluster_names_str)


cl_list = []
for i in cluster_names:
    if i not in cl_list:
        cl_list.append(i)
ENCODING = 'utf-8'
clusters = [[] for _ in range(len(cl_list))]
cl_list.sort()

for i, val in enumerate(cl_list):
    # Create a file for each
    # cluster and append the header...
    print('Cluster ' + str(val))
    with open(WORK_DIR + '/rstress/rcodes/indata/'
              'cluster_' + str(val) + '.csv', 'a') as of:
        writer = csv.DictWriter(of, fieldnames=["strike", "dip",
                                                "rake", "err"], delimiter=',')
        writer.writeheader()
    for k, l in enumerate(strike):
        if cluster_names[k] == val:
            with open(WORK_DIR + '/rstress/rcodes/indata/'
                      'cluster_'+str(val)+'.csv', 'a') as of:
                of.write('{}, {}, {}, {}\n'.format(strike[k],
                         dip[k], rake[k], error[k]))
# to copy paste in drive_inversion.q for running Rstress command...
for i, val in enumerate(cl_list):
    with open(WORK_DIR + '/rstress/'
              'Rstress_command.csv', 'a') as of:
        of.write('{}\n'.format('"cluster_' + str(val)+'.csv",'))

##############################
# Information on the details
# of the clusters
cluster_file = (WORK_DIR + '/matlab/output/'
                'cluster_details.csv')
cluster_det_csv_file = (cluster_file)
# Read cluster details
c_reader = csv.reader(open(cluster_det_csv_file, 'r'), delimiter=',')
Y_str = list(zip(*c_reader))[1]
c_reader = csv.reader(open(cluster_det_csv_file, 'r'), delimiter=',')
X_str = list(zip(*c_reader))[0]
c_reader = csv.reader(open(cluster_det_csv_file, 'r'), delimiter=',')
number_of_eqs_str = list(zip(*c_reader))[2]
c_reader = csv.reader(open(cluster_det_csv_file, 'r'), delimiter=',')
median_depths_str = list(zip(*c_reader))[3]

X = convert_str2float(X_str)
Y = convert_str2float(Y_str)
number_of_eqs = convert_str2int(number_of_eqs_str)
median_depths = convert_str2float(median_depths_str)


cluster_nam = []
for i, j in enumerate(number_of_eqs):
    print(i+1, j)
    cluster_nam.append(i+1)


PROJECTION_IN = Proj(init='epsg:4326')
PROJECTION_OUT = Proj(init='epsg:2193')
DEGREES_2_ROTATE = -54

orig = [0, 0]
rot = CoordRotator(orig, np.radians(DEGREES_2_ROTATE))

cent_lon = []
cent_lat = []
for i, j in enumerate(X):
    aX, aY = rot.inverse(X[i], Y[i])
    lon, lat = transform(PROJECTION_OUT, PROJECTION_IN, aX, aY)
    cent_lat.append(lat)
    cent_lon.append(lon)
    # to use it later with the rstress outputs
    with open(WORK_DIR + '/rstress/details/'
              'cluster_details.csv', 'a') as of:
        of.write('{}, {}, {}, {}, {}\n'.format(cluster_nam[i],
                 lon, lat, median_depths[i], number_of_eqs[i]))
    # to use it for GMT plot
    with open(WORK_DIR + '/gmt_files/'
              'box_names.dat', 'a') as of:
        of.write('{} {} {}\n'.format(lon, lat, cluster_nam[i]))
# Run rstress codes...


# More GMT stuff (plot bins colored according to the seismic moment release
quadtree2GMT(boxes_file=WORK_DIR + '/matlab/output/final_boxes.csv',
             values_file=WORK_DIR + '/matlab/output/final_values.csv',
             output_path=WORK_DIR + '/gmt_files')

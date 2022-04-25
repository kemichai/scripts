#!/usr/bin/env python
"""
Read quakeml catalog of earthquakes with focal mechanisms and
prepares input file for the Matlab quadtree code.

Converts earthquake location's projection from lat/lon to NZGD2000
(New Zealand Transverse Mercator 2000) and rotates the
earthquake locations to align with Alpine Fault's strike.

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

# Define working directory
WORK_DIR = os.getcwd()

# Read catalog with focal mechanisms (quakeML format)
# TODO: do the same for a text file
# COMMENTED; USE FOR FOCAL MECH CATALOG
# CAT_PATH = WORK_DIR + '/quakeml_file/subset_cat.xml'

CAT_PATH = WORK_DIR + '/quakeml_file/Southern_Alps_NZ_microseismicity_catalog_HDD.xml'

print('>>> Reading catalog...')
catalog = read_events(CAT_PATH)
print('>>> Done.')

# Create lists from quakeml (lon, lat, dep, etc.)
# COMMENTED; USE FOR FOCAL MECH CATALOG
# lon, lat, dep, moment, strike, dip, rake, error, ev_id = read_quakeml2create_lists(catalog)

# Create lists from quakeml (lon, lat, dep, etc.)
lon, lat, dep, ev_id = read_quakeml2create_lists_alt(catalog)

##############################
# Convert and rotate locations
# Define parameters
PROJECTION_IN = Proj(init='epsg:4326')
PROJECTION_OUT = Proj(init='epsg:2193')
DEGREES_2_ROTATE = -54
# Convert from lat/lon to NZGD2000 (New Zealand Transverse Mercator 2000)
nz_x, nz_y, z = transform(PROJECTION_IN, PROJECTION_OUT, lon, lat, dep)

# Define origin point for rotation
orig = [0, 0]
nz_x_orig = np.asarray(nz_x) + orig[0]
nz_y_orig = np.asarray(nz_y) + orig[1]

# Rotate parallel to AF's strike
rot = CoordRotator(orig, np.radians(-54))
nz_x_rot, nz_y_rot = rot.forward(nz_x_orig, nz_y_orig)

# Write text file (input on matlab code)
# COMMENTED; USE FOR FOCAL MECH CATALOG
# for i, dat in enumerate(ev_id):
#     with open(WORK_DIR + '/matlab/input/'
#                          'matlab_quadtree_input.dat', 'a') as f:
#         f.write('{} {} {} {} {} {} {} {} {}\n'.format(nz_x_rot[i], nz_y_rot[i],
#                                                       dep[i], moment[i], strike[i],
#                                                       dip[i], rake[i], error[i],
#                                                       ev_id[i]))

# Write text file (input on matlab code)
for i, dat in enumerate(ev_id):
    with open(WORK_DIR + '/matlab/input/'
                         'matlab_quadtree_input.dat', 'a') as f:
        f.write('{} {} {} {}\n'.format(nz_x_rot[i], nz_y_rot[i],
                                                      dep[i],
                                                      ev_id[i]))

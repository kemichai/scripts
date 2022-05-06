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
year = list(zip(*c_reader))[-2]
c_reader = csv.reader(open(cat_csv_file, 'r'), delimiter=',')
dep = list(zip(*c_reader))[-3]

eqz_year = convert_str2int(year)
eqz_dep = convert_str2float(dep)
cluster_names = convert_str2int(cluster_names_str)

cl_list = []
for i in cluster_names:
    if i not in cl_list:
        cl_list.append(i)
ENCODING = 'utf-8'
clusters = [[] for _ in range(len(cl_list))]
cl_list.sort()


for yr in range(2008, 2018, 1):
    for j, cl in enumerate(cl_list):
        eqz_within_box_yearly = []
        for i, eq_dep in enumerate(eqz_dep):
            if eqz_year[i] == yr and cluster_names[i] == cl:
                eqz_within_box_yearly.append(eq_dep)
        if len(eqz_within_box_yearly) > 15:
            med_dep = np.median(eqz_within_box_yearly)
            print(str(yr) + ', ' + str(med_dep) + ', ' + str(cl))

# Will need to work on this to store the outputs in a more useful format

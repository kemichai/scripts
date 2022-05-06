#!/usr/bin/env python
"""
VUW
July 2019
Konstantinos Michailos
"""

from pyproj import Proj, transform
from obspy import read_events
import numpy as np
import csv

# class used to rotate coordinates
class CoordRotator:
    def __init__(self, origin, angle):
        self.origin = np.asarray(origin, dtype=np.float64)
        self.angle = float(angle)

    def forward(self, x, y):
        rx = np.asarray(x, dtype=np.float64) - self.origin[0]
        ry = np.asarray(y, dtype=np.float64) - self.origin[1]
        ca = np.cos(self.angle)
        sa = np.sin(self.angle)
        xx = rx * ca + ry * sa
        yy = -rx * sa + ry * ca
        return xx, yy

    def inverse(self, xx, yy, z=None):
        ca = np.cos(self.angle)
        sa = np.sin(self.angle)
        rx = xx * ca + -yy * sa
        ry = xx * sa + yy * ca
        x = rx + self.origin[0]
        y = ry + self.origin[1]
        return x, y


def read_quakeml2create_lists(cat):
    """
    *** BEWARE OF HARD CODING BELOW ***
        reads last origin and first magnitude
    """
    dep = []
    moment = []
    mag = []
    lat = []
    lon = []
    strike = []
    dip = []
    rake = []
    error = []
    ev_id = []
    for ev in cat:
        orig = str(ev.origins[-1].time)
        ymd = orig[0:4] + orig[5:7] + orig[8:10]
        hmss = orig[11:13] + orig[14:16] + orig[17:19] + orig[20:22]
        time = ymd + '.' + hmss
        ev_id.append(time)
        dep.append(ev.origins[-1].depth/1000)
        lat.append(ev.origins[-1].latitude)
        lon.append(ev.origins[-1].longitude)

        magnitude = ev.magnitudes[0].mag
        # We calculate, Mo, seismic moment according to Hanks and
        # Kanamori, 1978 relationship in Nm.
        seismic_moment = 10**(1.5*magnitude + 9)
        # seismic_moment_log = np.log10(seismic_moment)
        moment.append(seismic_moment)
        mag.append(magnitude)
        strike.append(ev.focal_mechanisms[0].nodal_planes.nodal_plane_1.strike)
        dip.append(ev.focal_mechanisms[0].nodal_planes.nodal_plane_1.dip)
        rake.append(ev.focal_mechanisms[0].nodal_planes.nodal_plane_1.rake)
        error.append(str(float(ev.focal_mechanisms[0].comments[0].text)))

    return lon, lat, dep, moment, strike, dip, rake, error, ev_id


def read_quakeml2create_lists_alt(cat):
    """
    *** BEWARE OF HARD CODING BELOW ***
        reads last origin and first magnitude
    """
    dep = []
    moment = []
    mag = []
    lat = []
    lon = []
    ev_id = []

    for ev in cat:
        orig = str(ev.origins[-1].time)
        ymd = orig[0:4] + orig[5:7] + orig[8:10]
        hmss = orig[11:13] + orig[14:16] + orig[17:19] + orig[20:22]
        time = int(orig[0:4])
        ev_id.append(time)
        dep.append(ev.origins[-1].depth/1000)
        lat.append(ev.origins[-1].latitude)
        lon.append(ev.origins[-1].longitude)

        magnitude = ev.magnitudes[0].mag
        # We calculate, Mo, seismic moment according to Hanks and
        # Kanamori, 1978 relationship in Nm.
        seismic_moment = 10**(1.5*magnitude + 9)
        # seismic_moment_log = np.log10(seismic_moment)
        moment.append(seismic_moment)
        mag.append(magnitude)

    return lon, lat, dep, ev_id



def convert_str2float(input_list):
    output_list = [float(i) for i in input_list]
    return output_list


def convert_str2int(input_list):
    output_list = [int(i) for i in input_list]
    return output_list

def dist_calc(loc1, loc2):
    """
    Function to calculate the distance in km between two points.

    Uses the flat Earth approximation. Better things are available for this,
    like `gdal <http://www.gdal.org/>`_.

    :type loc1: tuple
    :param loc1: Tuple of lat, lon, depth (in decimal degrees and km)
    :type loc2: tuple
    :param loc2: Tuple of lat, lon, depth (in decimal degrees and km)

    :returns: Distance between points in km.
    :rtype: float
    :author: Calum Chamberlain
    """
    R = 6371.009  # Radius of the Earth in km
    dlat = np.radians(abs(loc1[0] - loc2[0]))
    dlong = np.radians(abs(loc1[1] - loc2[1]))
    ddepth = abs(loc1[2] - loc2[2])
    mean_lat = np.radians((loc1[0] + loc2[0]) / 2)
    dist = R * np.sqrt(dlat ** 2 + (np.cos(mean_lat) * dlong) ** 2)
    dist = np.sqrt(dist ** 2 + ddepth ** 2)
    return dist


def quadtree2GMT(boxes_file, values_file, output_path, degrees=-54,
                 inProject=Proj(init='epsg:4326'),
                 outProject=Proj(init='epsg:2193')):
    """
    Function to read the outputs from the quadtree codes.
    Rotate and convert back to lat/lon.
    Write an input file for GMT.
    GMT command to plot the boxes....
    psxy boxes_gmt.dat -R -J -W0.25p -O -K -L -Cseis.cpt -t0  >> $out

    Arguments:
        boxes_file {str}  -- Path to the file containing the edges of the boxes
        values_file {str} -- Path to the file with the values of each box
        degrees {float}   -- Degrees for the coordinates to be rotated
                             (-54 will make it parallel to the Alpine
                             Fault's strike)
        inProject {str}   -- Input projection. Default is (epsg:4326)
        outProject {str}  -- Output projection. Default is (epsg:2193)
    """

    # Read files
    boxes_csv_file = (boxes_file)
    values_csv_file = (values_file)

    c_reader = csv.reader(open(boxes_csv_file, 'r'), delimiter=',')
    X2 = list(zip(*c_reader))[1]
    c_reader = csv.reader(open(boxes_csv_file, 'r'), delimiter=',')
    X1 = list(zip(*c_reader))[0]
    c_reader = csv.reader(open(boxes_csv_file, 'r'), delimiter=',')
    Y1 = list(zip(*c_reader))[2]
    c_reader = csv.reader(open(boxes_csv_file, 'r'), delimiter=',')
    Y2 = list(zip(*c_reader))[3]

    # four edges of each box
    X1 = np.asarray([float(i) for i in X1])
    X2 = np.asarray([float(i) for i in X2])
    Y1 = np.asarray([float(i) for i in Y1])
    Y2 = np.asarray([float(i) for i in Y2])

    # values for each box
    c_reader = csv.reader(open(values_csv_file, 'r'), delimiter=',')
    avdep = list(zip(*c_reader))[0]
    avdep = np.asarray([float(i) for i in avdep])

    # Define the projections
    inProj = inProject
    outProj = outProject
    orig = [0, 0]
    aaa = CoordRotator(orig, np.radians(-54))

    for i, j in enumerate(X1):
        # print 'gmt psxy -R -J -Wthin,red -O -K  >> $out << END'
        aX, aY = aaa.inverse(X1[i], Y1[i])
        bX, bY = aaa.inverse(X2[i], Y2[i])
        cX, cY = aaa.inverse(X2[i], Y1[i])
        dX, dY = aaa.inverse(X1[i], Y2[i])

        cx, cy = transform(outProj, inProj, cX, cY)
        bx, by = transform(outProj, inProj, bX, bY)
        dx, dy = transform(outProj, inProj, dX, dY)
        ax, ay = transform(outProj, inProj, aX, aY)
        #############################################################  #
        C = [cx, cy, 0]
        B = [bx, by, 0]
        D = [dx, dy, 0]
        A = [ax, ay, 0]
        # Calculate the area of the bin, which is actually a rombus!
        area = ((dist_calc(C, A)*dist_calc(D,B))/2)
        ##############################################################
        if avdep[i] == 0:
            Mo = 0
            Mo_area = 0
        else:
            Mo = np.log10(avdep[i])
            Mo_area = np.log10(avdep[i]/area)
        # Write boxes in a format readable by gmt
        with open(output_path + '/boxes_gmt.dat', 'a') as of:
            of.write('{}{}\n\
                    {} {}\n\
                    {} {}\n\
                    {} {}\n\
                    {} {}\n\
                    {} {}\n'.format('>-Z', str(Mo), cx, cy, bx, by, dx,
                                    dy, ax, ay, cx, cy))
        with open(output_path + '/box_values.dat', 'a') as of:
            of.write('{} {} {} {}\n'.format(str(Mo),str(Mo_area), cx, cy))

        with open(output_path + '/boxes_gmt_2area.dat', 'a') as of:
            of.write('{}{}\n\
                    {} {}\n\
                    {} {}\n\
                    {} {}\n\
                    {} {}\n\
                    {} {}\n'.format('>-Z', str(Mo_area), cx, cy, bx, by, dx,
                                    dy, ax, ay, cx, cy))
    return

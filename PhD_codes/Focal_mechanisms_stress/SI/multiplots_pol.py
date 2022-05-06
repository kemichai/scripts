#!/usr/bin/env python
"""
Read quakeml catalog of earthquakes with focal mechanisms and
plot details.

May 2022
Konstantinos Michailos
"""
from obspy import read_events
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

# Define current directory
WORK_DIR = os.getcwd()

CAT_PATH = '/'.join(WORK_DIR.split('/')[0:-1]) + '/quakeml_file/Final_foc_mec_catalog.xml'

print('>>> Reading catalog...')
cat = read_events(CAT_PATH)
print('>>> Done.')

strike = []
dip = []
rake = []
azim = []
npol = []
unce = []
mag = []
dep = []
for ev in cat:
    strike.append(float(ev.focal_mechanisms[0].nodal_planes.nodal_plane_1.strike))
    dip.append(float(ev.focal_mechanisms[0].nodal_planes.nodal_plane_1.dip))
    rake.append(float(ev.focal_mechanisms[0].nodal_planes.nodal_plane_1.rake))
    azim.append(float(ev.focal_mechanisms[0].azimuthal_gap))
    npol.append(int(ev.focal_mechanisms[0].station_polarity_count))
    unce.append(float(ev.focal_mechanisms[0].comments[0].text))
    mag.append(float(ev.magnitudes[0].mag))
    dep.append(float(ev.origins[-1].depth)/1000)
    
font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]

fig_size[0] = 14
fig_size[1] = 17
plt.rcParams["figure.figsize"] = fig_size


ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=1)
ax1.scatter(mag, unce, facecolor='darkgrey', alpha=0.6, edgecolor='black', linewidth=1.)
ax1.set_ylabel('Focal Mechanism \n error ($^\circ$)', fontsize=20)
ax1.set_xlabel('Magnitudes (M$_L$)', fontsize=20)
plt.xticks(np.arange(1.5, 5.5, 1))
plt.yticks(np.arange(20,45, 5))

ax2 = plt.subplot2grid((3, 3), (0, 1), colspan=1)
ax2.scatter(npol, unce, facecolor='darkgrey', alpha=0.6, edgecolor='black', linewidth=1.)
plt.yticks(np.arange(20,45, 5))
plt.xticks(np.arange(7, 40, 6))

ax3 = plt.subplot2grid((3, 3), (0, 2), colspan=1)
ax3.scatter(azim, unce, facecolor='darkgrey', alpha=0.6, edgecolor='black', linewidth=1.)
plt.xticks(np.arange(60, 361, 60))
plt.yticks(np.arange(20,45, 5))
#
ax4 = plt.subplot2grid((3, 3), (1, 1), colspan=1)
ax4.scatter(npol, mag, facecolor='darkgrey', alpha=0.6, edgecolor='black', linewidth=1.)
ax4.set_ylabel('Magnitudes (M$_L$)', fontsize=20)
ax4.set_xlabel('Number of \n polarities', fontsize=20)
plt.xticks(np.arange(7, 40, 6))
plt.yticks(np.arange(1.5, 5.5, 1))
#
ax5 = plt.subplot2grid((3, 3), (1, 2), colspan=1)
ax5.scatter(azim, mag, facecolor='darkgrey', alpha=0.6, edgecolor='black', linewidth=1.)
plt.yticks(np.arange(1.5, 5.5, 1))
plt.xticks(np.arange(60, 361, 60))
#
ax6 = plt.subplot2grid((3, 3), (2, 2), colspan=1)
ax6.scatter(azim, npol, facecolor='darkgrey', alpha=0.6, edgecolor='black',
            linewidth=1.)
ax6.set_xlabel('Focal Mechanism \n azimuthal gap ($^\circ$)', fontsize=20)
ax6.set_ylabel('Number of \n polarities', fontsize=20)
plt.yticks(np.arange(7, 40, 6))
plt.xticks(np.arange(60, 361, 60))
plt.show()


##
ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=1)
ax1.scatter(mag, rake, facecolor='darkgrey', alpha=0.6, edgecolor='black',
            linewidth=1.)
ax1.set_ylabel('Rake ($^\circ$)', fontsize=20)
ax1.set_xlabel('Magnitudes (M$_L$)', fontsize=20)
ax1.set_ylim([-180, 180])
ax1.set_ylim([1.5, 5.5])
plt.xticks(np.arange(1.5, 5.5, 1))
plt.yticks(np.arange(-180,181, 60))
#
ax2 = plt.subplot2grid((3, 3), (0, 1), colspan=1)
ax2.scatter(dep, rake, facecolor='darkgrey', alpha=0.6, edgecolor='black',
            linewidth=1.)
plt.yticks(np.arange(-180,181, 60))
plt.xticks(np.arange(0, 31, 10))
ax2.set_ylim([-180, 180])
ax2.set_xlim([0, 30])
#
ax4 = plt.subplot2grid((3, 3), (1, 1), colspan=1)
ax4.scatter(dep,mag, facecolor='darkgrey', alpha=0.6, edgecolor='black',
            linewidth=1.)
ax4.set_ylabel('Magnitudes (M$_L$)', fontsize=20)
ax4.set_xlabel('Hypocentral \n depths (km)', fontsize=20)
ax4.set_xlim([0, 30])
plt.xticks(np.arange(0, 31, 10))
plt.yticks(np.arange(1.5, 5.5, 1))
plt.show()


font = {'family': 'normal',
        'weight': 'normal',
        'size': 18}
matplotlib.rc('font', **font)
# Set figure width to 12 and height to 9
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 9
fig_size[1] = 7
plt.rcParams["figure.figsize"] = fig_size

bins = np.arange(-0.5, 40.5, 0.5)
ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1)
ax1.hist(unce, bins, histtype='step', orientation='vertical',
         color='black',facecolor='grey', alpha=0.9, linewidth=2, 
         edgecolor='k',fill=True)
ax1.set_ylim([0, 100])
ax1.set_xlim([0, 40])
ax1.set_xlabel('Focal Mechanism error ($^\circ$)', fontsize=18)
ax1.set_ylabel(r'Number of events', fontsize=18)
# plt.axhline(np.median(dep), color='k', linestyle='dashed', linewidth=2, label='Mean')
plt.axvline(np.mean(unce), color='k', linestyle='dashed', linewidth=2, 
            label='Mean (' +str(round(np.mean(unce),1)) + '$^\circ$)' )
plt.axvline(np.mean(unce)-np.std(unce), color='k',
            linestyle='dotted', linewidth=2,
            label='Std (' +str(round(np.std(unce),1)) + '$^\circ$)' )
plt.axvline(np.mean(unce)+np.std(unce), color='k',
            linestyle='dotted', linewidth=2)            
plt.legend(loc="upper left", markerscale=1., scatterpoints=1, fontsize=16)
plt.show()
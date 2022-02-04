""""
Plot data continuity of seismic waveforms.

Original codes by CJC
:Location: Chavannes-pres-renens, CH
:Date: Dec 2021
:Author: K. Michailos
"""
import glob
import numpy as np
import matplotlib.dates as mpdates
import matplotlib.pyplot as plt
import datetime as dt

# Set paths
path = '/Volumes/GeoPhysics_09/users-data/chambeca/SAMBA_archive/day_volumes_S'
WIZpath = '/Volumes/GeoPhysics_05/users-data/michaiko/WIZARD'
DFDpath = '/Volumes/GeoPhysics_09/users-data/chambeca/DFDP_archive/day_volumes_S'
DFDP10path = '/Volumes/GeoPhysics_05/users-data/michaiko/DFDP10/'
alfa08path = '/Volumes/GeoPhysics_05/users-data/michaiko/ALFA08/'


# Give station names
stations = ['COSA', 'EORO', 'MTFO', 'WHAT2', 'WHYM', 'POCR2',
            'LABE', 'GOVA', 'FRAN', 'REYN', 'COVA', 'LARB',
            'SOLU', 'MTBA',
            'WVZ', 'RPZ', 'FOZ', 'LBZ', 'JCZ',
            'WZ01', 'WZ02', 'WZ03', 'WZ04', 'WZ05', 'WZ06', 'WZ07',
            'WZ08', 'WZ09', 'WZ10', 'WZ11', 'WZ12', 'WZ13', 'WZ14',
            'WZ15', 'WZ16', 'WZ17', 'WZ18', 'WZ19', 'WZ20', 'WZ21',
            'BLLO', 'NOLA', 'ROTO', 'MTFE',
            'WDSZ', 'WMSZ', 'WPSZ', 'WTSZ', 'WBSZ',
            'BLO', 'DRC', 'GCK', 'NOL', 'POE', 'WHB', 'BON', 'ERE',
            'GHU', 'ONE', 'VBV', 'WNQ',
            'TURI', 'MTHA2', 'CAMEL', 'HATT', 'BURA2', 'CASH',
            'FBLA2', 'UMAT']


fig = plt.figure(num=1, dpi=100, facecolor='w', edgecolor='b')
fig.suptitle('Data Continuity', fontsize=30)
fig.set_size_inches(40, 30)

# Set min and max dates
mindate = mpdates.date2num(dt.datetime.strptime('20081101', '%Y%m%d'))
maxdate = mpdates.date2num(dt.datetime.strptime('20170301', '%Y%m%d'))

x = range(int(mindate), int(maxdate))
xdates = [mpdates.num2date(xd) for xd in x]
kstas = np.zeros(len(xdates))

plotno = 1
for sta in stations:
    y = []
    for i in xrange(len(xdates)):
        if glob.glob(path+'/'+xdates[i].strftime('Y%Y/R%j.01')+'/'+sta+'*'):
            y += [1]
            kstas[i] += 1
        elif glob.glob(WIZpath+'/'+xdates[i].strftime('Y%Y/R%j.01')+'/''*.'+sta+'*'):
            y += [1]
            kstas[i] += 1
        elif glob.glob(DFDpath+'/'+xdates[i].strftime('Y%Y/R%j.01')+'/'+sta+'*'):
            y += [1]
            kstas[i] += 1
        elif glob.glob(DFDP10path+'/'+xdates[i].strftime('Y%Y/R%j.01')+'/''*.'+sta+'*'):
            y += [1]
            kstas[i] += 1
        elif glob.glob(alfa08path+'/'+xdates[i].strftime('Y%Y/R%j.01')+'/''*.'+sta+'*'):
            y += [1]
            kstas[i] += 1
        else:
            y += [0]
    ax = fig.add_subplot(len(stations), 1, plotno)
    ax.scatter(xdates, y, c='k')

    ax.set_ylim([0.5, 1.5])
    ax.set_xlim([min(x), max(x)])
    ax.set_ylabel(sta, rotation=0, fontsize=18)
    ax.get_yaxis().set_label_coords(-0.015, 0.1)
    ax.xaxis.set_major_locator(mpdates.YearLocator())
    ax.xaxis.set_major_formatter(mpdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_minor_locator(mpdates.MonthLocator())

    plotno += 1
fig.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
ax.fmt_xdata = mpdates.DateFormatter('%Y-%m-%d')
fig.autofmt_xdate()
plt.show()


fig.savefig("cont_plot.svg", dpi=500, facecolor='w', edgecolor='w',
            orientation='landscape', format='svg',
            pad_inches=0.1)

#plt.plot(xdates,kstas,'b')
#plt.show()

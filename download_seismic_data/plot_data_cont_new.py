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
import matplotlib

# Set paths
paths = ['/media/kmichall/SEISMIC_DATA/Data_archive/GANSSER',
         '/media/kmichall/SEISMIC_DATA/Data_archive/HIMNT',
         '/media/kmichall/SEISMIC_DATA/Data_archive/GIC_XE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/Hi_CLIMB',
         '/media/kmichall/SEISMIC_DATA/Data_archive/BPE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/BSN',
         '/media/kmichall/SEISMIC_DATA/Data_archive/BSN',
         '/media/kmichall/SEISMIC_DATA/Data_archive/IC',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NSNI',
         '/media/kmichall/SEISMIC_DATA/Data_archive/PIRE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NAMASTE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NK_KKN',
         '/media/kmichall/SEISMIC_DATA/Data_archive/EVK2',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NQ']
# TODO: Need to fix the loop through the days to make it loop through the wav paths as well!!!

# stah = []
# with open('/home/kmichall/Desktop/Codes/bitbucket/him_seismicity/maps/files/sta_NAMASTE.txt', 'r') as f:
#     for line in f:
#         if line.startswith('#'):
#             print(line)
#             continue
#         else:
#             ln = line.split()
#             stah.append(ln[1])

# give station names

# HIMNT
# stations = ['BIRA', 'BUNG', 'DINX', 'GAIG', 'HILE', 'DINX', 'HILE', 'ILAM', 'JANA', 'JIRI', 'LAZE', 'MAZA',
#             'MNBU', 'NAIL', 'NAMC', 'ONRN', 'PHAP', 'PHID', 'RBSH', 'RC14', 'RUMJ', 'SAGA', 'SAJA', 'SIND',
#             'SSAN', 'SUKT', 'THAK', 'TUML', 'XIXI', 'YALA',]#
# stations = ['ES41', 'ES42', 'ES43']

# NAMASTE
# stations =['NA010',  'NA020', 'NA030', 'NA040', 'NA050', 'NA060', 'NA070', 'NA080', 'NA090', 'NA100', 'NA110',
           # 'NA120', 'NA130', 'NA150', 'NA160', 'NA170', 'NA180', 'NA190', 'NA200', 'NA210', 'NA220', 'NA230',
           #  'NA240', 'NA250', 'NA260', 'NA270', 'NA280', 'NA290', 'NA300', 'NA310', 'NA320', 'NA330', 'NA340',
            # 'NA350', 'NA360']
# stations = ['NA370', 'NA380', 'NA390', 'NA400', 'NA410', 'NA420', 'NA430', 'NA450', 'NABKT',
            # 'NAHTG', 'NAKRM', 'NANST']


# GANSSER
# stations = ['BHE01', 'BHE02', 'BHE03', 'BHE04', 'BHE05', 'BHE06',
#              'BHE07', 'BHE08', 'BHE09', 'BHE10', 'BHE11', 'BHE12', 'BHE13', 'BHE14', 'BHE15', 'BHN01', 'BHN02',
#              'BHN03', 'BHN04', 'BHN05', 'BHN06', 'BHN07', 'BHN08', 'BHN09', 'BHN10', 'BHN11', 'BHTRO', 'BHW01',
#             'BHW02', 'BHW03', 'BHW04', 'BHW05', 'BHW06', 'BHW07', 'BHW08', 'BHW09', 'BHW10', 'BHW11', 'BHW12',
#              'BHW13', 'BHW14', 'BHW15', 'BHW16', 'BHW17', 'BHW18']

# BPE
# stations = ['BUMT', 'CHUK', 'DOCH', 'PARO', 'TASH']

# BSN
# stations = ['MONG', 'TRON', 'WANG']

# IC
# stations = ['LSA']

# IN
# stations = ['SHL']

#PIRE
# stations = ['JAFL', 'JAML']

# NQ
stations = ['KATNP', 'KNSET', 'KTNP2']

# NK
# stations = ['KKN']

# EvK2
# stations = ['EVN']

##########
# Hi-CLIMB
# stations=[  'H0010', 'H0020', 'H0030', 'H0040', 'H0050', 'H0060',
#             'H0070', 'H0080', 'H0090', 'H0100', 'H0120', 'H0130', 'H0150', 'H0160', 'H0170', 'H0180', 'H0190', 'H0200',
#              'H0210', 'H0220', 'H0230', 'H0240', 'H0250', 'H0260', 'H0270', 'H0280', 'H0290', 'H0310',  'H0330',
#              'H0340', 'H0350', 'H0360', 'H0370', 'H0380', 'H0390']

# stations=['H0400', 'H0410', 'H0420', 'H0430', 'H0440',
#              'H0460', 'H0470', 'H0480', 'H0490', 'H0500', 'H0510', 'H0520', 'H0530', 'H0540', 'H0550', 'H0560', 'H0570',
#              'H0580', 'H0590', 'H0600', 'H0610', 'H0620', 'H0630', 'H0640', 'H0641', 'H0650', 'H0655', 'H0660', 'H0670',
#             'H0680', 'H0690', 'H0700', 'H0710', 'H0720', 'H0730']

# stations=['H0740', 'H0750', 'H0760', 'H0770', 'H0780', 'H0790',
#              'H0800', 'H0810', 'H1000', 'H1010', 'H1020', 'H1030', 'H1040', 'H1050', 'H1060', 'H1070', 'H1071', 'H1080',
#             'H1090', 'H1100', 'H1110', 'H1120', 'H1130', 'H1140', 'H1150', 'H1160', 'H1170', 'H1180', 'H1190', 'H1200',
#              'H1210', 'H1220', 'H1230', 'H1240', 'H1250']

# stations=['H1260', 'H1270', 'H1280', 'H1290', 'H1300', 'H1310', 'H1320',
#              'H1330', 'H1340', 'H1350', 'H1360', 'H1370', 'H1380', 'H1390', 'H1400', 'H1405', 'H1410', 'H1415', 'H1420',
#              'H1421', 'H1422', 'H1423', 'H1425', 'H1430', 'H1440', 'H1450', 'H1460', 'H1470', 'H1480', 'H1490', 'H1500',
#              'H1510', 'H1520', 'H1530', 'H1540']

# stations= ['H1550', 'H1560', 'H1570', 'H1580', 'H1590', 'H1600', 'H1610', 'H1620',
#              'H1630', 'NBENS', 'NBIRA', 'NBUNG', 'NDOML', 'NG010', 'NG020', 'NG030', 'NG040', 'NG050', 'NG060', 'NGUMB',
#              'NHILE', 'NJANA', 'NKDMG', 'NNAMC', 'NP010', 'NP020', 'NP030', 'NP035', 'NP040', 'NP050', 'NP060', 'NP070',
#              'NP071', 'NP075', 'NP080']
# stations= ['NP082', 'NP085', 'NP090', 'NP100', 'NPHAP', 'NRUMJ', 'NSIND', 'NSUKT', 'NTHAK',
#              'NTUML']
##########


font = {'family': 'normal',
        'weight': 'normal',
        'size': 32}
matplotlib.rc('font', **font)

fig = plt.figure(num=1, dpi=100, facecolor='w', edgecolor='b')
fig.suptitle('Data Continuity', fontsize=50)
fig.set_size_inches(40, 60)

# set min and max dates
mindate = mpdates.date2num(dt.datetime.strptime('20120101', '%Y%m%d'))
maxdate = mpdates.date2num(dt.datetime.strptime('20200401', '%Y%m%d'))

x = range(int(mindate), int(maxdate))
xdates = [mpdates.num2date(xd) for xd in x]
kstas = np.zeros(len(xdates))

plotno = 1
for sta in stations:
    y = []
    for i in range(len(xdates)):
        if glob.glob(paths[1] +'/'+xdates[i].strftime('Y%Y/R%j.01')+'/''*.'+sta+'*'):
            y += [1]
            kstas[i] += 1
        elif glob.glob(paths[4] + '/' + xdates[i].strftime('Y%Y/R%j.01') + '/'+sta+'*'):
            y += [1]
            kstas[i] += 1
        elif glob.glob(paths[2] + '/' + xdates[i].strftime('Y%Y/R%j.01') + '/'+'*'+sta+'*'):
            y += [1]
            kstas[i] += 1
        elif glob.glob(paths[-1] + '/' + xdates[i].strftime('Y%Y/R%j.01') + '/'+'*'+sta+'*'):
            y += [1]
            kstas[i] += 1
        else:
            y += [0]
    ax = fig.add_subplot(len(stations), 1, plotno)
    ax.scatter(xdates, y, c='k')

    ax.set_ylim([0.5, 1.5])
    ax.set_xlim([min(x), max(x)])
    ax.set_ylabel(sta, rotation=0, fontsize=25)
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
# plt.show()


fig.savefig("cont_plot.png", dpi=100, facecolor='w', edgecolor='w',
            orientation='landscape', format='png',
            pad_inches=0.1)

#plt.plot(xdates,kstas,'b')
#plt.show()
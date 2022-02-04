"""
Get data from data centers without using mass-downloader.
This is only recommended for smaller chunks of data.

Author: CJC
Modified by K. Michailos
"""


def get_data(network, station, starttime, endtime, outdir='.'):
    """
    Function to download all data for a given network.

    Will download all stations available between start and end times given. \
    It will create day-long files, and will download data as day-long segments.

    :type network: str
    :param network: Network code
    :type starttime: UTCDateTime
    :param starttime: Time to begin donloading from, will use the date
    :type endtime: UTCDateTime
    :param endtime: Time to end download, will use the date
    :type outdir: str
    :param outdir: Path to write to, will write in Y????/R???.01 directories \
        within this head directory.

    .. note:: This function is slow and doesn't cope with errors, suggest \
        using the mass-downloader.
    """
    import obspy
    import os
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
    else:
        from obspy.fdsn import Client
    from obspy import UTCDateTime, Stream, read
    starttime = UTCDateTime(starttime)
    endtime = UTCDateTime(endtime)

    kdays = (endtime.datetime - starttime.datetime).days + 1
    for i in range(kdays):
        t1 = starttime + (86400 * i)
        t2 = t1 + 86400
        bulk_info = [(network, station, '01', '*', t1, t2)]

        client = Client('IRIS', debug=True)
        try:
            st = client.get_waveforms_bulk(bulk_info)
            # st.plot()
            # Now you should split the data into traces by station and channel and
            # Save them into Y????/R???.01 directories, named something useful,
            # See the SAMBA_archive for details.
            for tr in st:
                chan_seq = (tr.stats.station, tr.stats.network,
                              tr.stats.location, tr.stats.channel,
                              tr.stats.starttime.strftime('%Y.%j'))
                f_name = '.'.join(chan_seq)
                path = outdir + tr.stats.starttime.strftime('/Y%Y/R%j.01/')
                if not os.path.isdir(path):
                    os.makedirs(path)
                tr.write(path + f_name, format='MSEED')
        except Exception as e:
            print(e)
            print('No data for %s day, carry on.' % str(i))


# massdownloader doesn't work for some reason...
# Wrote this sort script to download data using get_wav instead....
starttime = '2019-01-01T00:00:00.0'
endtime = '2020-04-01T00:00:00.0'
network = 'NQ'
outdir = '/media/kmichall/SEISMIC_DATA/Data_archive/NQ/'
sta_file = ['/home/kmichall/Desktop/Codes/bitbucket/him_seismicity/maps/files/sta_NQ.txt']

sta_names = []
with open(sta_file[0], 'r') as f:
   for line in f:
       if line.startswith('#'):
           print(line)
           continue
       else:
           ln = line.split()
           print(ln)
           sta_names.append(ln[1])
# Running it at the moment only for the first 3 stations
for sta in sta_names[1:3]:
   get_data(network, sta, starttime, endtime, outdir)



#################################################
starttime = '2014-11-01T00:00:00.0'
endtime = '2018-01-01T00:00:00.0'
network = 'ZO'
outdir = '/media/kmichall/SEISMIC_DATA/Data_archive/HiKNet/'
sta = 'BJ01'

get_data(network, sta, starttime, endtime, outdir)


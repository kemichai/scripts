"""
Short script to detect earthquakes in continuous waveform data using
energy based triggering methods. Reads single channel day long
waveform data and using the recursive_sta_lta triggering routine
from ObsPy and returns a list of  network coincidence triggers.
Potential earthquake detections are stored in mseed files.

=============================================
Requirements:
    * ObsPy >= 1.1.0
    * numpy
    * datetime
=============================================

UNIL
:Date: Feb 2020
:Author: K. Michailos
"""

import obspy
import glob
import numpy as np
from obspy import read
import matplotlib.dates as mpdates
import datetime as dt
from obspy import UTCDateTime, Stream
from obspy.signal.trigger import coincidence_trigger

# Define paths where data are stored
data_paths = ['/media/kmichall/SEISMIC_DATA/Data_archive/GANSSER',
         '/media/kmichall/SEISMIC_DATA/Data_archive/HIMNT',
         '/media/kmichall/SEISMIC_DATA/Data_archive/GIC_XE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/Hi_CLIMB',
         '/media/kmichall/SEISMIC_DATA/Data_archive/BPE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/BSN',
         '/media/kmichall/SEISMIC_DATA/Data_archive/IC',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NSNI',
         '/media/kmichall/SEISMIC_DATA/Data_archive/PIRE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NAMASTE',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NK_KKN',
         '/media/kmichall/SEISMIC_DATA/Data_archive/EVK2',
         '/media/kmichall/SEISMIC_DATA/Data_archive/NQ']

# Define path to output waveforms of detections
outdir = '/home/kmichall/Desktop/outputs'


def read_data_from_wav_paths(date, waveform_paths):
    """
    Function to read all available data on the VERTICAL(*Z!!!) components of seismic sites
    for a given day.

    :type date: UTC datetime
    :param date: Day for which seismic data will be read
    :type waveform_paths: list of strings
    :param waveform_paths: paths showing the folders where data are stored

    :returns: Waveform data for given station/channels and days
    :rtype: Stream object
    """

    st_ = Stream()
    for i in range(len(waveform_paths)):
        if waveform_paths[i].split('/')[-1] == 'IC':
            try:
                st_ += read(waveform_paths[i] + '/' + UTCDateTime(date).strftime('Y%Y/R%j.01') + '/*E.*')
            except Exception as e:
                print(e)
                continue
        else:
            try:
                st_ += read(waveform_paths[i] + '/' + UTCDateTime(date).strftime('Y%Y/R%j.01') + '/*Z.*')
            except Exception as e:
                print(e)
                continue

    for tr in st_:
        if tr.stats.network == 'YL' and tr.stats.sampling_rate != 40.0:
            tr.resample(40.0)
    for tr in st_:
        tr.filter('bandpass', freqmin=5.0, freqmax=25.0)

    try:
        st_.merge()
        # st_ = st.split().detrend('simple').merge(fill_value=0)
    except Exception as e:
        print(e)
    return st_


# Define min and max dates for which to run the trigger routine
mindate = mpdates.date2num(dt.datetime.strptime('20011119', '%Y%m%d'))
maxdate = mpdates.date2num(dt.datetime.strptime('20011120', '%Y%m%d'))
x = range(int(mindate), int(maxdate))
days = [mpdates.num2date(xd) for xd in x]

# Create a list of lists of daily coincidence triggers
coinc_triggers = []
for day in days:
    print(str(day.year) + '-' + str(day.month) + '-' + str(day.day))
    st = read_data_from_wav_paths(day, data_paths)
    # Run coincidence trigger for each day
    day_coinc_triggers = coincidence_trigger(trigger_type="recstalta", thr_on=3.5,
                                             thr_off=3.0, stream=st,
                                             sta=0.3, lta=10, thr_coincidence_sum=5,
                                             details=True, delete_long_trigger=False)
    # Don't write empty list if no triggers on that day
    if day_coinc_triggers != []:
        print('Adding triggers...')
        coinc_triggers.append(day_coinc_triggers)
# Now that the daily coincidence triggers are created
# we go through them and create mseed files of all the
# available waveforms around the trigger time
for day_trig in coinc_triggers:
    day = day_trig[-1]['time']
    all_wav_files = []
    for j in range(len(data_paths)):
        all_wav_files += glob.glob(data_paths[j]+'/'+day.strftime('Y%Y/R%j.01') +
                                   '/*' + day.strftime('%Y.%j'))
    all_wav_files.sort()

    for trig in day_trig:
        print('Writing trigger : ' + str(trig['time']))
        trig_time = trig['time']
        st1 = Stream()
        for wav in all_wav_files:
            st1 += read(wav, starttime=trig_time - 10,
                        endtime=trig_time + 20)
            # Write mseed file
            st1.write("%s/%s.MSEED" % (outdir,
                      str(trig_time - 10).split('.')[0]),
                      format='MSEED')
        # Plot all the waveforms around the trigger for a
        # quick inspection (comment line if not needed)
        # st1.plot(equal_scale=False, size=(800, 600))

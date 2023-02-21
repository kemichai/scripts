"""
:author: Konstantinos Michailos
:date: 13/09/2018

"""
from obspy import read_events
import numpy as np


def read_cat4vars(cat_path, output_filename='parameters.npz'):
    """
    Function to read *.xml earthquake catalogs and create lists of the eq parameters (e.g. Lat, Lon, Dep, etc.)
    and save these values in a dictionary.

    :param cat_path: Path to the quakeml file.
    :type cat_path: str
    """

    print('Reading earthquake catalog ' + cat_path.split('/')[-1] + '...')
    cat = read_events(cat_path)
    print('Done reading earthquake catalog')

    # Earthquake location details
    event_id = []
    dattime = []
    lat = []
    lon = []
    dep = []
    mag= []
    for ev in cat:
        # Earthquake location details
        event_id.append(ev.preferred_origin_id)
        o = ev.preferred_origin()
        lat.append(o.latitude)
        lon.append(o.longitude)
        dep.append(o.depth)
        mag.append(float(ev.magnitudes[0].mag))
        dattime.append(o.time)
    print('Writing dictionary ...')
    # Write a dictionary
    dictionary = {"eq_ids": event_id, "longitudes": lon, "latitudes": lat,
                  "depths": dep, "magnitudes": mag, "datetimes":dattime}
    np.savez(output_filename, **dictionary)

    return


# Create the npz file
a = read_cat4vars('/home/kmichailos/Desktop/codes/github/codes/PhD_codes/Focal_mechanisms_stress/quakeml_file/subset_cat.xml')

# Read the npz file
parameters = np.load('parameters.npz')
parameters.items()
lats = parameters["latitudes"]


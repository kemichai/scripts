"""
:author: Konstantinos Michailos
:date: 13/09/2018

"""
from obspy import read_events
from progressbar import ProgressBar
import numpy as np


def read_cat4vars(cat_path):
    """
    Function to read *.xml earthquake catalogs and create lists of the eq parameters (e.g. Lat, Lon, Dep, etc.)
    and save these values in a dictionary.

    :param cat_path: Path to the quakeml file.
    :type cat_path: str

    Example:
    To read the dictionary use:
    >>> parameters = np.load("parameters.npz")
    >>> parameters.items()
    # This will output the names of the lists...
    >>> parameters.keys()
    >>> lats = parameters["latitude"]
    # etc

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




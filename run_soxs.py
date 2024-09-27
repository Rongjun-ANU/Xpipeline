#!/usr/bin/env python3

# Import necessary libraries
import sys
from IPython.display import set_matplotlib_formats
from scipy import interpolate
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from decimal import *
# from tabulate import tabulate
from astropy.table import QTable
from astropy.io import fits
from astropy.io import ascii
from astropy.stats import sigma_clip
import math
import numpy as np
np.float = np.float64
np.int = np.int64
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import scipy.stats as st
from astropy import constants as c
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
# %matplotlib inline
# set_matplotlib_formats('svg')
import matplotlib.cm
import glob
import h5py
import os
import time as ostime
import scipy.interpolate as interpolate
from multiprocessing import Process, Queue
import multiprocessing
import matplotlib.cm     as cm
from scipy import optimize
from scipy import stats
import datetime
import argparse
import unyt
from unyt import unyt_array
import yt
# yt.toggle_interactivity()
from yt.units import kpc
from yt.units import dimensions
import run_pyxsim
from pyxsim.source_models import CIESourceModel
import soxs


def simulate_instrument(events_file, tag, exp_time, exp_label, dist):
    # Determine filenames based on the events file
    base_filename = events_file.replace('_events.h5', '')  # Remove '_events.h5' for base
    evt_filename = f"{base_filename}_evt_{exp_label}.fits"
    img_filename = f"{base_filename}_img_{exp_label}.fits"
    spec_filename = f"{base_filename}_spec_{exp_label}.pha"

    # Simulate the instrument's response
    soxs.instrument_simulator(
        events_file,
        evt_filename,
        exp_time,
        "chandra_acisi_cy22",  # Replace with "chandra_aciss_cy22" if needed
        [30., 45.],
        overwrite=True,
        ptsrc_bkgnd=False,
        instr_bkgnd=False,
        foreground=False
    )

    # Write the image
    soxs.write_image(evt_filename, img_filename, 
                     emin=0.4, emax=2.0, overwrite=True)

    # Write the spectrum
    soxs.write_spectrum(evt_filename, spec_filename, 
                        format='fits', emin=0.4, emax=2.0, overwrite=True)

    print(f"Processed {tag}: events {events_file} with exposure {exp_label} to EVT, IMG, and PHA.")

def main(events_files, exp_times, dist):
    for events_file in events_files:
        # Extract the tag from the file name for logging
        parts = events_file.split('_')
        tag = '_'.join(parts[:-2])  # Reconstruct tag without '_events.h5'
        for exp_time, exp_label in exp_times:
            simulate_instrument(events_file, tag, exp_time, exp_label, dist)

if __name__ == "__main__":
    # Example exposure times and distance - replace with actual values as needed
    exp_times = [
        (yt.YTQuantity(1000, "Ms"), "1000Ms"),
        (yt.YTQuantity(100, "Ms"), "100Ms"),
        (yt.YTQuantity(10, "Ms"), "10Ms"),
        (yt.YTQuantity(1, "Ms"), "1Ms"),
        (yt.YTQuantity(0.1, "Ms"), "0.1Ms")
    ]
    dist = yt.YTQuantity(500, "kpc")  # kpc
    events_files = sys.argv[1:]  # Expects file names like 'combined_*_events.h5'
    start_time = ostime.time()
    main(events_files, exp_times, dist)
    end_time = ostime.time()
    print(f"Elapsed time taken for running run_soxs_hires.py: {end_time - start_time} seconds")

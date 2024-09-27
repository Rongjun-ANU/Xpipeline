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

# plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [12, 8]

import concurrent.futures

# Define the main function
def main():
    start_time = ostime.time()
    # Check if the npz filename is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: source pyxsim.py <filename>.npz")
        sys.exit(1)

    # The first command-line argument is the script name, so we need the second one.
    npz_filename = sys.argv[1]

    # Extract the base name from the npz filename
    # Assumes the file is named like 'plt1555000_Ofive.npz' and you want 'plt1555000' as the base name
    base_name = os.path.basename(npz_filename).split('.')[0]

    # Constants
    proton_mass = yt.YTQuantity(1.6726219e-24, "g")  # Proton mass in grams
    solar_mass = c.M_sun.to(u.g).value * yt.units.g  # Ensure this line matches your context

    # Adjusting the dimensions
    dimension_x, dimension_y, dimension_z = 128*4, 128*4, 1024*4
    dx = dy = dz = yt.YTQuantity(2.3578125e+19/4, 'cm')

    # Read the data from the npz file
    data_from_npz = np.load(npz_filename)

    # Create a data dictionary using the data from the file
    data = {
        ('gas', 'density'): yt.YTArray(data_from_npz['Density'], "g/cm**3"),
        ('gas', 'temperature'): yt.YTArray(data_from_npz['Temperature'], "K"),
        ('gas', 'mass'): yt.YTArray(data_from_npz['Mass'], "g"),
        ('gas', 'emission_measure'): yt.YTArray(data_from_npz['Emission_Measure'], "cm**-3"),
        ('gas', 'Energy'): yt.YTArray(data_from_npz['GasEnergy'], "erg"),
        ('gas', 'InternalEnergy'): yt.YTArray(data_from_npz['GasInternalEnergy'], "erg/cm**3"),
        ('gas', 'scalar_0'): yt.YTArray(data_from_npz['Scalar_0'], "dimensionless"),
        ('gas', 'scalar_1'): yt.YTArray(data_from_npz['Scalar_1'], "dimensionless"),
        ('gas', 'scalar_2'): yt.YTArray(data_from_npz['Scalar_2'], "dimensionless"),
        ('gas', 'x_GasMomentum'): yt.YTArray(data_from_npz['x_GasMomentum'], "g*cm/s"),
        ('gas', 'y_GasMomentum'): yt.YTArray(data_from_npz['y_GasMomentum'], "g*cm/s"),
        ('gas', 'z_GasMomentum'): yt.YTArray(data_from_npz['z_GasMomentum'], "g*cm/s"),
        ('gas', 'velocity_x'): yt.YTArray(np.zeros((dimension_x, dimension_y, dimension_z)), "cm/s"),
        ('gas', 'velocity_y'): yt.YTArray(np.zeros((dimension_x, dimension_y, dimension_z)), "cm/s"),
        ('gas', 'velocity_z'): yt.YTArray(np.zeros((dimension_x, dimension_y, dimension_z)), "cm/s"),
        ('gas', 'metallicity'): yt.YTArray(data_from_npz['Zmet'], "dimensionless")
    }

    # Adjusting the bounding box to accommodate the new dimensions
    bbox = np.array([[0, dimension_x * dx], 
                     [0, dimension_y * dy], 
                     [-dimension_z * dz / 2, dimension_z * dz / 2]])

    # Create a yt dataset
    ds = yt.load_uniform_grid(data, (dimension_x, dimension_y, dimension_z), 
                              length_unit="cm", 
                              bbox=bbox, 
                              nprocs=128, default_species_fields="ionized")

    # Further analysis and visualization with yt can go here

    print(f"Dataset from {npz_filename} loaded successfully into yt.")

    # Constants
    dimension_z = 1024*4
    dz = ds.domain_width[2] / dimension_z  # Assuming ds.domain_width[2] gives z width

    # Center slice 'D' modifications for 64*dz
    z_center = ds.domain_center[2]
    D_start = z_center - 32*4 * dz  # Central slice starts 32*dz before the center
    D_end = z_center + 32*4 * dz  # Central slice ends 32*dz after the center
    D = ds.region(center=[dimension_x*dx/2, dimension_y*dy/2, z_center],
                left_edge=[0, 0, D_start],
                right_edge=[dimension_x*dx, dimension_y*dy, D_end])
    print(f"Finish slicing of Region D in {base_name}")

    # Above slices 'A1' to 'A5', adjustments for 96*dz
    for i in range(1, 6):
        locals()[f'A{i}'] = ds.region(center=[dimension_x*dx/2, dimension_y*dy/2, D_end + (48*4 + 96*4 * (i - 1)) * dz],
                                    left_edge=[0, 0, D_end + 96*4 * (i - 1) * dz],
                                    right_edge=[dimension_x*dx, dimension_y*dy, D_end + 96*4 * i * dz])
        print(f"Finish slicing of Region A{i} in {base_name}")
        # print(locals()[f'A{i}'].left_edge[2], locals()[f'A{i}'].right_edge[2])
    # A5.right_edge[2] = ds.domain_right_edge[2]

    # Below slices 'B1' to 'B5', adjustments for 96*dz
    for i in range(1, 6):
        locals()[f'B{i}'] = ds.region(center=[dimension_x*dx/2, dimension_y*dy/2, D_start - (48*4 + 96*4 * (i - 1)) * dz],
                                    left_edge=[0, 0, D_start - 96*4 * i * dz],
                                    right_edge=[dimension_x*dx, dimension_y*dy, D_start - 96*4 * (i - 1) * dz])
        print(f"Finish slicing of Region B{i} in {base_name}")
    # Adjust the left edge of B5 to be ds.domain_left_edge in the z-direction
    # B5.left_edge[2] = ds.domain_left_edge[2]


    # Load all_data() into a variable for further analysis
    ad = ds.all_data()

    
    # Emission model for X-ray photons
    T_min = np.min(ad['gas', 'temperature'])  # minimum temperature for which we'll emit X-rays
    T_max = np.max(ad['gas', 'temperature'])  # maximum temperature for which we'll emit X-rays
    e_min = (c.k_B*T_min*u.K).to("keV")
    e_max = (c.k_B*T_max*u.K).to("keV")

    # source_model = pyxsim.CIESourceModel("spex", e_min, e_max, nbins=1024, Zmet=1.0, binscale="log")
    source_model = run_pyxsim.CIESourceModel("apec", 0.4, 2.0, 
                                        nbins=256, 
                                        Zmet=("gas", "metallicity"), 
    #                                      Zmet=1, 
                                        binscale="log") # hot gas
    

    # Key parameters here!!! 
    exp_time = yt.YTQuantity(1000, "Ms")  # exposure time
    area = yt.YTQuantity(np.pi*0.6**2, "m**2")  # collecting area yt.YTQuantity(np.pi*0.6**2, "m**2")
    # redshift = 0.0001 # 0.4Mpc
    dist=yt.YTQuantity(500, "kpc") # kpc
    print(f"Key parameters in {base_name} are Exposure time: {exp_time}, Collecting area: {area}, Redshift: 0, Distance: {dist}")


    # Make photons

    # For entire region
    entire_region = ds.all_data()
    n_photons, n_cells = pyxsim.make_photons(
        f"{base_name}_photons",  # Use the base_name variable here 
        entire_region, redshift=0, dist=dist, area=area, exp_time=exp_time, source_model=source_model
    )
    print(f"Finish making photons of entire region in {base_name}")

    # For central region 'D'
    n_photons_D, n_cells_D = run_pyxsim.make_photons(
        f"{base_name}_D_photons",  # Use the base_name variable here 
        D, redshift=0, dist=dist, area=area, exp_time=exp_time, source_model=source_model
    )
    print(f"Finish making photons of D in {base_name}")

    # For above slices 'A1' to 'A5'
    for i in range(1, 6):
        locals()[f'n_photons_A{i}'], locals()[f'n_cells_A{i}'] = run_pyxsim.make_photons(
            f"{base_name}_A{i}_photons", 
            locals()[f'A{i}'], redshift=0, dist=dist, area=area, exp_time=exp_time, source_model=source_model
        )
        print(f"Finish making photons of A{i} in {base_name}")

    # For below slices 'B1' to 'B5'
    for i in range(1, 6):
        locals()[f'n_photons_B{i}'], locals()[f'n_cells_B{i}'] = run_pyxsim.make_photons(
            f"{base_name}_B{i}_photons", 
            locals()[f'B{i}'], redshift=0, dist=dist, area=area, exp_time=exp_time, source_model=source_model
        )
        print(f"Finish making photons of B{i} in {base_name}")

    # Project photons
        
    # For entire region 
    n_events = pyxsim.project_photons(
        f"{base_name}_photons",
        f"{base_name}_events",
        "y",
        [30., 45.],
        absorb_model=None,
        nH=None
    )
    print(f"Finish projecting photons of entire region in {base_name}")

    # For central region 'D'
    n_events_D = run_pyxsim.project_photons(
        f"{base_name}_D_photons",
        f"{base_name}_D_events",
        "y",
        [30., 45.],
        absorb_model=None,
        nH=None
    )
    print(f"Finish projecting photons of D in {base_name}")

    # For above slices 'A1' to 'A5'
    for i in range(1, 6):
        locals()[f'n_events_A{i}'] = run_pyxsim.project_photons(
            f"{base_name}_A{i}_photons",
            f"{base_name}_A{i}_events",
            "y",
            [30., 45.],
            absorb_model=None,
            nH=None
        )
        print(f"Finish projecting photons of A{i} in {base_name}")

    # For below slices 'B1' to 'B5'
    for i in range(1, 6):
        locals()[f'n_events_B{i}'] = run_pyxsim.project_photons(
            f"{base_name}_B{i}_photons",
            f"{base_name}_B{i}_events",
            "y",
            [30., 45.],
            absorb_model=None,
            nH=None
        )
        print(f"Finish projecting photons of B{i} in {base_name}")


####################################################################################################################                                
    # Only need simput files if SOXS version is before 4.3.0
    
    # Write the photon lists to simput files
    
    # # Write the photon lists to simput files
    # # For central region 'D'
    # events_D = pyxsim.EventList(f"{base_name}_D_events.h5")
    # events_D.write_to_simput(f"{base_name}_D", overwrite=True)
    # print(f"Finish creating photon list of D in {base_name}")

    # # For above slices 'A1' to 'A5'
    # for i in range(1, 6):
    #     locals()[f'events_A{i}'] = pyxsim.EventList(f"{base_name}_A{i}_events.h5")
    #     locals()[f'events_A{i}'].write_to_simput(f"{base_name}_A{i}", overwrite=True)
    #     print(f"Finish creating photon list of A{i} in {base_name}")

    # # For below slices 'B1' to 'B5'
    # for i in range(1, 6):
    #     locals()[f'events_B{i}'] = pyxsim.EventList(f"{base_name}_B{i}_events.h5")
    #     locals()[f'events_B{i}'].write_to_simput(f"{base_name}_B{i}", overwrite=True)
    #     print(f"Finish creating photon list of B{i} in {base_name}")

    # # For entire region
    # events = pyxsim.EventList(f"{base_name}_events.h5")
    # events.write_to_simput(f"{base_name}", overwrite=True)
    # print(f"Finish creating photon list of entire region in {base_name}")

    # # print("Simput files created successfully.")

####################################################################################################################                                
        
    end_time = ostime.time()
    print(f"Elapsed time taken of running pyxsim_hires.py for of {base_name}: {end_time - start_time} seconds")

if __name__ == "__main__":
    main()

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

start_time = ostime.time()

# Import values
gamma = 5/3
boltzmann_constant_cgs = (c.k_B.to(u.erg / u.K)).value*unyt.erg/unyt.K
m_p_cgs = (c.m_p.to(u.g)).value*unyt.g

def main():
    # Check if the dataset name is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: source read.py <dataset_name>")
        sys.exit(1)

    # The first command-line argument is the script name, so we need the second one.
    dataset_name = sys.argv[1]

    # Load the dataset using yt
    ds = yt.load(dataset_name, default_species_fields='ionized')

    # Add any additional processing you want here
    print(f"Dataset {dataset_name} loaded successfully.")
    
        # Define the level of refinement. For base level, use level=0.
    level = 0

    # Find the domain dimensions of the dataset at the specified level of refinement.
    dims = ds.domain_dimensions * ds.refine_by**level

    # Use the covering_grid method with the dataset's domain left and right edges.
    covering_grid = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)
    gasDensity = covering_grid[("boxlib", "gasDensity")]
    gasInternalEnergy = covering_grid[("boxlib", "gasInternalEnergy")]

    time = ds.current_time.to('Myr')

    def _gasDensity(field, data):
        density = gasDensity*unyt.g/unyt.cm**3
        return density

    ds.add_field(('gas', 'density'), 
                            function=_gasDensity, units='g/cm**3', sampling_type='cell', force_override=True)

    def _gasInternalEnergy(field, data):
        InternalEnergy = gasInternalEnergy*unyt.erg
        return InternalEnergy

    ds.add_field(('gas', 'gasInternalEnergy'), 
                            function=_gasInternalEnergy, units='erg', sampling_type='cell', force_override=True)


    # Define the temperature field
    print(f"Start to calculate the temperature field of {dataset_name}.")

    #Grackle Cooling Table
    file = 'grackle_data_files/input/CloudyData_UVB=HM2012.h5'
    grackle = h5py.File(file)
    array = grackle['CoolingRates/Primordial/MMW'][()]
    table = array[:,0,:]
    table_nH   = np.logspace(-6, 4, array.shape[0])
    table_temp = np.logspace(1,  9, array.shape[2])

    bins = 100
    egas_arr = np.logspace(-16., -5., bins)
    nH_arr   = np.logspace(-6.0, 4.0, int(bins))
    logrho_arr = np.log10(nH_arr[:-1])
    logEgas_arr = np.log10(egas_arr[:-1])

    #Set up the interpolator

    T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))
    i=0
    for egas in egas_arr:
        j=0
        for nH in nH_arr:
            C = (gamma - 1.) * egas / (boltzmann_constant_cgs.v*nH)
            minT = C*np.amin(table)
            maxT = C*np.amax(table)
            def func(T):
                mu = interpolate.interp2d(table_temp, table_nH, table,\
                                kind='linear', copy=True, bounds_error=False, fill_value=None)
                return C*mu(T,nH)[0] - T

            T[i,j] = optimize.toms748(func, minT, maxT)
            j+=1
        i+=1

    # import sys
    # import os
    # import argparse

    # # Mock the command line arguments
    # sys.argv = ['ipykernel_launcher.py', '--file', os.getcwd()]


    lev = 0
    home = os.getcwd()

    # parser = argparse.ArgumentParser(description='Optional app description')
    # parser.add_argument('--file', type=str, help='Filename to be analysed')
    # args = parser.parse_args()

    class Data:
        fac = 1
        lev = 0
        file = ''
        dom_min = [0.0, 0.0, 0.0]
        def getData(file):
            dom_min = [0.0, 0.0, 0.0]
            ds   = yt.load(file)
            data = ds.covering_grid(level=lev, left_edge=dom_min, dims=ds.domain_dimensions * fac)
            density = np.array(data['gasDensity'])
            time = ds.current_time.to('Myr')
            Egas = np.array(data["gasInternalEnergy"])
            return density, Egas, time

    def getdomain(file):
        infile = open(file)
        lines = infile.readlines()
        dom_range = np.zeros((2,3))
        ncell = np.zeros(3)
        dom_min = [0.0,0.0,0.0]
        dom_min[0] = float(lines[3].split()[2])
        dom_min[1] = float(lines[3].split()[3])
        dom_min[2] = float(lines[3].split()[4])
        
        dom_max = [0.0,0.0,0.0]
        dom_max[0] = float(lines[4].split()[2])
        dom_max[1] = float(lines[4].split()[3])
        dom_max[2] = float(lines[4].split()[4])

        ncell[0]=int(lines[15].split()[2])
        ncell[1]=int(lines[15].split()[3])
        ncell[2]=int(lines[15].split()[4])
        
        return dom_min, dom_max, ncell

    # Extracted information from the log
    domain_left_edge = ds.domain_left_edge.to('kpc')
    domain_right_edge = ds.domain_right_edge.to('kpc')
    domain_dimensions = ds.domain_dimensions
    domain_center = ds.domain_center.to('kpc')

    # Compute ranges
    x_range = (domain_left_edge[0], domain_right_edge[0])
    y_range = (domain_left_edge[1], domain_right_edge[1])
    z_range = (domain_left_edge[2], domain_right_edge[2])

    # Extract number of cells
    n_cell_x = domain_dimensions[0]
    n_cell_y = domain_dimensions[1]
    n_cell_z = domain_dimensions[2]
    n_cell_total = n_cell_x * n_cell_y * n_cell_z

    timestep = ds.current_time.to('Myr')

    rho = gasDensity
    egas0 = gasInternalEnergy

    rho0 = rho/m_p_cgs

    logrho_arr = np.log10(nH_arr[:-1])
    logrho = np.log10(rho0)
    delta_rho = logrho_arr[1] - logrho_arr[0]
    idxrho = (np.floor((logrho - np.amin(logrho_arr))/delta_rho)).astype('int')

    logEgas_arr = np.log10(egas_arr[:-1])
    logEgas = np.log10(egas0)
    delta_egas = logEgas_arr[1] - logEgas_arr[0]
    idxegas = (np.floor((logEgas-np.amin(logEgas_arr))/delta_egas)).astype('int')


    wgt_rho  = (logrho - (np.amin(logrho_arr) + delta_rho*idxrho))/delta_rho
    wgt_egas = (logEgas - (np.amin(logEgas_arr) + delta_egas*idxegas))/delta_egas


    temp = (1.-wgt_rho)*(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
            wgt_rho *    wgt_egas * T[tuple(idxegas+1), tuple(idxrho+1)] +\
        (1. -wgt_rho)*    wgt_egas * T[tuple(idxegas+1), tuple(idxrho)]   +\
            wgt_rho *(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)] 

    def _temp_field(field, data):
    #     reshaped_temp = temp.reshape(ds.domain_dimensions[::1])
        reshaped_temp = temp
        return data.ds.arr(reshaped_temp, "K")

    ds.add_field(("gas", "temperature"), function=_temp_field, units="K", sampling_type='cell', 
                            force_override=True)

    ad_plt1555000 = ds.all_data()
    print(f"Finish calculating the temperature field of {dataset_name}.")


    # Now do emssion measure
    print(f"Start to calculate the emission measure field of {dataset_name}.")
    # Calculate the size of the domain in each dimension
    domain_size = domain_right_edge - domain_left_edge

    # Calculate the cell width in each dimension
    cell_widths = domain_size / domain_dimensions

    # Assuming cubic cells, take the cell width from one dimension
    dx = (cell_widths[0]).to(unyt.cm)
    dy = (cell_widths[1]).to(unyt.cm)
    dz = (cell_widths[2]).to(unyt.cm)

    def _gasMass(field, data):
        mass = ad_plt1555000[('gas', 'density')]*dx*dy*dz
        return mass
    ds.add_field(('gas', 'mass'), 
                            function=_gasMass, units='g', sampling_type='cell', force_override=True)
    # ad_plt1555000[('gas', 'mass')]

    ad_plt1555000 = ds.all_data()
    def _emission_measure(field, data):
        nH = ad_plt1555000[('gas', 'density')]/m_p_cgs
    #     ne = 1.2 * nH
        return nH **2 * dx * dy * dz  # .d extracts the value as a float

    ds.add_field(('gas', 'emission_measure'), 
                            function=_emission_measure, units='cm**-3', sampling_type='cell', force_override=True)
    print(f"Finish calculating the emission measure field of {dataset_name}.")


    # Now extract the data 
    # Solar mass
    solar_mass = c.M_sun.to(u.g).value * unyt.g

    # Create a covering grid for the entire domain at the finest resolution level
    level = 0  # Change this if you need a different refinement level
    dims = ds.domain_dimensions * 2**level
    covering_grid = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)

    # Extract the properties from the dataset without flattening
    density_data = covering_grid[('gas', 'density')].d
    temperature_data = covering_grid[('gas', 'temperature')].d
    mass_data = covering_grid[('gas', 'mass')].d
    emission_measure_data = covering_grid[('gas', 'emission_measure')].d

    gasEnergy_data = covering_grid[('boxlib', 'gasEnergy')].d
    gasInternalEnergy_data = covering_grid[('boxlib', 'gasInternalEnergy')].d
    scalar_0_data = covering_grid[('boxlib', 'scalar_0')].d
    scalar_1_data = covering_grid[('boxlib', 'scalar_1')].d
    scalar_2_data = covering_grid[('boxlib', 'scalar_2')].d
    x_GasMomentum_data = covering_grid[('boxlib', 'x-GasMomentum')].d
    y_GasMomentum_data = covering_grid[('boxlib', 'y-GasMomentum')].d
    z_GasMomentum_data = covering_grid[('boxlib', 'z-GasMomentum')].d

    # Extract the base name from the dataset name
    base_name = os.path.basename(dataset_name)  # This removes the directory path
    base_name = os.path.splitext(base_name)[0]  # This removes the file extension, if it exists

    # Create a new filename for the output file

    print(f"Start to calculate metallicity field of {dataset_name} and extract the data.")

    from concurrent.futures import ThreadPoolExecutor

    # Function for calculating and saving metallicity data with specific file names
    def calculate_and_save_metallicity(suffix, Zmet_value, base_name, density_data, temperature_data, mass_data, emission_measure_data, gasEnergy_data, gasInternalEnergy_data, scalar_0_data, scalar_1_data, scalar_2_data, x_GasMomentum_data, y_GasMomentum_data, z_GasMomentum_data):
        # Calculate the metallicity
        Zmet_data = Zmet_value + scalar_2_data / density_data * (solar_mass.value) / 8
        # Construct file name
        file_name = f"{base_name}_{suffix}.npz"
        # Save the data arrays
        np.savez_compressed(file_name, 
                            Density=density_data, Temperature=temperature_data, Mass=mass_data, 
                            Emission_Measure=emission_measure_data, GasEnergy=gasEnergy_data, 
                            GasInternalEnergy=gasInternalEnergy_data, Scalar_0=scalar_0_data, 
                            Scalar_1=scalar_1_data, Scalar_2=scalar_2_data, 
                            x_GasMomentum=x_GasMomentum_data, y_GasMomentum=y_GasMomentum_data, 
                            z_GasMomentum=z_GasMomentum_data, Zmet=Zmet_data)
        print(f"Data written to {file_name}")

    # List of metallicity values and corresponding file suffixes
    tasks = [
        ("Zero", 0.0), # 0 background metallicity
        ("Ofive", 0.5), # 0.5 background metallicity
        ("Otwo", 0.2), # 0.2 background metallicity
        ("One", 1.0), # 1.0 background metallicity
        ("Two", 2.0), # 2.0 background metallicity
    ]

    # Using ThreadPoolExecutor to run calculations in parallel
    with ThreadPoolExecutor(max_workers=3) as executor:
        for suffix, Zmet_value in tasks:
            executor.submit(calculate_and_save_metallicity, suffix, Zmet_value, base_name, density_data, temperature_data, mass_data, emission_measure_data, gasEnergy_data, gasInternalEnergy_data, scalar_0_data, scalar_1_data, scalar_2_data, x_GasMomentum_data, y_GasMomentum_data, z_GasMomentum_data)

    # Force the uniform solar metallicity 
    Zmet_data = np.ones_like(density_data) # uniform solar metallicity
    # Save the data arrays using np.savez_compressed
    np.savez_compressed(f'{base_name}_Uni.npz', 
                        Density=density_data, Temperature=temperature_data, Mass=mass_data, 
                        Emission_Measure=emission_measure_data, GasEnergy=gasEnergy_data, 
                        GasInternalEnergy=gasInternalEnergy_data, Scalar_0=scalar_0_data, 
                        Scalar_1=scalar_1_data, Scalar_2=scalar_2_data, 
                        x_GasMomentum=x_GasMomentum_data, y_GasMomentum=y_GasMomentum_data, 
                        z_GasMomentum=z_GasMomentum_data, Zmet=Zmet_data)
    print(f"Data written to {base_name}_Uni.npz")

    end_time = ostime.time()
    print(f"Elapsed time taken of running read.py for of {dataset_name}: {end_time - start_time:.2f} seconds.")


if __name__ == "__main__":
    main()

# End of file
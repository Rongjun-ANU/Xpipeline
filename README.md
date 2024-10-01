## Xpipeline: A [QUOKKA](https://github.com/quokka-astro/quokka) Simulation Post-Processing and X-ray Analysis Pipeline 

**Xpipeline** is a Python-based pipeline designed for the post-processing and analysis of QUOKKA simulation data in the soft X-ray band (0.4-2.0 keV). The repository includes a set of scripts that read, process, and simulate synthetic X-ray emissions, allowing researchers to explore various physical properties and simulate instrument responses.

### Key Features
- **Data Ingestion & Processing**: Reads simulation data and computes derived fields such as gas density, temperature, and emission measure.
- **Synthetic Photon Simulation**: Creates synthetic X-ray photons using thermal sources in Collisional Ionization Equilbrium (CIE) and projects them into events files.
- **Instrument Response Simulation**: Simulates X-ray instrument responses (e.g., Chandra ACIS-I) to generate realistic event, image, and spectrum files for visualization.

### Workflow Overview
1. **`read.py`**: Reads the QUOKKA simulation dataset, calculates the missing properties, and saves the data into a yt object in a `.npz` file. 
   - **Temperature Field**: Interpolated using the internal energy and density values.
   - **Emission Measure**: Calculated as the square of the number density (approximation).
   - **Metallicity Variants**: Saves data files with varying background metallicities to explore their effect on X-ray emissions.

2. **`run_pyxsim.py`**: Uses the processed `.npz` files to simulate X-ray photons using the `pyxsim` library.
   - **Region Definitions**: Central regions (`D`), above slices (`A1` to `A5`), and below slices (`B1` to `B5`).
   - **Photon Creation**: Uses source models like APEC to simulate hot gas emissions.
   - **Projection**: Creates synthetic X-ray events for various regions.

3. **`run_soxs.py`**: Takes the photon event files and simulates instrument responses using `SOXS`.
   - **Multiple Exposure Times**: Simulates the response for various exposure times (e.g., 1 Ms, 100 Ms).
   - **Output Formats**: Generates FITS files for events, images, and spectra.

### Getting Started
1. **Installation**: Clone the repository and install the necessary dependencies:
   ```bash
   git clone https://github.com/your_username/Xpipeline.git
   cd Xpipeline
   pip install -r requirements.txt
   ```

2. **Running the Scripts**:
   - **`read.py`**: Reads and processes the dataset.
     ```bash
     python3 read.py <dataset_name>
     ```
   - **`run_pyxsim.py`**: Generates synthetic photon events.
     ```bash
     python3 run_pyxsim.py <processed_file>_<background_metallicity>.npz
     ```
   - **`run_soxs.py`**: Simulates instrument response for the photon events.
     ```bash
     python3 run_soxs.py <processed_file>_<background_metallicity>_events.h5
     ```

### Dependencies

The following Python packages are required to run the **Xpipeline**. Ensure that these packages are installed in your environment:

```bash
attrs==20.2.0
backcall==0.2.0
certifi==2020.6.20
cycler==0.10.0
Cython==0.29.21
decorator==4.4.2
hypothesis==5.33.1
iniconfig==0.0.0
ipython==7.18.1
ipython-genutils==0.2.0
jedi==0.17.2
kiwisolver==1.2.0
matplotlib==3.3.1
more-itertools==8.5.0
numpy==1.19.1
packaging==20.4
parso==0.7.1
pexpect==4.8.0
pickleshare==0.7.5
Pillow==7.2.0
pip==24.2
pluggy==0.13.1
prompt-toolkit==3.0.7
ptyprocess==0.6.0
py==1.9.0
Pygments==2.6.1
pyparsing==2.4.7
pytest==6.0.1
python-dateutil==2.8.1
scipy==1.5.2
setuptools==50.3.0
six==1.15.0
sortedcontainers==2.2.2
toml==0.10.1
traitlets==5.0.4
wcwidth==0.2.5
```

### Example Workflow
1. **Read the Simulation Data**:
   ```bash
   python3 read.py plt9640000
   ```
   Outputs multiple `.npz` files with varying metallicities. E.g., for half solar background, it will be `plt9640000_Ofive.npz`. 

2. **Simulate X-ray Photons**:
   ```bash
   python3 run_pyxsim.py plt9640000_Ofive.npz
   ```
   Generates synthetic X-ray photon events for different regions. Outputs multiple `*_event.h5` files as the projection photon lists.
   You may choose to do the slicing. E.g., `plt9640000_Ofive_event.h5` for the entire region, or `plt9640000_Ofive_D_event.h5` for the central disc region "D".

4. **Simulate Instrument Response**:
   ```bash
   python3 run_soxs.py plt9640000_Ofive_events.h5
   ```
   Produces realistic event, image, and spectrum files for the simulated observation.

### Warning

1. **Grackle Cooling Table Requirement**: 
   - When calculating the temperature field in `read.py`, it is crucial to use the [Grackle cooling table](https://github.com/grackle-project/grackle) for accurate interpolation. Ensure that the Grackle library is properly set up and the necessary data files are available (`CloudyData_UVB=HM2012.h5` or similar) to avoid any runtime errors.

2. **Memory Usage Considerations**:
   - Running the entire pipeline can be *extremely computationally expensive*, especially for large exposure times such as 1000Ms (~31 years). A full 1000Ms observation for the entire region or the central disk region may require up to **3TB** of memory during the `project_photon` step in `run_pyxsim.py` and subsequently in `run_soxs.py`. Even smaller exposure times, like 100Ms, may require around **500GB** of memory. 

   - **Suggested Approach**: Consider breaking `run_pyxsim.py` into two separate parts:
     1. **Photon Generation**: Run up to the `make_photons` step to generate the `*photons.h5` files.
     2. **Memory Estimation**: Check the size of the generated `*photons.h5` file. If the file size is `XXXGB`, expect the subsequent projection step (`project_photons`) to need approximately three times `XXXGB` memory (e.g., a `200GB` `plt9640000_Ofive_D_photons.h5` file may require `600GB` of memory to generate the projection photon list `plt9640000_Ofive_D_events.h5`).

   - If memory resources are a concern, consider reducing the simulation region or exposure time or collection area, or using higher performance computing facilities. You may be looking for [Intel Optane DC Persistent Memory](https://software.intel.com/content/www/us/en/develop/articles/quick-start-guide-configure-intel-optane-dc-persistent-memory-on-linux.html).  

### Acknowledgements
This pipeline makes use of several scientific libraries and tools, including:

- **`QUOKKA`**: For simulation data generation. See [QUOKKA GitHub Repository](https://github.com/quokka-astro/quokka).
- **`yt`**: For data ingestion and visualization. See [yt Project](https://yt-project.org/).
- **`pyxsim`**: For synthetic X-ray photon simulation. See [pyXSIM Documentation](http://www.ascl.net/1608.002).
- **`SOXS`**: For simulating X-ray instrument responses. See [SOXS Documentation](http://ascl.net/2301.024).
- **`numpy`**: For numerical computations. See [NumPy Official Website](https://numpy.org).
- **`matplotlib`**: For creating visualizations and plots. See [Matplotlib Official Website](https://matplotlib.org/).

Great thanks for [Dr. Aditi Vijayan](https://github.com/aditivijayan) and [Prof. Mark Krumholz](https://github.com/markkrumholz) for providing assistance in developing this pipeline. 

### Contact
For more information or to report issues, please contact [Rongjun.Huang@alumni.anu.edu.au](mailto:Rongjun.Huang@alumni.anu.edu.au) 
or [Astro@Rongjun-Huang.com](mailto:Astro@Rongjun-Huang.com).
```

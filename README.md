## Xpipeline: A [QUOKKA](https://github.com/quokka-astro/quokka) Simulation Post-Processing and X-ray Analysis Pipeline 

**Xpipeline** is a Python-based pipeline designed for the post-processing and analysis of QUOKKA simulation data in the soft X-ray band (0.4-2.0 keV). The repository includes a set of scripts that read, process, and simulate synthetic X-ray emissions, allowing researchers to explore various physical properties and simulate instrument responses.

### Key Features
- **Data Ingestion & Processing**: Reads simulation data and computes derived fields such as gas density, temperature, and emission measure.
- **Synthetic Photon Simulation**: Creates synthetic X-ray photons using various source models and projects them into events files.
- **Instrument Response Simulation**: Simulates X-ray instrument responses (e.g., Chandra ACIS) to generate realistic event, image, and spectrum files for visualization.

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

### Acknowledgements
This pipeline makes use of several scientific libraries and tools, including:
- **`yt`**: For data ingestion and visualization.
- **`pyxsim`**: For synthetic X-ray photon simulation.
- **`SOXS`**: For simulating X-ray instrument responses.

### Contact
For more information or to report issues, please contact [Rongjun.Huang@alumni.anu.edu.au](mailto:Rongjun.Huang@alumni.anu.edu.au) 
or [Astro@Rongjun-Huang.com](mailto:Astro@Rongjun-Huang.com).
```

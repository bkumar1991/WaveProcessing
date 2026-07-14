Supplementary_Code Package README
==================================

Title:
------
Automated Surface-Wave Waveform Extraction and Quality Control for Array Analysis

Author:
-------
Bhupender Kumar and Dr Aibing Li

Description:
------------
This supplementary package contains Python scripts used for:
1. Downloading and mapping seismic station metadata
2. Retrieving and filtering seismic events
3. Preparing datasets for waveform correlation and seismic analysis

The workflow is developed for preprocessing broadband seismic waveforms to generate high-quality surface wave ( both Rayleigh- and Love-wave) datasets for a regional seismic array. The method systematically isolates vertical (Z) and transverse (T) components, retaining only stations with complete tri-axial recordings and stable amplitude behaviors. Noisy or incomplete traces are automatically excluded through frequency-dependent filtering and correlation-based selection, ensuring coherent surface-wave arrivals across multiple frequency bands. The processed SAC files are organized by epicentral distance, enabling rapid visualization and spatial analysis of consistent, high signal-to-noise waveforms. Statistical evaluation of arrival-time attributes confirms improved signal coherence and reproducibility. This automated workflow delivers high-fidelity surface-wave observations, providing a robust foundation for phase-velocity inversion, tomographic modeling, and regional-scale imaging of subsurface structures.

Software Dependencies:
----------------------
Python 3.x with the following packages:

- obspy
- cartopy
- matplotlib
- pandas

Install dependencies using:

pip install obspy cartopy matplotlib pandas

------------------------------------------------------------
SCRIPT 1: Station Metadata Download and Mapping
------------------------------------------------------------

File:
-----
GetStations.py

Description:
------------
This script downloads broadband seismic station metadata from the IRIS FDSN web service for Alaska and surrounding regions during the year 2020. The script extracts station information for LH? channels and generates both station metadata files and a map of station distribution.

Processing Steps:
-----------------
1. Queries IRIS FDSN station service for LH? channel data
2. Time window: 2017-08-01 to 2018-12-31
3. Spatial coverage:
      Latitude : 56°N to 72°N
      Longitude: 170°W to 132°W
4. Extracts station metadata (network, station, latitude, longitude)
5. Writes station list to text file
6. Saves full station inventory in STATIONXML format
7. Generates a map of station distribution using Cartopy

Outputs:
--------
AEStations.txt
    Station list containing:
    network code, station code, longitude, latitude

AESta.xml
    Station metadata in STATIONXML format

Station distribution figure (Alaska region)

------------------------------------------------------------
SCRIPT 2: Event Selection and Filtering
------------------------------------------------------------

File:
-----
GetEvents.py

Description:
------------
This script retrieves seismic events from the IRIS FDSN event catalog and filters them based on magnitude, epicentral distance, and depth relative to a reference station.

Processing Steps:
-----------------
1. Queries IRIS event catalog for magnitude ≥ 6.0
2. Time window: 2017-08-01 to 2018-12-31
3. Defines reference station location:
      Latitude : 64.0°N
      Longitude: 152.0°W
4. Computes epicentral distance between station and each event
5. Applies selection criteria:
      - Distance: 30° to 120°
      - Depth: < 150 km
6. Builds filtered event catalog
7. Saves filtered catalog in QUAKEML and ZMAP formats
8. Generates global event distribution map

Outputs:
--------
AE2017-18.xml
    Filtered event catalog in QuakeML format

AE2017-18.txt
    Filtered event catalog in ZMAP format

Global event distribution figure

------------------------------------------------------------
SCRIPT 3: Event-Based Waveform Download and Data Assembly
------------------------------------------------------------

File:
------------------------------------------------------------
DataDownload_by_Event_Area.py

Description:
------------
This script downloads seismic waveform data for selected earthquake events and associated station metadata using the IRIS FDSN services. Events are filtered from a pre-processed catalog, and waveforms are
retrieved for each event within a specified time window and geographic region.

The workflow uses ObsPy’s MassDownloader module for efficient bulk waveform retrieval and organizes the output in an event-based directory structure for further processing and analysis. 

Data Sources:
------------
- IRIS FDSN Event Catalog (filtered ZMAP format)
- IRIS FDSN Waveform and Station Services

Mass Download Tool:
-------------------
This workflow uses ObsPy’s MassDownloader implementation for bulk data retrieval:

ObsPy MassDownloader:
https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.mass_downloader.html

(Implementation reference used in this study:http://gist.github.com/NeptuneProjects/MassDownloader.py)

Time Window:
-----------
2017-08-01 to 2018-12-31

Reference Station:
------------------
Latitude : 64.0°N
Longitude: 152.0°W

Spatial Coverage (Station Search):
----------------------------------
Latitude : 56°N to 72°N
Longitude: 170°W to 132°W

Processing Steps:
-----------------
1. Reads filtered earthquake catalog (AE2017-18.txt in ZMAP format)
2. Iterates through all events in the catalog
3. Computes epicentral distance between each event and reference station
4. Applies event selection based on time and distance
5. Assigns waveform window length depending on epicentral distance:
      - dist ≤ 60°   → 3000 s
      - 60–80°       → 3600 s
      - 80–100°      → 4200 s
      - > 100°       → 4800 s
6. Defines rectangular geographic domain for station selection
7. Uses ObsPy MassDownloader to retrieve waveform data (miniSEED)
8. Downloads corresponding StationXML metadata
9. Organizes outputs into event-based directories
10. Stores event metadata in CSV format for tracking and reproducibility

Output Structure:
-----------------
./AE2017-18/
    ├── catalog-*.csv
    ├── YYYY-MM-DD-HH-MM/
    │     ├── waveforms/
    │     │     └── miniSEED files
    │     └── stations/
    │           └── StationXML files

Generated Files:
---------------
catalog-*.csv
    Event metadata table containing:
    - origin time
    - latitude
    - longitude
    - depth (km)
    - magnitude
    - epicentral distance

Waveform Data:
--------------
MiniSEED files downloaded using ObsPy MassDownloader

Station Data:
------------
StationXML metadata corresponding to each event and station selection

Required Python Packages:
-------------------------
ObsPy 1.4.0
Matplotlib 
Pandas

Key ObsPy Modules Used:
-----------------------
- FDSN client (event queries)
- MassDownloader (bulk waveform retrieval)
- RectangularDomain (spatial filtering)
- Restrictions (data selection rules)

Usage:
------
python event_waveform_downloader.py

Notes:
------
- Internet connection is required for waveform and station downloads.
- MassDownloader enables efficient bulk retrieval of seismic data.
- Data are automatically structured by event origin time.
- This workflow ensures reproducible event-based waveform extraction.

SCRIPT 4: Waveform Processing and SAC File Generation
=================================================================

File:
-----
Data_Process_1.py

Description:
------------
The script extracts and organize the vertical and transverse components.  This algorithm reads the CSV file for event information and applies the same processing for each event.  It removes instrument responses, trim the data for the same window, and applies a zero phase bandpass filter of 0.002 to 0. 1 Hz. Waveform data are retained only for stations that record complete tri-axial components (E, N, and Z). Station information is extracted from the StationXML file, which is needed for determining event-station geometry and coordinate transformation.  Horizontal components are rotated from the NE to the RT frame. Then Z and T components are stored for further Rayleigh and Love wave processing.  

Data Inputs:
------------
- Event catalog (CSV format generated from previous workflow)
- MiniSEED waveform data (event-based directory structure)
- StationXML metadata files

Processing Steps:
-----------------
1. Reads event catalog containing:
      - origin time
      - latitude
      - longitude
      - depth
      - magnitude
      - epicentral distance

2. Iterates over each event directory:
      - loads waveform data (MiniSEED)
      - loads station metadata (StationXML)

3. Groups data by station and selects three-component records (ZNE)

4. Applies preprocessing:
      - instrument response removal (ObsPy remove_response)
      - spline detrending
      - bandpass filtering (0.002–0.1 Hz)
      - amplitude scaling

5. Performs data quality control:
      - removes traces with zero amplitude
      - checks amplitude ratios between components

6. Computes event–station geometry:
      - epicentral distance
      - backazimuth using IRIS distance/azimuth tools

7. Rotates horizontal components:
      - NE → RT coordinate system
      - generates Radial (R) and Transverse (T) components

8. Converts processed traces into SAC format and enriches SAC headers with:
      - event coordinates (evla, evlo, evdp, magnitude)
      - station coordinates (stla, stlo, stel)
      - epicentral distance (gcarc)
      - backazimuth (baz)

9. Separates waveforms into:
      - Rayleigh wave components (Radial/Z)
      - Love wave components (Transverse)

10. Sorts traces by epicentral distance and writes SAC files in ordered sequence for seismic section plotting

Output Structure:
-----------------
AE2017-18/
    ├── YYYY-MM-DD-HH-MM/
    │     ├── LHZ/
    │     │     └── D001.*.SAC
    │     ├── LHT/
    │     │     └── D001.*.SAC
    │     ├── waveforms/
    │     └── stations/

Generated Data:
---------------
LHZ directory:
    SAC files containing Rayleigh wave (vertical component)

LHT directory:
    SAC files containing Love wave (transverse component)

SAC Header Fields Added:
------------------------
- evla, evlo : event latitude and longitude
- evdp       : event depth (km)
- mag        : magnitude
- stla, stlo : station coordinates
- stel       : station elevation
- gcarc      : epicentral distance (degrees converted to km scale)
- baz        : backazimuth
- o          : origin time reference

Required Python Packages:
-------------------------
ObsPy (1.4.0)
NumPy
Matplotlib
Pandas

Key ObsPy Modules Used:
-----------------------
- Stream and Trace processing
- Instrument response removal
- Filtering and detrending
- Rotation (NE → RT)
- SACTrace conversion

Usage:
------
python waveform_processing_and_rotation.py

Notes:
------
- Only three-component stations are used for analysis.
- Data quality control is performed using amplitude ratio checks.
- Waveforms are organized by increasing epicentral distance.
- Output SAC files are structured for seismic section plotting.
- This script is intended for Rayleigh and Love wave separation studies.

SCRIPT 5: Waveform Extraction and Retention
===============================================================================

File:
-----
Data_Process_2.py

Description:
------------
This script performs frequency-dependent surface-wave extraction and quality control for processed SAC waveform data output from Wavefoem_Processing_1.py The workflow separates Rayleigh and Love wave energy, applies narrow-band filtering, identifies signal windows using waveform envelopes, removes noisy traces, and generates filtered datasets waveform analysis.

Input Data:
-----------
- Event catalog in CSV format
- Processed SAC waveform files from Waveform_Processing_1.py
- Rayleigh-wave directories (LHZ)
- Love-wave directories (LHT)

Processing Workflow:
--------------------
1. Reads event catalog and event metadata
2. Loads processed SAC files for each event
3. Separates Rayleigh (Z/R) and Love (T) wave datasets
4. Defines target center frequencies for narrow-band analysis
5. Applies bandpass filtering around each selected frequency
6. Estimates theoretical phase-arrival times using reference phase velocities
7. Computes waveform envelopes for each trace
8. Identifies maximum envelope amplitudes and corresponding arrival times
9. Determines adaptive signal windows using local envelope minima
10. Applies time-window tapering to isolate surface-wave energy
11. Removes pre-signal and post-signal noise
12. Performs signal-to-noise quality control using:
      - Amax/Aavg amplitude ratios
      - phase-arrival consistency
      - timing residual thresholds
13. Removes outlier traces based on statistical deviation
14. Saves filtered SAC waveforms for accepted traces
15. Generates section plots for:
      - all filtered traces
      - accepted high-quality traces

Frequency Analysis:
-------------------
Example frequency used in this workflow:

0.0225 Hz

Bandpass Window:
----------------
freqmin = center_frequency - 0.005 Hz
freqmax = center_frequency + 0.005 Hz

Reference Surface-Wave Velocities:
----------------------------------
Rayleigh wave velocity:
    3.95 km/s

Love wave velocity:
    4.27 km/s

Quality-Control Criteria:
-------------------------
Rayleigh and Love wave traces are retained only if:
- envelope amplitude ratio (Amax/Aavg) exceeds threshold
- observed arrival time is close to predicted phase arrival
- statistical outliers are removed using:
      |dtt - mean(dtt)| < 2.5 × std(dtt)

Signal Windowing:
-----------------
The signal window is determined automatically using:
- envelope maxima
- local minima surrounding the maximum amplitude
- adaptive tapering

Output Structure:
-----------------
EventDirectory/
    ├── LHZ/
    │     ├── f-0.0225/
    │     │     ├── filtered SAC files
    │     │     ├── section plots
    │     │     └── retained trace plots
    │
    ├── LHT/
    │     ├── f-0.0225/
    │     │     ├── filtered SAC files
    │     │     ├── section plots
    │     │     └── retained trace plots

Generated Files:
---------------
Filtered SAC files:
    Narrow-band filtered waveforms for accepted traces

SVG Figures:
------------
f-<frequency>.svg
    Section plot for all filtered traces

f-<frequency>kp.svg
    Section plot for accepted high-quality traces

Output Log:
-----------
catalog-1.output
    Processing summary including:
    - event information
    - wave type
    - frequency
    - number of retained traces

Key Signal-Processing Methods:
------------------------------
- Narrow-band bandpass filtering
- Envelope analysis
- Adaptive signal-window selection
- Time-domain tapering
- Surface-wave phase prediction
- Statistical outlier rejection

Required Python Packages:
-------------------------
ObsPy 1.4.0
NumPy
Matplotlib
Pandas

Key ObsPy Modules Used:
-----------------------
- Stream and Trace processing
- SAC waveform handling
- Envelope computation
- Trigger and signal utilities
- Filtering and plotting

Usage:
------
python Waveform_Processing_2.py

Notes:
------
- Filtered traces are grouped by frequency and wave type.
- Directories containing fewer than 20 accepted traces are automatically removed.
- Output SAC files retain original metadata and distance ordering.
- The script performs automated quality control for large seismic datasets.

SCRIPT 6: Waveform Quality control and Final Trace selection
===============================================================================

File:
-----
Data_Process_Final.py

Description:
------------
The script evaluates the traces of vertical or transverse components for each event at a given frequency as a group and identifies the extremely anomalous traces based on the coherence of arrival times and phase shifts. It iterates through the event catalog, the components, and the frequencies. For each identified folder, the associated SAC files are copied into a mirrored structure to ensure that both raw and processed data remain structurally consistent. 

Input Data:
-----------
- Processed SAC waveform files from Waveform_Processing_2.py
- Rayleigh-wave directories (LHZ) having frequency dependent subfolders (example f0.0225) containing the Processed SAC waveform files from Waveform_Processing_2.py
- Love-wave directories (LHT) having frequency dependent subfolders (example f0.0225) containing the Processed SAC waveform files from Waveform_Processing_2.py

Processing Workflow:
--------------------
The script:

1. Reads filtered SAC waveform files from frequency-dependent folders (example f-0.015) generated as output in Data_Process_2.py
2. Computes analytic signal envelopes of each trace recorded by different seismic stations for each individual events using the Hilbert transform 
3. Determines envelope peak amplitudes and peak arrival times
4. The script performs weighted least-squares matrix inversion to estimate the expected trend between trace offset distance and envelope peak arrival time. A linear model is fitted to the    data to identify anomalous traces having high standard deviation values
5. Residuals between observed and predicted peak times are computed. calculates residual mean and standard deviation, estimates a robust standard deviation and automatically determines outlier thresholds. Identifies traces with abnormal envelope arrival behavior. Detected outlier traces are removed automatically. 
6. Performs waveform polarity/coherence analysis using cross-correlation. The algorithm compares neighboring traces, tracks waveform coherence trends, identifies polarity-reversed traces, detects incoherent waveform segments and allows new coherent waveform segments when persistent polarity changes occur. Traces classified as polarity reversed or incoherent are removed automatically.
7. Produces cleaned waveform datasets and diagnostic figures

Output Structure:
-----------------
AE2017-18Processed/
│
├── EventDirectory/
│   ├── LHZ/
│   │   └── f-0.0225/
│   │       ├── D001.AK.STA.SAC
│   │       ├── D002.AK.STA.SAC
│   │       └── ...
│   │
│   └── LHT/
│       └── f-0.0225/
│           ├── D001.AK.STA.SAC
│           └── ...
│
└──EventDirectory/

Generated Files:
---------------
For each processed frequency folder, the script generates:

File	                                  Description
Lags.txt / lags.txt	    Correlation and polarity analysis results
final_traces.tif	            Final cleaned section plot
trace_with_outlier.tif	    Section plot before outlier removal
Residual_Trace.tif	    Section plot after residual filtering
Correlation Vs Offset.tif	    Correlation coefficient versus offset
peak_time.tif	            Observed vs predicted envelope peak times

Required Python Packages:
-------------------------
ObsPy (1.4.0)
NumPy
Matplotlib
Pandas

Key ObsPy Modules Used:
-----------------------

- waveform I/O (`read`, `Stream`, `Trace`), 
-SAC file handling (`SACTrace`), signal-processing utilities, 
-event/station metadata tools, and IRIS/FDSN client services for distance calculations and seismic data management.

Usage:
------

python Waveform_Processing_Final.py

The script will Search all event folders, Process LHZ/f-0.0225 and LHT/f-0.0225, Remove residual outlier and polarity-reversed traces and finally, Saves cleaned datasets and figures into the output directory

Notes:
------
-Thresholds for coherence and outlier rejection may be adjusted depending on network geometry and signal quality.
-Diagnostic figures are saved at high-quality resolution (dpi=600).







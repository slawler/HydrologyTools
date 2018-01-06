# HydrologyTools
#### Notebooks for working with USGS stream gage records & USACE modelling software. 
*Source code written in Python & R*

This repository was created to document and share the development of open source tools used for the 
analysis and modeling of riverine flood hazards. 

Use cases included:
1. Downloading historical USGS gage data (instantaneous, daily & peak) for multiple gages
using the R package [dataRetrieval](https://github.com/USGS-R/dataRetrieval) developed by the USGS.

2. Accessing, analyzing and plotting data directly from [HEC-RAS](http://www.hec.usace.army.mil/software/hec-ras/) model simulations.
3. __*More to come...*__

## Notebooks
##### Data Tools
+  [Download USGS Gage Data](nbs/GageGrabber.ipynb) (R)
+  [Explore USGS Gage Data](nbs/GageGrabber.ipynb) (Python)

##### Discharge Hydrograph Development
+  [Use historic gage data to develop a return period hydrograph](nbs/MethodologyOverview.ipynb)

##### Levee Breach Hydrograph Development
+  [Example # 1. Hypothetical Breach Hydrograph using 1D Steady State model results](nbs/Lisle_BreachHydro.ipynb)
+  [Example # 2. Hypothetical Breach Hydrograph using 1D Steady State model results](nbs/WP_BreachHydro.ipynb)
+  [Example # 3. Hypothetical Breach Hydrograph using data from HEC-RAS Unsteady Simulation](nbs/Deposit_BreachHydro.ipynb)

## Requirements
To use the Notebooks available in this repository, the following software is required:
- Python
- R
- Jupyter

All of these open source tools are available at no cost via the [Anaconda Distribution](https://www.anaconda.com/distribution/).


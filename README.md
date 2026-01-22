# FWI: Full Waveform Inversion for Marine Data

## Features

- The software offers the possibility to invert a p-wave velocity model from marine field data.
- HPC is implemented in the code.

## Prerequisites

Before installing, ensure you have the following dependencies:
- Fortran Compilers (gfortran/mpif90/mpirun).
- Seismic Unix (SU): The software handles data in `.su` format.
- Make: For build automation.

## Installation
- Navigate to the source directory and compile:
- cd FWI/src
- make

## Inputs
- The software requires seismic marine field data (WAS,MCS) in Seismic Unix (.su) format.
- Navigation: An ASCII file containing the geometry/navigation information.
- Parameter File: A configuration file (read during execution) where the user specify the input parameters, velocity model, grid dimensions, etc. An input file "parfile_fwi" is included in this folder as an example.

## Run
- Example (run with 40 MPI processes):
- mpirun -np 40 psdm_run parfile_fwi
 
## Documentation and testing
- Manual: A detailed PDF manual is currently under construction and will be added to the repository soon.
- Test Data: Sample datasets for testing the code will be hosted on the Zenodo database (link forthcoming).

## Author & Acknowledge
- Author: Clara Estela Jiménez Tejero.
- Institution: Barcelona Center for Subsurface Imaging, ICM-CSIC.
- While this software is freely available, we would be grateful if you could notify us of its use. Please send a brief email to ejimenez@icm.csic.es to acknowledge your usage, which helps us justify continued development and support.

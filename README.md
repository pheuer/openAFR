# openAFR
 Open-source analysis code for angular filter refractometry analysis
 
This code is very much a work in progress - please contact me with suggestions bugs.

# Installation (development mode)

The best way to install the package is in development mode. Clone the
repository onto your computer (ideally using GitHub Desktop) and then follow the 
instructions below. In this development mode, any changes
made to the local files (local edits or updates by pulling from github) will be
automatically reflected in the package.

## Using Anaconda
1. Open an Anaconda prompt (Windows) or a regular terminal prompt (Mac, Linux)
2. Install the package dependencies by running this command from the anaconda terminal: 'conda install matplotlib numpy scipy astropy h5py'
4. Run the command 'conda develop _' where _ is the package directory path. eg.
'conda develop /Users/Username/Documents/GitHub/openAFR'
5. Open Anaconda (close and restart if necessary) and test that the package has been
added to the path by importing something, eg. "from openAFR.afrfilter import AFRFilter
6. If necessary, packages can be uninstalled with 'conda develop -u _'

## Using pip
(Instructions untested)
(Anaconda is recommened - only use these instructions if you don't have Anaconda). 
1. Open a command prompt
2. Install package dependencies 'pip install matplotlib numpy scipy astropy h5py'
3. Run the command 'pip install -e _' where _ is the path to the package directory, eg. 
'conda develop /Users/Username/Documents/GitHub/openAFR'
 
# Downloading and installing h5toh4convert

The openAFR routines require HDF5 files, while the current output format for datafiles on the OMEGA EP facility is HDF4. 

openAFR includes a function `ensure_hdf5` that use the HDFGroup's h4toh5convert function to convert hdf4 files () to hdf5. 

1) Download and install the lastest version of the conversion utility from
   the HDF group: https://portal.hdfgroup.org/display/support/h4h5tools%202.2.5#files
   
2) Add the \bin file to your path. This should be something like 
   C:\Program Files\HDF_Group\H4TOH5\2.2.2\bin

3) Test by running "h4toh5convert -h" in the command line

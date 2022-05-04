# openAFR
 Open-source analysis code for angular filter refractometry analysis
 
 * This is very much a work in progress * ! 
 
# Downloading and installing h5toh4convert

The openAFR routines require HDF5 files, while the current output format for datafiles on the OMEGA EP facility is HDF4. 

openAFR includes a function `ensure_hdf5` that use the HDFGroup's h4toh5convert function to convert hdf4 files () to hdf5. 

1) Download and install the lastest version of the conversion utility from
   the HDF group: https://portal.hdfgroup.org/display/support/h4h5tools%202.2.5#files
   
2) Add the \bin file to your path. This should be something like 
   C:\Program Files\HDF_Group\H4TOH5\2.2.2\bin

3) Test by running "h4toh5convert -h" in the command line

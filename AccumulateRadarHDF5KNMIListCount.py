#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python script
# Name: AccumulateRadarHDF5KNMIListCount.py
# Part of MapRAINLINK: https://github.com/overeem11/MapRAINLINK
#
#
## Version 1.1
## Copyright (C) 2025 Aart Overeem
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#
# Description: Script to accumulate KNMI HDF5 radar rainfall images where the file names are provided as input or to accumulate KNMI HDF5 radar rainfall images for a given path. 
#              Read KNMI HDF5 radar file, extract rainfall depth, accumulate rainfall depth over provided input files.
#              Output is a KNMI HDF5 radar file with the accumulated rainfall depth.
#              Use a conversion factor of 1 when input files contain accumulated rainfall (mm).
#              Use a conversion factor of 0.25 when input files contain rainfall intensity in mm per hour.
#              nan and nodata are set to 0, so all data are accumulated irrespective of their value (e.g. missing data), but
#              the needed minimum availability per radar pixel is employed in the end to decide whether a radar pixel should have data or nodata. 
#              Also in case of missing files an output file is made. If data criterion is not satisfied nodata values are present in the output file. 
#              Note that we specifically deal with NaN values by setting them to nodata.
#              The maximum number of radar groups occurring in the list of files is stored, together with the radar IDs (if available) and radar locations.
#              calibration_missing_data is often 65534 or 65535. When template file with 65535 is used, it will become 65535 in case of nodata, i.e., 0 images with data values.
#              No thresholding is applied to maximize intensity at e.g. 100 mm/h. No isolated pixels are removed. No 5-point median clutter filter is applied.
#              Lower threshold is applied though, e.g. to select >= 0.1 mm/h, i.e., 0.1/12 mm per 5 minutes.
# Usage: python AccumulateRadarHDF5KNMIListCount.py [output filename] [input file names or path with files] [conversion factor] 
# [needed minimum number of images with data] [input file which is used to construct output file in case all input files are not valid] [path or file names?: choose "path" or "files"] [product_group_name]
# [threshold value in mm or mm/h, i.e., 0.1, but also take into account scale & offset, i.e. apply these to the threshold which is supplied!]
# Example: python AccumulateRadarHDF5KNMIListCount.py "RAD_NL25_RAC_MFBS_15min_201109102045_NL.h5" "RAD_NL25_RAC_MFBS_5min_201109102035_NL.h5  RAD_NL25_RAC_MFBS_5min_201109102045_NL.h5 RAD_NL25_RAC_MFBS_5min_201109102040_NL.h5" 1 1 RAD_NL25_RAC_RT_202107141500.h5 files RAD_NL25_RAC_MFBS_15min 0


# Load Python packages:
import sys
import os
import numpy as np
import shutil
from pathlib import Path
import warnings
#warnings.filterwarnings("ignore")
import h5py
import natsort


# Parameters from command line:    
OutputFileName = sys.argv[1]
InputFileNames = sys.argv[2]
ConversionFactor = float(sys.argv[3])
MinImages = int(sys.argv[4])
InputFileNameNodata = sys.argv[5]
PathOrFileNames = sys.argv[6]
product_group_name = sys.argv[7]
threshold = float(sys.argv[8])


if PathOrFileNames=="files":
    # Split string with input file names in string with multiple lines:
    pathlist_temp = InputFileNames.split()
    # Sort file names:
    pathlist = natsort.natsorted(pathlist_temp)

if PathOrFileNames=="path":
    pathlist_temp = Path(InputFileNames).glob('**/*.h5')
    # Sort file names:
    pathlist = natsort.natsorted(pathlist_temp)

if PathOrFileNames!="files" and PathOrFileNames!="path":
    print("Accumulation cannot be performed. Please specify whether a list of files (files) or a directory path with files (path) is supplied!")
    sys.exit(0)



#############################################################################
# Read HDF5 radar files (KNMI HDF5 format) and accumulate to rainfall depth.#
#############################################################################


i = 0
DATAFIELD_NAME = '/image1/image_data'
radardata = radardata_temp = Count = []
for path in pathlist:
    if os.path.getsize(path) > 0:               # Is used to check the size of specified path. It returns the size of specified path in bytes.
    						# The method raises OSError if the file does not exist or is somehow inaccessible.
       if h5py.is_hdf5(path):                   # Check that a file is a valid HDF5 file.
          i = i + 1
          print(path)
          if i==1:
             # Open file:
             f = h5py.File(path, "r")
             # Read data:
             nodata = f['/image1/calibration'].attrs['calibration_out_of_image']
             missingdata = f['/image1/calibration'].attrs['calibration_missing_data']             
             radardata = np.int_(f[DATAFIELD_NAME][:])
             # Set missing data to 0:             
             truth_table = radardata==missingdata
             indices = np.where(truth_table)
             radardata[indices] = 0             
             # Replace nan with nodata, and set nodata to 0:
             radardata[np.isnan(radardata)] = nodata
             truth_table = radardata==nodata
             indices = np.where(truth_table)
             radardata[indices] = 0
             # Set radardata values below threshold to 0:
             truth_table_freq = radardata < threshold
             indices_freq = np.where(truth_table_freq)
             radardata[indices_freq] = 0              
             startdatetime = f['/overview'].attrs['product_datetime_start']
             enddatetime = f['/overview'].attrs['product_datetime_end']             
             # Count for each radar pixel the number of valid values, i.e. not equal to nodata or nan:
             Count = np.int_(f[DATAFIELD_NAME][:])
             Count[np.isnan(Count)] = nodata
             truth_table = Count!=nodata
             indices = np.where(truth_table)
             Count[indices] = 1
             truth_table = Count==nodata
             indices = np.where(truth_table)
             Count[indices] = 0
             # Determine number of radars:
             number_radar_groups = f['/overview'].attrs['number_radar_groups']
             InputFileName = str(path)
             f.close()
          else:
             # Open file:
             f = h5py.File(path, mode='r')
             # Read data:  
             nodata = f['/image1/calibration'].attrs['calibration_out_of_image']
             missingdata = f['/image1/calibration'].attrs['calibration_missing_data']              
             radardata_temp = []  
             radardata_temp = np.int_(f[DATAFIELD_NAME][:])
             # Set missing data to 0:             
             truth_table = radardata_temp==missingdata
             indices = np.where(truth_table)
             radardata_temp[indices] = 0               
             # Replace nan with nodata, and set nodata to 0:
             radardata_temp[np.isnan(radardata_temp)] = nodata
             truth_table = radardata_temp==nodata
             indices = np.where(truth_table)
             radardata_temp[indices] = 0
             # Set radardata values below threshold to 0:
             truth_table_freq = radardata_temp < threshold
             indices_freq = np.where(truth_table_freq)
             radardata_temp[indices_freq] = 0               
             radardata = np.add(radardata,radardata_temp)
             enddatetime = f['/overview'].attrs['product_datetime_end']
             # Count for each radar pixel the number of valid values, i.e. not equal to nodata or nan:
             Count_temp = np.int_(f[DATAFIELD_NAME][:])
             Count_temp[np.isnan(Count_temp)] = nodata
             truth_table = Count_temp!=nodata
             indices = np.where(truth_table)
             Count[indices] = Count[indices] + 1
             # Determine number of radars and if larger than in previous file(s), store the maximum value:
             number_radar_groups_temp = f['/overview'].attrs['number_radar_groups']             
             if number_radar_groups_temp > number_radar_groups:
                number_radar_groups = number_radar_groups_temp
                InputFileName = str(path)
             f.close()


print(i)

##############################################
# Make output radar file in ODIM HDF5 format.#
##############################################
if i > 0:        # Implies that at least 1 radar image can be read, but can contain nodata only.
   # Apply conversion factor:
   radardata = radardata * ConversionFactor
   # Apply data availability criterion separately per pixel. If availability too low, set radardata to nodata:
   truth_table = Count<MinImages
   indices = np.where(truth_table)    
   radardata[indices] = nodata
   # Remove output file:
   try:
       os.remove(OutputFileName)
   except OSError:
       pass
   # Copy input file to HDF5 output file:       
   shutil.copy(InputFileName, OutputFileName)  
   # Note that metadata is copied from a file in the accumulation period.   
   # Make output file:
   hf = h5py.File(OutputFileName, "a")
   # Remove data field:	
   del hf[DATAFIELD_NAME]
   # Create data field including attributes and write to HDF5 output file:
   dset = hf.create_dataset(DATAFIELD_NAME, data=radardata, compression="gzip", compression_opts=6)
   dset.attrs["CLASS"] = np.string_("IMAGE ")
   dset.attrs["IMAGE_VERSION"] = np.string_("1.2 ")
   dset = hf["image1"]
   dset.attrs.modify('image_product_name',product_group_name)
   # Modify attributes:
   dset = hf["overview"]
   # If all files are available the start and end date & time will be correct.
   # Note that for the start date and time the first file which could be read is used, whereas for the end date and time the last file which could be read is used.
   # So in case of missing files, the start and/or end date & time can become different. In case of a missing file between start and end date & time, this cannot be seen
   # in the attributes enddate, endtime, startdate, and starttime.
   dset.attrs.modify('product_datetime_start',startdatetime)
   dset.attrs.modify('product_datetime_end',enddatetime)
   dset.attrs.modify('product_group_name',product_group_name) 
   #
   hf.close()


if i==0:                 # If none of the input files can be read, make an output file with metadata and the radardata field being entirely nodata:
   print("All radar files are empty. Hence, an accumulated radar rainfall image could not be produced, but a file with nodata entries is produced.")
   if os.path.getsize(InputFileNameNodata) > 0:               # Is used to check the size of specified path. It returns the size of specified path in bytes.
                                                              # The method raises OSError if the file does not exist or is somehow inaccessible.
      if h5py.is_hdf5(InputFileNameNodata):                   # Check that a file is a valid HDF5 file.
         # Open file:
         f = h5py.File(InputFileNameNodata, "r")
         # Read data:
         nodata = f['/image1/calibration'].attrs['calibration_out_of_image']
         radardata = np.int_(f[DATAFIELD_NAME][:])
         truth_table = radardata!=nodata
         indices = np.where(truth_table)
         radardata[indices] = nodata
         f.close()
         # Remove output file:
         try:
             os.remove(OutputFileName)
         except OSError:
             pass
         # Copy input file to HDF5 output file:       
         shutil.copy(InputFileNameNodata, OutputFileName) 
         # Note that metadata is copied from a file from another period, where some attributes may not be right, but this is not really relevant given the nodata field.   
         # Make output file:
         hf = h5py.File(OutputFileName, "a")
         # Remove data field:	
         del hf[DATAFIELD_NAME]
         # Create data field including attributes and write to HDF5 output file:
         dset = hf.create_dataset(DATAFIELD_NAME, data=radardata, compression="gzip", compression_opts=6)
         dset.attrs["CLASS"] = np.string_("IMAGE ")
         dset.attrs["IMAGE_VERSION"] = np.string_("1.2 ") 
         dset = hf["image1"]
         dset.attrs.modify('image_product_name',product_group_name)
         # Modify attributes:
         dset = hf["overview"]                  
         # Since the date & time metadata are generally not used, the date & time from the original file are kept, which are usually from a different period.
         # Note that this only occurs if all files for a given accumulation interval cannot be opened.
         dset.attrs.modify('product_group_name',product_group_name)
         #
         hf.close()
      else:
         print("Nodata input file cannot be opened! No output file is constructed!")
   else:
      print("Nodata input file cannot be opened! No output file is constructed!")







#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python script
# Name: VisualizeCMLsGaugesRadars_CartopyOSM_GM.py
# Part of MapRAINLINK: https://github.com/overeem11/MapRAINLINK
#
#
## Version 1.11
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
# Description: Python script to visualize rain gauge, radar, and commercial microwave link locations and their rainfall estimates on a map.
#	       Script to visualize commercial microwave link (CML) and/or rain gauge locations on a map (or only plot a map).
#              Can also visualize rainfall for rain gauge locations, CML paths, and/or CML maps.
#              Instead of CML maps, radar precipitation maps can be visualized from either data in KNMI-HDF5 format or data in ODIM-HDF5 format.
#              This implies that Dutch gridded radar precipitation data or OPERA gridded radar precipitation data can be visualized.
#              Rain gauge and CML data are read from text files.
#              All settings should be supplied by using a Python configuration script as argument to this program.
#              This Python script is more widely applicable. It can be used to plot lines and/or points (and their values) on a map, not necessarily commercial microwave link paths 
#              and rain gauge locations, and not necessarily rainfall intensities or accumulations. It can also be used to plot other data fields.
#              Plotting 1 map takes typically 5-10 seconds, or up to 15 seconds when plotting over entire Europe. 
# Usage: python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py [additional settings can be added, which override the configuration script]
# Example: python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py
# Example, where additional settings are supplied, which override those in "ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py" is provided below. 
# This can be convenient when maps for different intervals need to be plotted. Then there is no need to make a configuration file for each interval, 
# but the settings that change only need to be changed in the command line, so typically input and output filenames, and extra text.
# To apply this to many files, a user could generate command lines, where the same configuration file is used, but only a couple of settings are changed in the command line.
# Example: python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py "ColorLinkPaths = 'darkgreen'" "ScaleType = 'YellowRed'" "ColorGaugeNetwork3 = 'brown'" "OutputFileName = 'test.jpg'"



# Load Python packages:
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import cartopy.mpl.ticker as cticker
from matplotlib.colors import BoundaryNorm
from matplotlib.collections import LineCollection
import numpy as np
from numpy import loadtxt
import copy
import h5py
import pyproj
import pandas as pd
# Added this static backend (https://matplotlib.org/stable/users/explain/figure/backends.html#static-backends) to avoid that command line is possibly not cleared after running this script in command line in remote sessions:
mpl.use('Cairo')


# Loading Python configuration script with parameter settings:
exec(open(sys.argv[1]).read())


# Additional arguments to script, which will override settings in configuration script
# (useful when visualizations are performed for different time intervals, so with different input & output files):
if len(sys.argv)>2:
   with open("tempInput.py", "w") as f:
      f.write("#!/usr/bin/env python3\n")
      f.write("# -*- coding: utf-8 -*-\n")
      for i in range(2,len(sys.argv)):
         f.write(sys.argv[i])
         f.write("\n")
   f.close()
   # Loading Python configuration script parameter settings of separate arguments used with this script and remove it:
   exec(open("tempInput.py").read())
   os.remove("tempInput.py")


#############################################
# 1. Read locations of CMLs and rain gauges.#
#############################################


if PlotCML=="yes":
   # Read locations of CMLs:
   startlon, startlat, endlon, endlat = loadtxt(FileNameCMLLocations, comments="#", delimiter=" ", usecols=(0, 1, 2, 3), unpack=True, skiprows=1)

if PlotGaugeNetwork1=="yes":
   station_nameNetwork1Temp = loadtxt(FileNameGaugeLocationsNetwork1, comments="#", delimiter=",", usecols=(2), dtype=str, skiprows=1)
   station_nameNetwork1 = [x[:NumberOfCharactersNamesNetwork1] for x in station_nameNetwork1Temp]
   if PlotValuesGaugesNetwork1=="yes":
      # Import coordinates of weater stations (longitude and latitude in degrees) and rainfall for network 1:
      lon_gaugesNetwork1, lat_gaugesNetwork1, GaugeValues1 = loadtxt(FileNameGaugeLocationsNetwork1, comments="#", delimiter=",", usecols=(0, 1, 3), unpack=True, skiprows=1)
      GaugeValues1 = GaugeValues1 * ConversionFactorGaugesNetwork1
      if DoColorSetUnder!="yes":
         GaugeValues1[GaugeValues1 < LowestValue] = np.nan      
   else:
      # Import coordinates of weater stations (longitude and latitude in degrees) for network 1:
      lon_gaugesNetwork1, lat_gaugesNetwork1 = loadtxt(FileNameGaugeLocationsNetwork1, comments="#", delimiter=",", usecols=(0, 1), unpack=True, skiprows=1)

if PlotGaugeNetwork2=="yes":
   station_nameNetwork2Temp = loadtxt(FileNameGaugeLocationsNetwork2, comments="#", delimiter=",", usecols=(2), dtype=str, skiprows=1)
   station_nameNetwork2 = [x[:NumberOfCharactersNamesNetwork2] for x in station_nameNetwork2Temp]
   if PlotValuesGaugesNetwork2=="yes":
      # Import coordinates of weater stations (longitude and latitude in degrees) and rainfall values for network 2:
      lon_gaugesNetwork2, lat_gaugesNetwork2, GaugeValues2 = loadtxt(FileNameGaugeLocationsNetwork2, comments="#", delimiter=",", usecols=(0, 1, 3), unpack=True, skiprows=1)
      GaugeValues2 = GaugeValues2 * ConversionFactorGaugesNetwork2
      if DoColorSetUnder!="yes":
         GaugeValues2[GaugeValues2 < LowestValue] = np.nan       
   else:
      # Import coordinates of rain gauges (longitude and latitude in degrees) for network 2:
      lon_gaugesNetwork2, lat_gaugesNetwork2 = loadtxt(FileNameGaugeLocationsNetwork2, comments="#", delimiter=",", usecols=(0, 1), unpack=True, skiprows=1)

if PlotGaugeNetwork3=="yes":
   station_nameNetwork3Temp = loadtxt(FileNameGaugeLocationsNetwork3, comments="#", delimiter=",", usecols=(2), dtype=str, skiprows=1)
   station_nameNetwork3 = [x[:NumberOfCharactersNamesNetwork3] for x in station_nameNetwork3Temp]
   # Import coordinates of rain gauges (longitude and latitude in degrees) for network 3:
   lon_gaugesNetwork3, lat_gaugesNetwork3 = loadtxt(FileNameGaugeLocationsNetwork3, comments="#", delimiter=",", usecols=(0, 1), unpack=True, skiprows=1)
          


########################################
# 2. Read CML data and obtain polygons.#
########################################

if PlotDataField=="yes" and PlotDataFieldRadarGrid!="yes":
   # Import longitude & latitude and other features of interpolation grid:
   Lon, Lat, IndexLon, IndexLat, GridCML = loadtxt(InterpolationGrid, comments="#", delimiter=",", usecols=(0, 1, 2, 3, 4), unpack=True)
   # Import interpolated field of path-averaged rainfall intensities from CMLs (mm/h):
   data = loadtxt(InputFileNameDataField, comments="#", delimiter=",", usecols=(0), unpack=True, skiprows=1)

   Lon = Lon - GridRes/2
   Lat = Lat - GridRes/2
   # The grid on which each field of values is presented is, e.g., a 0.01x0.01 degrees lat./lon. (Cylindrical Equal Distance) array of points.
   IndexLon = IndexLon.astype(int)
   IndexLat = IndexLat.astype(int)

   LonArray = [[np.nan for x in range(max(IndexLon)+1)] for x in range(max(IndexLat)+1)]
   LatArray = [[np.nan for x in range(max(IndexLon)+1)] for x in range(max(IndexLat)+1)]
   RArray = [[np.nan for x in range(max(IndexLon)+1)] for x in range(max(IndexLat)+1)]
   for k in range(0,Lon.shape[0]):
      # Coordinate of pixel represents corner of pixel.
      LonArray[IndexLat[k]][IndexLon[k]] = Lon[k]
      LatArray[IndexLat[k]][IndexLon[k]] = Lat[k]
      if GridCML[k]==1:
         RArray[IndexLat[k]][IndexLon[k]] = data[k]
      if GridCML[k]==0:
         RArray[IndexLat[k]][IndexLon[k]] = 0
         # The dimensions of X and Y should be one greater than those of C. Alternatively, X, Y and C may have equal dimensions, in which case the last row and column of C will be ignored.
         # This is the case here.

   # No data values are made "not a number".
   RArray = np.array(RArray) * ConversionFactorCMLDataField
   if DoColorSetUnder!="yes":
      RArray[RArray < LowestValue] = np.nan

      
if PlotDataField=="yes" and PlotDataFieldRadarGrid=="yes":
   # Import interpolated field of path-averaged rainfall intensities from CMLs (mm/h):
   RArrayTemp = loadtxt(InputFileNameDataField, comments="#", delimiter=",", usecols=(0), unpack=True, skiprows=1)   
   RArrayTemp = np.array(RArrayTemp) * ConversionFactorCMLDataField


if PlotCMLTimeInterval=="yes":
   # Read locations of CMLs (if Run.R from RAINLINK is followed, NAs do not exist (see na.omit in Run.R); RainfallDepthPath = NA may be a problem here):
   startlon, startlat, endlon, endlat, RainfallDepthPath = loadtxt(FileNameCMLLocationsTimeInterval, comments="#", delimiter=" ", usecols=(3, 4, 5, 6, 1), unpack=True, skiprows=1)
   RainfallDepthPath = RainfallDepthPath * ConversionFactorCMLPathAverage



#############################################################
# 3. Read HDF5 radar file (KNMI format) and obtain polygons.#
#############################################################

if PlotKNMIRadar=="yes" or PlotDataFieldRadarGrid=="yes":
   FILE_NAME = KNMIRadarInputFileName
   f = h5py.File(FILE_NAME, mode='r')
   # Read metadata:    
   xscale = f['/geographic'].attrs['geo_pixel_size_x'][0]
   yscale = f['/geographic'].attrs['geo_pixel_size_y'][0]
   xoffset = f['/geographic'].attrs['geo_column_offset'][0]
   yoffset = f['/geographic'].attrs['geo_row_offset'][0]
   proj4str = f['/geographic/map_projection'].attrs['projection_proj4_params']
   Ncols = f['/geographic'].attrs['geo_number_columns']
   Nrows = f['/geographic'].attrs['geo_number_rows']

   my_string = str(f[DATAFIELD_NAMECAL].attrs['calibration_formulas'])
   string_temp = [str(x) for x in my_string.split('=')][1]
   zscale = float([str(x) for x in string_temp.split('*')][0])
   string_temp = [str(x) for x in my_string.split('+')][1]
   try:
      zoffset = float(string_temp.replace("'",""))
   except:
      pass
   try:
      zoffset = float(string_temp.replace("']",""))
   except:
      pass
   nodata = zoffset + zscale * f[DATAFIELD_NAMECAL].attrs['calibration_out_of_image']
   missingdata = zoffset + zscale * f[DATAFIELD_NAMECAL].attrs['calibration_missing_data']

   # Read data:
   dset = f[DATAFIELD_NAME]
   RArray = zoffset + zscale * dset[:]


   # Convert image coordinates to longitude & latitude in degrees.
   p = pyproj.Proj(proj4str.decode('utf-8'))
   LonArray = [[np.nan for x in range(Ncols[0]+1)] for x in range(Nrows[0]+1)]
   LatArray = [[np.nan for x in range(Ncols[0]+1)] for x in range(Nrows[0]+1)]
   RArrayCML = [[np.nan for x in range(Ncols[0]+1)] for x in range(Nrows[0]+1)]   
   # Obtain image coordinates of surrounding radar pixels needed to compute the coordinates of the radar polygons.
   Nrow = Nrows[0]+1
   Ncol = Ncols[0]+1
   for j in range(0,Nrow):
      for i in range(0,Ncol):
         # Coordinate of radar pixel represents corner of pixel.
         Xref = int((i+xoffset)*xscale)
         Yref = int((j+yoffset)*yscale)
         LonArray[j][i], LatArray[j][i] = p(Xref,Yref,inverse=True)
         # Start with a rainfall field filled with nans:
         RArrayCML[j][i] = np.nan
         # The dimensions of X and Y should be one greater than those of C. 


   if PlotDataFieldRadarGrid=="yes":
      # Import row and column numbers of interpolation grid, i.e., connected to the Dutch radar grid:
      RowNr, ColNr = loadtxt(FileNameRowColNumbersRadarGrid, comments="#", delimiter=",", usecols=(0, 1), unpack=True, dtype ='int')
      RArray = np.array(RArrayCML)
      # Fill array with CML rainfall field values for row and column numbers provided in "FileNameRowColNumbersRadarGrid":
      for k in range(0,RowNr.shape[0]):
         RArray[RowNr[k]][ColNr[k]] = RArrayTemp[k]
      


   # Set data below LowestValue in [scale numbers] to "not available":
   if DoColorSetUnder!="yes":
      RArray[RArray < LowestValue] = np.nan
   if PlotKNMIRadar=="yes" and PlotDataFieldRadarGrid!="yes":
      # No data & missing data values are made "not available".
      RArray[RArray == nodata] = np.nan
      RArray[RArray == missingdata] = np.nan
      RArray = RArray * ConversionFactorKNMIRadar

   # Probably this still holds for this script & dataset:
   # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.pcolormesh.html:
   # The dimensions of X and Y should be one greater than those of C. Alternatively, X, Y and C may have equal dimensions, in which case the last row and column of C will be ignored.
   # Here the dimensions of X and Y are one greater than those of C. 



#############################################################
# 4. Read HDF5 files (ODIM HDF5 format) and obtain polygons.#
#############################################################

if PlotOPERARadar=="yes":
   DATAFIELD_NAME = DatasetNr + '/data1/data'
   FILE_NAME = OPERARadarInputFileName
   f = h5py.File(FILE_NAME, mode='r')
   # Read metadata:    
   Ncols = int(f['/where'].attrs['xsize'])
   Nrows = int(f['/where'].attrs['ysize'])

   ATTR_NAME = DatasetNr + '/what'
   zscale = f[ATTR_NAME].attrs['gain']
   zoffset = f[ATTR_NAME].attrs['offset']
   nodata = f[ATTR_NAME].attrs['nodata']
   undetect = f[ATTR_NAME].attrs['undetect']

   # Read data:
   dset = f[DATAFIELD_NAME]
   RArray = zoffset + zscale * dset[:]

   # Read file with coordinates OPERA radar grid:
   Grid = np.array(pd.read_csv("CoordinatesHDF5ODIMWGS84.dat", delimiter = " ", dtype="float",header=None))
   Xcoor = Grid[:,0]
   Ycoor = Grid[:,1]


   # Image coordinates to longitude & latitude in degrees:
   LonArray = [[np.nan for x in range(Ncols)] for x in range(Nrows)]
   LatArray = [[np.nan for x in range(Ncols)] for x in range(Nrows)]
   # Obtain image coordinates of surrounding pixels (center of grid cells):
   Nrow = Nrows
   Ncol = Ncols
   for j in range(0,Nrow):
      for i in range(0,Ncol):
         # Coordinate of pixel represents center of pixel.  
         LonArray[j][i] = Xcoor[i+j*Ncols]
         LatArray[j][i] = Ycoor[i+j*Ncols]

   # Set data below LowestValue in [scale numbers] to "not available":
   RArray[np.isnan(RArray)] = nodata
   if DoColorSetUnder!="yes":
      RArray[RArray < LowestValue] = np.nan
   # No data & undetect data values are made "not available".
   RArray[RArray == nodata] = np.nan
   RArray[RArray == undetect] = np.nan



##################
# 5. Make a plot.#
##################

# Map settings (e.g. projection and extent of area):
plt.rcParams["font.family"] = "serif"
plt.rcParams.update({'font.size': FontSizeLegendLabel}) 
plt.close('all')
transform = ccrs.PlateCarree()
fig = plt.figure(figsize=(16, 8))


# To use OpenStreetMap:
if TypeBackGroundMap=="OSM":
   request = cimgt.OSM()
   # Set map:
   ax = plt.axes(projection=request.crs)
   ax.set_extent(extent)
   ax.add_image(request, ValueRequest)
# To use Google Maps:
if TypeBackGroundMap=="GM":
   request = cimgt.GoogleTiles(style=style)
   # Set map:
   ax = plt.axes(projection=request.crs)
   ax.set_extent(extent)
   ax.add_image(request, ValueRequest)
# To use Natural Earth map:
if TypeBackGroundMap=="NE":
   # Map settings (e.g. projection and extent of area):
   if Projection=="yes":
      ax = plt.axes(projection=projection)
   else:
      ax = plt.axes(projection=transform)
   ax.set_extent(extent)


# Plot commercial microwave link paths:
if PlotCML=="yes" or PlotCMLTimeInterval=="yes":
   if (PlotValuesCMLPathAverage!="yes"):
      # Plot locations CMLs:
      for i in range(0,len(startlon)):
         xs = []
         ys = []
         xs.append(startlon[i])
         ys.append(startlat[i])
         xs.append(endlon[i])
         ys.append(endlat[i])
         plt.plot(xs,ys,linewidth=CMLLineWidth,alpha=alphaPath,c=ColorLinkPaths, zorder=zorderPath, transform=transform)
   

# Plot locations of rain gauges from network 1:
if PlotGaugeNetwork1=="yes":
   # Plot locations rain gauges as points:
   if PlotValuesGaugesNetwork1=="no":
      plt.scatter(lon_gaugesNetwork1, lat_gaugesNetwork1, s=SizeMarkerGaugesNetwork1, color=ColorGaugeNetwork1, zorder=zorderPoint, transform=transform, marker=MarkerGaugeNetwork1, alpha=alphaPoint)
      if LabelGaugeNetwork=="yes":
         # Plot legend:
         plt.scatter(LongitudeLabelGaugesNetwork1, LatitudeLabelGaugesNetwork1, s=SizeMarkerGaugesNetwork1, color=ColorGaugeNetwork1, zorder=zorderPoint, transform=transform, marker=MarkerGaugeNetwork1, alpha=alphaPoint)
         plt.text(LongitudeLabelGaugesNetwork1+LongitudeOffsetLabelNameAndNamesGaugesNetwork1, LatitudeLabelGaugesNetwork1+LatitudeOffsetLabelNameAndNamesGaugesNetwork1, LabelNameGaugesNetwork1, fontsize=FontSizeLabelNameGaugesNetwork1, transform=transform, color=ColorLabelNameGaugesNetwork1)
   if NamesGaugeNetwork=="yes":
      for i, label in enumerate(station_nameNetwork1):
         plt.text(lon_gaugesNetwork1[i]+LongitudeOffsetLabelNameAndNamesGaugesNetwork1, lat_gaugesNetwork1[i]+LatitudeOffsetLabelNameAndNamesGaugesNetwork1, label, zorder=zorderPoint, transform=transform, fontsize=FontSizeNamesGaugesNetwork1, color=ColorNamesGaugesNetwork1)


# Plot locations of rain gauges from network 2:
if PlotGaugeNetwork2=="yes":
   # Plot location hourly rain gauges as points:
   if PlotValuesGaugesNetwork2=="no":
      plt.scatter(lon_gaugesNetwork2, lat_gaugesNetwork2, s=SizeMarkerGaugesNetwork2, color=ColorGaugeNetwork2, zorder=zorderPoint, transform=transform, marker=MarkerGaugeNetwork2, alpha=alphaPoint)
      if LabelGaugeNetwork=="yes":
         # Plot legend:
         plt.scatter(LongitudeLabelGaugesNetwork2, LatitudeLabelGaugesNetwork2, s=SizeMarkerGaugesNetwork2, color=ColorGaugeNetwork2, zorder=zorderPoint, transform=transform, marker=MarkerGaugeNetwork2, alpha=alphaPoint)
         plt.text(LongitudeLabelGaugesNetwork2+LongitudeOffsetLabelNameAndNamesGaugesNetwork2, LatitudeLabelGaugesNetwork2+LatitudeOffsetLabelNameAndNamesGaugesNetwork2, LabelNameGaugesNetwork2, fontsize=FontSizeLabelNameGaugesNetwork2, transform=transform, color=ColorLabelNameGaugesNetwork2)
   if NamesGaugeNetwork=="yes":
      for i, label in enumerate(station_nameNetwork2):
         plt.text(lon_gaugesNetwork2[i]+LongitudeOffsetLabelNameAndNamesGaugesNetwork2, lat_gaugesNetwork2[i]+LatitudeOffsetLabelNameAndNamesGaugesNetwork2, label, zorder=zorderPoint, transform=transform, fontsize=FontSizeNamesGaugesNetwork2, color=ColorNamesGaugesNetwork2)


# Plot locations of rain gauges from network 3:
if PlotGaugeNetwork3=="yes":
   # Plot location hourly rain gauges as points:
   plt.scatter(lon_gaugesNetwork3, lat_gaugesNetwork3, s=SizeMarkerGaugesNetwork3, color=ColorGaugeNetwork3, zorder=zorderPoint, transform=transform, marker=MarkerGaugeNetwork3, alpha=alphaPoint)
   if LabelGaugeNetwork=="yes":
      # Plot legend:
      plt.scatter(LongitudeLabelGaugesNetwork3, LatitudeLabelGaugesNetwork3, s=SizeMarkerGaugesNetwork3, color=ColorGaugeNetwork3, zorder=zorderPoint, transform=transform, marker=MarkerGaugeNetwork3, alpha=alphaPoint)
      plt.text(LongitudeLabelGaugesNetwork3+LongitudeOffsetLabelNameAndNamesGaugesNetwork3, LatitudeLabelGaugesNetwork3+LatitudeOffsetLabelNameAndNamesGaugesNetwork3, LabelNameGaugesNetwork3, fontsize=FontSizeLabelNameGaugesNetwork3, transform=transform, color=ColorLabelNameGaugesNetwork3)
   if NamesGaugeNetwork=="yes":
      for i, label in enumerate(station_nameNetwork3):
         plt.text(lon_gaugesNetwork3[i]+LongitudeOffsetLabelNameAndNamesGaugesNetwork3, lat_gaugesNetwork3[i]+LatitudeOffsetLabelNameAndNamesGaugesNetwork3, label, zorder=zorderPoint, transform=transform, fontsize=FontSizeNamesGaugesNetwork3, color=ColorNamesGaugesNetwork3)
            


if PlotDataField=="yes" or PlotValuesGaugesNetwork1=="yes" or PlotValuesGaugesNetwork2=="yes" or PlotValuesCMLPathAverage=="yes" or PlotKNMIRadar=="yes" or PlotOPERARadar=="yes":
   # Obtain color scale:
   if ScaleType!="YellowRed" and ScaleType!="CbF":
      levels = levels
      cmap = mpl.colormaps.get_cmap(ScaleType)
   if ScaleType=="YellowRed":
      levels = levels
      cmap = mpl.colors.ListedColormap(['#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026'])
      colorSetOver = "black"
      colorSetUnder = "white"
   if ScaleType=="CbF":
      levels = levels
      cmap = mpl.colors.ListedColormap(['#DBEED3','#9CD5C4','#71B5C7','#858AC1','#A2569C','#96344E'])
      colorSetOver = "black"
      colorSetUnder = "white"
 
   norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)


   # Plot data field:
   if PlotDataField=="yes":
      if alpha<1:
         CS3 = plt.pcolormesh(np.asarray(LonArray),np.asarray(LatArray),RArray,alpha=alpha,cmap=cmap,norm=norm,linewidth=0.001,zorder=zorderMap,transform=transform) 
      if alpha==1:
         CS3 = plt.pcolormesh(np.asarray(LonArray),np.asarray(LatArray),RArray,alpha=alpha,cmap=cmap,norm=norm,zorder=zorderMap,transform=transform) 


   if PlotKNMIRadar=="yes":
      # Plot gridded KNMI radar data as colored polygons:
      CS3 = plt.pcolormesh(np.asarray(LonArray),np.asarray(LatArray),RArray,alpha=alpha,cmap=cmap,norm=norm,zorder=zorderMap,transform=transform) 
    

   if PlotOPERARadar=="yes":
      # Plot gridded OPERA radar data as colored polygons:
      CS3 = plt.pcolormesh(np.asarray(LonArray),np.asarray(LatArray),RArray,alpha=alpha,cmap=cmap,norm=norm,zorder=zorderMap,transform=transform,shading='nearest')


   # Plot point values for gauges:
   try:
      CS3 
   except NameError:
      if PlotGaugeNetwork1=="yes" and PlotValuesGaugesNetwork1=="yes":
         # Plot locations rain gauges as points:
         CS3 = plt.scatter(lon_gaugesNetwork1, lat_gaugesNetwork1, s=SizeMarkerGaugesNetwork1, c=GaugeValues1, cmap=cmap, zorder=zorderPoint, transform=transform, alpha=alpha, norm=norm, facecolors="none", edgecolors="black")
      if PlotGaugeNetwork2=="yes" and PlotValuesGaugesNetwork2=="yes":
         # Plot locations rain gauges as points:
         CS3 = plt.scatter(lon_gaugesNetwork2, lat_gaugesNetwork2, s=SizeMarkerGaugesNetwork2, c=GaugeValues2, cmap=cmap, zorder=zorderPoint, transform=transform, alpha=alpha, norm=norm, facecolors="none", edgecolors="black")
   else:
      if PlotGaugeNetwork1=="yes" and PlotValuesGaugesNetwork1=="yes":
         # Plot locations rain gauges as points:
         plt.scatter(lon_gaugesNetwork1, lat_gaugesNetwork1, s=SizeMarkerGaugesNetwork1, c=GaugeValues1, cmap=cmap, zorder=zorderPoint, transform=transform, alpha=alpha, norm=norm, facecolors="none", edgecolors="black")
      if PlotGaugeNetwork2=="yes" and PlotValuesGaugesNetwork2=="yes":
         # Plot locations rain gauges as points:
         plt.scatter(lon_gaugesNetwork2, lat_gaugesNetwork2, s=SizeMarkerGaugesNetwork2, c=GaugeValues2, cmap=cmap, zorder=zorderPoint, transform=transform, alpha=alpha, norm=norm, facecolors="none", edgecolors="black")
   

   # Plot path-average values for commercial microwave links:
   if PlotValuesCMLPathAverage=="yes":
      lags_norm = norm(RainfallDepthPath)
      # Plot locations CMLs:
      if DoColorSetUnder!="yes":
         for i in range(0,len(startlon)):
            if RainfallDepthPath[i] >= LowestValue:
               xs = []
               ys = []
               xs.append(startlon[i])
               ys.append(startlat[i])
               xs.append(endlon[i])
               ys.append(endlat[i])
               try:
                  CS3 
               except NameError:
                  # Plot point which is not visible. This is done to acquire the right color bar. Otherwise program will crash when "CS3.cmap.set_over" is used.
                  CS3 = plt.scatter(startlon[i],startlat[i], s=0.0001, c=RainfallDepthPath[i], cmap=cmap, zorder=zorderPath, transform=transform, alpha=alpha, norm=norm, facecolors="none")
                  plt.plot(xs,ys,linewidth=CMLLineWidth,alpha=alpha,c=cmap(lags_norm[i]), zorder=zorderPath, transform=transform)
               else:
                  CS3.cmap.set_over(colorSetOver, alpha=alpha)
                  plt.plot(xs,ys,linewidth=CMLLineWidth,alpha=alpha,c=cmap(lags_norm[i]), zorder=zorderPath, transform=transform)
      else:
         for i in range(0,len(startlon)):
            xs = []
            ys = []
            xs.append(startlon[i])
            ys.append(startlat[i])
            xs.append(endlon[i])
            ys.append(endlat[i])
            try:
               CS3 
            except NameError:
               # Plot point which is not visible. This is done to acquire the right color bar. Otherwise program will crash when "CS3.cmap.set_over" is used.
               CS3 = plt.scatter(startlon[i],startlat[i], s=0.0001, c=RainfallDepthPath[i], cmap=cmap, zorder=zorderPath, transform=transform, alpha=alpha, norm=norm, facecolors="none")
               plt.plot(xs,ys,linewidth=CMLLineWidth,alpha=alpha,c=cmap(lags_norm[i]), zorder=zorderPath, transform=transform)
            else:
               CS3.cmap.set_over(colorSetOver, alpha=alpha)
               CS3.cmap.set_under(colorSetUnder, alpha=alpha)
               plt.plot(xs,ys,linewidth=CMLLineWidth,alpha=alpha,c=cmap(lags_norm[i]), zorder=zorderPath, transform=transform)


   # Set highest class to chosen color "colorSetOver":
   CS3.cmap.set_over(colorSetOver, alpha=alpha)


   # Set proper location of label:
   if DoColorSetUnder=="yes":
      CS3.cmap.set_under(colorSetUnder, alpha=alpha)
      if ColorBar=="yes":
         # Plot color bar:
         font = mpl.font_manager.FontProperties(size=FontSizeLegend)         
         cbar = plt.colorbar(pad=DistanceBetweenMapAndColorBar,shrink=SizeLegend,extend="both")
         cbar.set_label(LabelName)
         text = ax.yaxis.label
         text.set_font_properties(font)
         cbar.ax.tick_params(labelsize=FontSizeLegend)         
   else:
      if ColorBar=="yes":
         # Plot color bar:
         font = mpl.font_manager.FontProperties(size=FontSizeLegend)         
         cbar = plt.colorbar(pad=DistanceBetweenMapAndColorBar,shrink=SizeLegend)
         cbar.set_label(LabelName)
         text = ax.yaxis.label
         text.set_font_properties(font)
         cbar.ax.tick_params(labelsize=FontSizeLegend)         


if TypeBackGroundMap=="NE":
   # Add natural earth features and borders in case of Natural Earth map:
   ax.add_feature(cartopy.feature.LAND, facecolor=ColorLand)
   ax.add_feature(cartopy.feature.OCEAN, facecolor=ColorOceanRiverLakes)
   ax.add_feature(cartopy.feature.LAKES, facecolor=ColorOceanRiverLakes, linewidth=0.00001,zorder=1)
   if DrawRivers=="yes":
      ax.add_feature(cartopy.feature.RIVERS, edgecolor=ColorOceanRiverLakes, linewidth=1.2, zorder=2)
   if DrawProvinces=="yes":
      ax.add_feature(cartopy.feature.STATES.with_scale("10m"), linewidth=0.3, zorder=2, edgecolor="gray")
   if DrawCountries=="yes":    
      ax.add_feature(cartopy.feature.BORDERS, linestyle="-", linewidth=LineWidthCoastLinesCountriesLakeLines, zorder=2, color=ColorCoastLinesCountriesLakeLines)
   if DrawCoastlines=="yes":
       ax.coastlines(resolution="10m", linewidth=LineWidthCoastLinesCountriesLakeLines, zorder=2, color=ColorCoastLinesCountriesLakeLines)
   if DrawLakelines=="yes":
      ax.add_feature(cartopy.feature.LAKES, edgecolor=ColorCoastLinesCountriesLakeLines, linewidth=LineWidthCoastLinesCountriesLakeLines, facecolor="none",zorder=2)

 
# Draw parallels & meridians: 
if DrawParallelsMeridians=="yes":
    if PlotParallelsMeridiansInFront=="yes":
       ax.set_axisbelow(False)
    ax.set_xticks(LongitudesMeridiansLevels, crs=transform)
    ax.set_xticklabels(LongitudesMeridiansLevels, fontsize=FontSizeTickLabels)
    ax.set_yticks(LatitudesParallelsLevels, crs=transform)
    ax.set_yticklabels(LatitudesParallelsLevels, fontsize=FontSizeTickLabels)
    ax.yaxis.tick_right()
    lon_formatter = cticker.LongitudeFormatter()
    lat_formatter = cticker.LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.grid(linewidth=LineWidth, color=ColorParallelsMeridians, alpha=alphaParallelsMeridians, linestyle=LineStyleParallelsMeridans)


# Plot North arrow: 
if DrawNorthArrow=="yes":
   ax.annotate('', xy=(LongitudeNorthArrowEnd, LatitudeNorthArrowEnd), xytext=(LongitudeNorthArrowStart,LatitudeNorthArrowStart),xycoords=transform._as_mpl_transform(ax),size=SizeNorthArrow,arrowprops=dict(facecolor="black",arrowstyle="fancy"))
   ax.text(LongitudeNorthArrowLabel, LatitudeNorthArrowLabel, "N", fontsize=FontSizeNorthArrowLabel, transform=transform)
   
   
# Plot extra text:
if PlotExtraText=="yes":
   ax.text(LonExtraText, LatExtraText, ExtraText, color=ColorExtraText, size=FontSizeExtraText, transform=transform)


# Plot title:
plt.title(TitlePlot,fontsize=FontSizeTitle) 


# Plot scale bar:
if DrawScaleBar=="yes":
   scale_bar(ax, LengthScaleBar)  


# Save figure:
plt.savefig(OutputFileName, bbox_inches = "tight", dpi = dpi)

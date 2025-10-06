#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python configuration script.
# Name: ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py
# Part of MapRAINLINK: https://github.com/overeem11/MapRAINLINK
#
#
## Version 1.12
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
# Settings (of parameters) needed for running "VisualizeCMLsGaugesRadars_CartopyOSM_GM.py",
# the Python script to visualize rain gauge, radar, and commercial microwave link locations and their rainfall estimates on a map.
# Mainly graphical settings for map written to output file.
# Available color scales: https://matplotlib.org/examples/color/colormaps_reference.html for parameters such as ColorLand.
# In publications you could state that Python package Cartopy has been used:
# @MANUAL{Cartopy,
# author = {{Met Office}},
# title = {Cartopy: a cartographic python library with a Matplotlib interface},
# year = {2010 - 2015},
# address = {Exeter, Devon },
# url = {https://scitools.org.uk/cartopy}
# }


#####################
# 1. Background map.#
#####################
#              Please use proper copyright in figure captions:
#              - For OSM: OpenStreetMap; https://www.openstreetmap.org/copyright or map data ©OpenStreetMap contributors.
#              - For GM: Map data ©2022 Google
#              - For NE: Maps made with Natural Earth. Free vector and raster map data ©naturalearthdata.com.

# 1.1 Type of background map ("OSM" for OpenStreetMap, "GM" for GoogleMaps, "NE" for Natural Earth.
TypeBackGroundMap = "NE"

# 1.2 Style of background map (only for GoogleMaps): "satellite", "street".
style = "street"

# 1.3 Value for request: determines resolution of background map / size of text (e.g. city names on map; for OpenStreetMap and GoogleMaps only).
ValueRequest = 9
# When zooming in use larger values. Too large values may result in an error. Then use lower values.

# 1.4 Map options for Natural Earth background map.
ColorLand = "lightgrey"                  	# Color of land surface.
ColorOceanRiverLakes = "lightblue"      	# Color of oceans, seas, rivers and lakes.
DrawProvinces = "no"     			# "yes" for drawing borders of departments/provinces. 
DrawRivers = "no"        			# "yes" for drawing rivers.
Projection = "yes"                               # "yes" to use projection with EPSG code given below:
projection = ccrs.epsg(3035)     		# epsg:3035 = ETRS89 / ETRS-LAEA (suited for Europe). See https://epsg.io/ for EPSG codes per region.

# 1.5 Map options for all background maps (note that borders can deviate a bit from the OpenStreetMap or GoogleMaps borders).
DrawCoastlines = "yes"    			# "yes" for drawing coastlines. Note that coastlines and country borders seem less accurate compared to OpenStreetMap and GoogleMap.
DrawCountries = "yes"	  			# "yes" for drawing country borders.
DrawLakelines = "yes"     			# "yes" for drawing lake lines.
LineWidthCoastLinesCountriesLakeLines = 0.3     # 0.3 works wel for European scale, 1.2 works well for northwestern Europe scale.
ColorCoastLinesCountriesLakeLines = "black"     # "black" can work well for European scale, "dimgray" can work well for northwestern Europe scale.



##############################
# 2. General settings of map.#
##############################

# 2.1 Switches for plotting.
DrawNorthArrow = "yes"           		# "yes" for drawing a north arrow.
DrawParallelsMeridians = "no"   		# "yes" for drawing parallels & meridians.
DrawScaleBar = "yes"		             	# "yes" for drawing a scale bar.
# Note that EPSG code will not work when parallels & meridians are plotted (this can only happen for Natural Earth background map). To prevent this:
if Projection=="yes" and TypeBackGroundMap=="NE":
   DrawParallelsMeridians = "no"


# 2.2 Output filename (map).
OutputFileName = "Netherlands.jpg"		# Name of output file with map. Can at least be "jpg", "pdf" or "png", which automatically chooses the file format.


# 2.3 Miscellaneous settings.
alpha = 1					# The alpha blending value, between 0 (transparent) and 1 (opaque) for plotting rainfall maps, values for path-averages and point values.
alphaPath = 0.4				        # The alpha blending value, between 0 (transparent) and 1 (opaque) for plotting locations of microwave links (but not used for their values).
						# Note that in case of a full-duplex link, i.e., 2 sub-links per link path, both are plotted and the value from only one of them will be visible
						# in case of an alpha value of 1, or their values are probably not distinguishable in case of an alpha value below 1.
alphaPoint = 1					# The alpha blending value, between 0 (transparent) and 1 (opaque) for plotting locations of rain gauges (but not used for their values).
dpi = 600 					# Resolution of output figure in dpi.
Coordinates = "[3.274, 7.273, 50.58, 53.593]"	# Plotting area: minimum longitude, maximum longitude, minimum latitude, maximum latitude. 
# http://www.openstreetmap.org/export can be used to select an appropriate area (bounding box). 
#Coordinates = "[79.64, 82.01, 5.86, 9.89]"	# Plotting area for Sri Lanka.
#Coordinates = "[-10, 30, 32.7, 73]"            # Plotting area for Europe.
extent = list(map(float, Coordinates.strip('[]').split(',')))	# Automatically obtain list of floats of values for bounding box.
# The order of plotting is determined by the zorder attribute. This can be chosen for part of the plot. Example: a higher zorder value for points means that these are plotted in front of rainfall maps,
# so that the underlying values from the rainfall map are not visible at the gauge locations, provided that alphaPoint = 1.
# It seems that if zorder values are the same, the (type of) variable which is plotted at the end, is plotted in front of all the other (type of) variables.
zorderMap = 2					# The zorder value for plotting rainfall maps.
zorderPath = 3				        # The zorder value for plotting locations of microwave links (or their values).
zorderPoint = 3					# The zorder value for plotting locations of rain gauges (or their values) and labels.


# 2.4 Settings for parallels & meridians.
alphaParallelsMeridians = 0.5			# The alpha blending value, between 0 (transparent) and 1 (opaque) for parallels & meridians.
ColorParallelsMeridians = "black"		# Color of parallels & meridians.
FontSizeTickLabels = 15				# Font size of x and y tick labels.
LatitudesParallels = "[51, 52, 53]"		# Latitudes of parallels to be plotted.
LatitudesParallelsLevels = list(map(float, LatitudesParallels.strip('[]').split(',')))		# Automatically obtain list of floats of values for latitudes of meridians.
LineStyleParallelsMeridans = "--"		# The line style for parallels & meridians.
LineWidth = 2					# Line width of parallels & meridians.
LongitudesMeridians = "[4, 5, 6, 7]"		# Longitudes of meridians to be plotted.
LongitudesMeridiansLevels = list(map(float, LongitudesMeridians.strip('[]').split(',')))	# Automatically obtain list of floats of values for longitudes of meridians.
PlotParallelsMeridiansInFront = "yes"           # "yes" for drawing parallels & meridians in front of all other visualizations.


# 2.5 Settings for drawing north arrow.
FontSizeNorthArrowLabel = 25			# Font size of north arrow.
LatitudeNorthArrowEnd = 51.1			# Latitude of end of north arrow.
LongitudeNorthArrowEnd = 6.9			# Longitude of end of north arrow.
LatitudeNorthArrowStart = 50.7			# Latitude of start of north arrow.
LongitudeNorthArrowStart = 6.9			# Longitude of start of north arrow.
LatitudeNorthArrowLabel = 51.1			# Latitude of north arrow label ("N").
LongitudeNorthArrowLabel = 6.95			# Longitude of north arrow label ("N").
SizeNorthArrow = 25				# Size of north arrow.


# 2.6. Settings for scale bar (a visual indication of distance and feature size on the map).
FontSizeScaleBar = 19				# Font size of accompanying text "km" (= kilometers) and values of scale bar.
LengthScaleBar = 50				# Length of scale bar in kilometers.
LineWidthScaleBar = 3				# Line width of scale bar.
XPositionScaleBar = 0.4				# Horizontal position of scale bar. Should be value between 0-1.
YPositionScaleBar = 0.13			# Vertical position of scale bar. Should be value between 0-1.



################################################
# 3. Label, legend (color bar), text and title.# 
################################################

# 3.1 Switches for plotting.
ColorBar = "yes"				# "yes" for drawing a color bar (legend).
DoColorSetUnder = "no"				# "yes" for drawing a class below the lowest value of the supplied scale "LevelsScale".
PlotExtraText = "yes"				# "yes" to plot extra text.


# 3.2 Options for title.
FontSizeTitle = 29				# Font size of title.
TitlePlot = "The Netherlands"		  	# Title of map.


# 3.3 Options for color bar (legend).
DistanceBetweenMapAndColorBar = 0.02		# Increase this value to make distance between map and color bar larger. Different values may be needed for different sizes of plotted areas.
						# These values could be tested as starting point: use larger value, e.g., 0.12, in case DrawParallelsMeridians = "yes", otherwise 0.02.
colorSetOver = "lightgray"			# Color for highest class.
colorSetUnder = "indigo" 			# Color for lowest class.
FontSizeLegend = 31				# Font size of legend.
FontSizeLegendLabel = 29			# Font size of legend label.
LabelName = "15-min rainfall depth (mm)"	# Label name of legend. Also Latex mathematical notation is accepted, e.g.: "Rainfall intensity (mm h$^{-1}$)"
LevelsScale = "[0.1,1.6,3.1,4.6,6.1,7.6,9.1]"	# Values for legend scale. Simply add or remove values to increase or decrease the number of classes 
                                                # if "ScaleType" is not equal to "YellowRed" and "CbF".
# LevelsScale = "[1,10,20,30,40,50,60]"	# Use this scale for plotting radar file "RAD_NL25_RAC_MFBS_24H_201805300800_NL.h5" with 24-h precipitation accumulation.
levels = list(map(float, LevelsScale.strip('[]').split(',')))	# Automatically obtain list of floats of values for legend scale.
LowestValue = float(LevelsScale.strip('[]').split(',')[0])      # Automatically obtain lowest value of legend scale.
ScaleType = "CbF"				# Color scheme/scale for legend, often colorblind friendly.
# Simply choose a colormap for the color bar from https://matplotlib.org/stable/gallery/color/colormap_reference.html?highlight=colormaps
# "Blues" can work out well with colorSetOver = "indigo" & colorSetUnder = "lightgray".
# 4 classes for "Blues" does not work well, then you get two white classes. In general, 4 classes seems not to work with get_cmap, so use at least 5 classes.
# Lowest class, e.g. 0.3 - 1.0 mm, includes the 0.3 values. Values below 1 are plotted if "DoColorSetUnder = "yes"" in the color specified by "colorSetUnder".
# The highest value of the highest class is plotted in the color specified by colorSetOver, i.e. as of 32.0 if this is the last number in "LevelsScale".
# So for each class its lowest value, i.e. the lowest value at the tick mark, belongs to that class and is plotted in that color, 
# whereas its highest value (the next tick mark) belongs to the next class.
# In addition, you can choose "CbF" or "YellowRed" for ScaleType, which are developed in case "LevelsScale" contains 7 values, i.e., 6 classes and a class for colorSetOver
#  and optionally a class for colorSetUnder. This does not work for a larger number of classes. Note that "colorSetOver = "black"" and "colorSetUnder = "white"" for these two color schemes.
SizeLegend = 0.8				# Size of legend scale.


# 3.4 Options for extra text.
ColorExtraText = "black"			# Color of extra text.
ExtraText = "End time: 201109102045"		# Extra text.
FontSizeExtraText = 15				# Font size of extra text.
LatExtraText = 50.65				# Latitude of extra text. 
LonExtraText = 3.4				# Longitude of extra text.



##################################################################################################################
# 4. Commercial microwave link path-averaged rainfall - or path values for any other sensor, variable or dataset.#
##################################################################################################################

# 4.1 Switches for plotting.
PlotCML = "no"					# "yes" for plotting link paths from file "FileNameCMLLocations".
PlotCMLTimeInterval = "yes"			# "yes" for plotting link paths from file "FileNameCMLLocationsTimeInterval", i.e. those from the specific time interval. 
# If "PlotCMLTimeInterval" is "yes" the link paths from file "FileNameCMLLocations" are not plotted:
if PlotCMLTimeInterval=="yes":
   PlotCML = "no"
PlotValuesCMLPathAverage = "no"			# "yes" for plotting path-average rainfall values according to color scale. Supersedes "ColorLinkPaths".
# If "PlotCMLTimeInterval" is "no", path-average rainfall values cannot be plotted:
if PlotCMLTimeInterval!="yes":
   PlotValuesCMLPathAverage = "no"


# 4.2 Input filenames.
FileNameCMLLocations = "CMLLocations_SriLanka_RmeanAvailable_MinMax.dat"	# File name for file with commercial microwave link locations.
FileNameCMLLocationsTimeInterval  = "linkdata_201109102045.dat" 	# File name for file with commercial microwave link locations for specific time interval.


# 4.3 Other options.
ColorLinkPaths = "black"			# Color of link paths.	
ConversionFactorCMLPathAverage = 1		# Conversion factor in case rainfall accumulations or intensities need to be converted. For path-average rainfall.
						# E.g. use 4 to convert 15-min rainfall accumulations to intensities in millimeters per hour.
						# 0.25 is the standard conversion needed for RAINLINK mean path-averaged rainfall intensity to obtain 15-min rainfall accumulations.						
						# Note that when "Run.R" script of RAINLINK is used to store these data, these are already converted to 15-min rainfall accumulations.
						# Hence, when files created by "Run.R" are used, the conversion has to be 1.
CMLLineWidth = 1				# Line width of link paths.



##########################################################################################################
# 5. Commercial microwave link rainfall maps - or data fields from any other sensor, variable or dataset.#
##########################################################################################################

# 5.1 Switch for plotting.
PlotDataField = "yes"				# "yes" for plotting commercial microwave link rainfall maps.


# 5.2 Other options.
ConversionFactorCMLDataField = 0.25		# Conversion factor in case rainfall accumulations or intensities need to be converted. For interpolated rainfall data.
						# E.g. use 4 to convert 15-min rainfall accumulations to intensities in millimeters per hour.
						# 0.25 is the standard conversion needed for RAINLINK interpolated fields of rainfall intensities to obtain 15-min rainfall accumulations.
GridRes = 0.02					# Resolution of grid in degrees (e.g. 0.02). Has to match resolution of grid in "InterpolationGrid".
InterpolationGrid = "InterpolationGrid.dat"	# Filename of file with grid coordinates.
# The coordinates for this field represent the middle of the grid cell and the resolution is on a regular grid with resolution "GridRes" in degrees.
# So the resolution in degrees is fixed. 
# The grid on which each field of values is presented is, e.g., a 0.02x0.02 degrees lat./lon. (Cylindrical Equal Distance) array of points.
InputFileNameDataField = "linkmap_201109102045.dat"	# Filename of file with interpolated rainfall map data.


# 5.3 Plot interpolated CML data on the Dutch radar grid, i.e., the default "InterpolationGrid.dat" provided by RAINLINK and used in https://doi.org/10.5194/amt-9-2425-2016.
# This actually reads a radar file, which is considered as template, where the radar input filename is given below "7. KNMI-HDF5 gridded radar precipitation data.".
# Note that this option overrules the plotting of radar data at "7. KNMI-HDF5 gridded radar precipitation data.".
# This option sets the radar data to 0 and fills the radar grid cells with the CML rainfall map values for the row and column numbers provided in file "FileNameRowColNumbersRadarGrid".
FileNameRowColNumbersRadarGrid = "InterpolationGridRowColNrsRadarGrid.dat"
PlotDataFieldRadarGrid = "yes"
if PlotDataField=="no" and PlotDataFieldRadarGrid=="yes":
   print("Error! Set PlotDataField to yes in order to plot the interpolated CML data on the Dutch radar grid!")
   
   


############################################################################################################
# 6. Rain gauges or point values for any other (opportunistic) sensor and variable (e.g., radar locations).#
############################################################################################################
# Locations for 3 networks and values for 2 networks can be plotted.

# 6.1 Switches for plotting.
LabelGaugeNetwork = "yes"			# "yes" to plot legend for markers of rain gauge networks. It could be wise to not plot it in case "PlotValuesGauges" is "yes".
NamesGaugeNetwork = "yes"			# "yes" to plot (abbreviated) names/IDs of rain gauge locations near markers.
PlotGaugeNetwork1 = "no"			# "yes" for plotting locations of gauge network 1 by markers (also needed to plot label and names).
PlotGaugeNetwork2 = "no"			# "yes" for plotting locations of gauge network 2 by markers (also needed to plot label and names).
PlotGaugeNetwork3 = "no"			# "yes" for plotting locations of gauge network 3 by markers (also needed to plot label and names).
PlotValuesGaugesNetwork1 = "no"        		# "yes" for plotting point rainfall values according to color scale for gauge network 1. These are plotted as dots with size "SizeMarkerGauges1".
PlotValuesGaugesNetwork2 = "no"        		# "yes" for plotting point rainfall values according to color scale for gauge network 2. These are plotted as dots with size "SizeMarkerGauges2".
# For the third network no rainfall values can be plotted.
# In case no rain gauge networks need to be plotted, set "PlotValuesGaugesNetwork1" and/or "PlotValuesGaugesNetwork2" to "no".
if PlotGaugeNetwork1=="no":
   PlotValuesGaugesNetwork1 = "no"
if PlotGaugeNetwork2=="no":
   PlotValuesGaugesNetwork2 = "no"


# 6.2 Input filenames and conversion factors.
ConversionFactorGaugesNetwork1 = 1                     # Conversion factor in case rainfall accumulations or intensities need to be converted. For gauge network 1.
ConversionFactorGaugesNetwork2 = 1                     # Conversion factor in case rainfall accumulations or intensities need to be converted. For gauge network 2.
						# E.g. use 12 to convert from 5-minute accumulations to millimeters per hour.
FileNameGaugeLocationsNetwork1 = "StationMetadataNetwork1.dat"	# Filename with gauge locations for network 1.
FileNameGaugeLocationsNetwork2 = "StationMetadataNetwork2.dat"	# Filename with gauge locations for network 2.
FileNameGaugeLocationsNetwork3 = "StationMetadataNetwork3.dat"	# Filename with gauge locations for network 3.



# 6.3 Settings for markers.
ColorGaugeNetwork1 = "red"			# Color of marker for gauge network 1.
ColorGaugeNetwork2 = "purple"			# Color of marker for gauge network 2.
ColorGaugeNetwork3 = "darkgreen"		# Color of marker for gauge network 3.
MarkerGaugeNetwork1 = "^"			# Symbol for marker of gauge network 1.
MarkerGaugeNetwork2 = "D"			# Symbol for marker of gauge network 2.
MarkerGaugeNetwork3 = "o"			# Symbol for marker of gauge network 3. Note that not every symbol will work, at least "+" does not work.
# See https://matplotlib.org/stable/api/markers_api.html for a list of possible markers for the gauge locations.
SizeMarkerGaugesNetwork1 = 100				# Size of marker of gauge network 1.
SizeMarkerGaugesNetwork2 = 100				# Size of marker of gauge network 2.
SizeMarkerGaugesNetwork3 = 100				# Size of marker of gauge network 3.
LineWidthMarkerGauges = 0.1                             # Line width of point values (line width of outer circle in case of a circle/dot).


# 6.4 Settings for label name of markers and for plotting (abbreviated) names/IDs near markers.
ColorLabelNameGaugesNetwork1 = "black"		# Color of label name for gauge network 1.
ColorLabelNameGaugesNetwork2 = "black"		# Color of label name for gauge network 2.
ColorLabelNameGaugesNetwork3 = "black"		# Color of label name for gauge network 3.
ColorNamesGaugesNetwork1 = "black"              # Color of (abbreviated) names/IDs of rain gauge locations near locations for network 1.
ColorNamesGaugesNetwork2 = "black"              # Color of (abbreviated) names/IDs of rain gauge locations near locations for network 2.
ColorNamesGaugesNetwork3 = "black"              # Color of (abbreviated) names/IDs of rain gauge locations near locations for network 3.
FontSizeLabelNameGaugesNetwork1 = 15		# Font size of label name for gauge network 1.
FontSizeLabelNameGaugesNetwork2 = 15		# Font size of label name for gauge network 2.
FontSizeLabelNameGaugesNetwork3 = 15		# Font size of label name for gauge network 3.
FontSizeNamesGaugesNetwork1 = 15		# Font size of (abbreviated) names/IDs for gauge network 1.
FontSizeNamesGaugesNetwork2 = 15		# Font size of (abbreviated) names/IDs for gauge network 2.
FontSizeNamesGaugesNetwork3 = 15		# Font size of (abbreviated) names/IDs for gauge network 3.
LabelNameGaugesNetwork1 = "Radar locations"	# Label name for gauge network 1.
LabelNameGaugesNetwork2 = "Network 2"		# Label name for gauge network 2.
LabelNameGaugesNetwork3 = "Network 3"		# Label name for gauge network 3.
LatitudeLabelGaugesNetwork1 = 53.4		# Longitude of label for gauge network 1 (degrees).
LatitudeLabelGaugesNetwork2 = 53.2		# Longitude of label for gauge network 2 (degrees).
LatitudeLabelGaugesNetwork3 = 53.0		# Longitude of label for gauge network 3 (degrees).
LatitudeOffsetLabelNameAndNamesGaugesNetwork1 = 0.06	# Offset for longitude for plotting of label name and names of gauges for network 1.
LatitudeOffsetLabelNameAndNamesGaugesNetwork2 = 0.06	# Offset for longitude for plotting of label name and names of gauges for network 2.
LatitudeOffsetLabelNameAndNamesGaugesNetwork3 = 0.06	# Offset for longitude for plotting of label name and names of gauges for network 3.
LongitudeLabelGaugesNetwork1 = 3.4		# Latitude of label for gauge network 1 (degrees).
LongitudeLabelGaugesNetwork2 = 3.4		# Latitude of label for gauge network 2 (degrees).
LongitudeLabelGaugesNetwork3 = 3.4		# Latitude of label for gauge network 3 (degrees).
LongitudeOffsetLabelNameAndNamesGaugesNetwork1 = -0.028	# Offset for latitude for plotting of label name and names of gauges for network 1.
LongitudeOffsetLabelNameAndNamesGaugesNetwork2 = -0.028	# Offset for latitude for plotting of label name and names of gauges for network 2.
LongitudeOffsetLabelNameAndNamesGaugesNetwork3 = -0.028	# Offset for latitude for plotting of label name and names of gauges for network 3.
NumberOfCharactersNamesNetwork1 = 100           # Number of characters of names/IDs of rain gauge locations to plot for network 1. Set to a large value to plot entire names/IDs. 
NumberOfCharactersNamesNetwork2 = 4             # Number of characters of names/IDs of rain gauge locations to plot for network 2. Set to a large value to plot entire names/IDs.
NumberOfCharactersNamesNetwork3 = 4             # Number of characters of names/IDs of rain gauge locations to plot for network 3. Set to a large value to plot entire names/IDs.



#################################################
# 7. KNMI-HDF5 gridded radar precipitation data.#
#################################################

# 7.1 Switch for plotting.
PlotKNMIRadar = "no"				# "yes" for plotting 2D KNMI radar image.
# Prevent plotting of data field of microwave link rainfall maps in case KNMI radar data are plotted:
if PlotDataField=="yes" and PlotKNMIRadar=="yes":
   PlotDataField = "no"

# 7.2 Other options.
ConversionFactorKNMIRadar = 1			# Conversion factor in case rainfall accumulations or intensities need to be converted.
						# E.g. use 12 to convert 5-min rainfall accumulations to intensities in millimeters per hour.
                                                # Note that the KNMI precipitation radar data provide accumulations in millimeters.
DATAFIELD_NAME = "/image1/image_data"           # Path of radar image to be plotted.
DATAFIELD_NAMECAL = "/image1/calibration"       # Location of metadata to convert precipitation to millimeters and to obtain "nodata" & "undetect" values.
KNMIRadarInputFileName = "RAD_NL25_RAC_MFBS_24H_201805300800_NL.h5"	# Filename of file with gridded radar rainfall data.



#######################################################
# 8. ODIM-HDF5 gridded OPERA radar precipitation data.#
#######################################################
# OPERA radar data: https://www.eumetnet.eu/activities/observations-programme/current-activities/opera/

# 8.1 Switch for plotting.
PlotOPERARadar = "no"				# "yes" for plotting 2D OPERA radar image.
# Prevent plotting of data field of microwave link rainfall maps and KNMI radar data in case OPERA radar data are plotted:
if PlotDataField=="yes" and PlotOPERARadar=="yes":
   PlotDataField = "no"
if PlotKNMIRadar=="yes" and PlotOPERARadar=="yes":
   PlotKNMIRadar = "no"

# 8.2 Other options.
ConversionFactorOPERARadar = 1			# Conversion factor in case rainfall accumulations or intensities need to be converted.
						# E.g. use 0.25 to convert 15-min rainfall intensities (millimeters per hour) to 15-min rainfall accumulations.
                                                # Note that the OPERA precipitation radar data provide accumulations in millimeters in case of hourly accumulations or the 1-h and 24-h
                                                # EURADCLIM datasets, but OPERA also provides a product of instantaneous surface rainfall intensities.
DatasetNr = "/dataset1"           		# Path of radar image to be plotted. 
OPERARadarInputFileName = "RAD_OPERA_24H_RAINFALL_ACCUMULATION_201305311400.h5"	# Filename of file with gridded radar rainfall data.



############################################
# 9. Load function to plot kilometer ruler.#
############################################

if DrawScaleBar=="yes":
   # Function for plotting scale bar (taken from https://stackoverflow.com/questions/32333870/how-can-i-show-a-km-ruler-on-a-cartopy-matplotlib-plot)
   def scale_bar(ax, length=None, location=(XPositionScaleBar, YPositionScaleBar), linewidth=LineWidthScaleBar):
       """
       ax is the axes to draw the scalebar on.
       length is the length of the scalebar in km.
       location is center of the scalebar in axis coordinates.
       (ie. 0.5 is the middle of the plot)
       linewidth is the thickness of the scalebar.
       """
       #Get the limits of the axis in lat long
       llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
       #Make tmc horizontally centred on the middle of the map,
       #vertically at scale bar location
       sbllx = (llx1 + llx0) / 2
       sblly = lly0 + (lly1 - lly0) * location[1]
       tmc = ccrs.TransverseMercator(sbllx, sblly, approx=True)  # Aart Overeem (KNMI): added ", approx=True" because of warning.
       #Get the extent of the plotted area in coordinates in metres
       x0, x1, y0, y1 = ax.get_extent(tmc)
       #Turn the specified scalebar location into coordinates in metres
       sbx = x0 + (x1 - x0) * location[0]
       sby = y0 + (y1 - y0) * location[1]

       #Calculate a scale bar length if none has been given
       #(Theres probably a more pythonic way of rounding the number but this works)
       if not length: 
          length = (x1 - x0) / 5000 #in km
          ndim = int(np.floor(np.log10(length))) #number of digits in number
          length = round(length, -ndim) #round to 1sf
          #Returns numbers starting with the list
          def scale_number(x):
             if str(x)[0] in ['1', '2', '5']: return int(x)        
             else: return scale_number(x - 10 ** ndim)
          length = scale_number(length) 

       #Generate the x coordinate for the ends of the scalebar
       bar_xs = [sbx - length * 500, sbx + length * 500]
       #Plot the scalebar
       ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
       #Plot the scalebar label
       ax.text(sbx, sby, str(length) + ' km', transform=tmc,
       horizontalalignment='center', verticalalignment='bottom', size=FontSizeScaleBar)


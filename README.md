# MapRAINLINK

# Introduction
Python script to visualize rain gauge, radar, and commercial microwave link (CML) locations and their rainfall estimates on a map. This Python script is more widely applicable. It can be used to plot lines and/or points (and their values) on a map, not necessarily CML paths and rain gauge locations, and not necessarily rainfall intensities or accumulations. It can also be used to plot other data fields. The maps are highly customizable. For instance, the background map can be chosen from GoogleMaps, Natural Earth, or OpenStreetMap. Plotting a map typically takes 5 up to 15 seconds.

This script is meant as a replacement of the visualization by the open-source R package RAINLINK for Dutch radar data and commercial microwave link data. RAINLINK (https://github.com/overeem11/RAINLINK) is a retrieval algorithm for rainfall mapping from microwave links in a cellular communication network. MapRAINLINK works with output data from RAINLINK, irrespective of the RAINLINK version. MapRAINLINK can be of general use for opportunistic sensing data, such as those from CMLs and personal weather stations (PWSs). Note that some Python libraries need to be installed: cartopy, copy, h5py, matplotlib, numpy, os, pandas, pyproj, sys.

# How to start?
Just clone this repository or simply copy all files to a local directory. With the default settings, the locations of CML paths and the interpolated CML rainfall estimates are plotted on a map of the Netherlands for one 15-min time interval. The interpolated values have been obtained by running RAINLINK with the Dutch sample dataset, which is part of RAINLINK. The map "Netherlands.jpg" is made by running the vizualisation script "VisualizeCMLsGaugesRadars_CartopyOSM_GM.py" using Python version 3:
```
python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py
```
<img src="Netherlands.jpg" alt="drawing" width="500"/>

Map made with Natural Earth. Free vector and raster map data &copy naturalearthdata.com. Note that this reproduces Figure 5 (left) from https://doi.org/10.5194/amt-9-2425-2016 (but that was obtained with a different plotting package and different colors and classes).

# Detailed explanation of plotting
All the plotting options and input and output files need to be specified in the configuration script "ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py", although it is also possible to provide additional arguments on the command line, which will overrule the settings in the configuration script. These interpolated CML values from the above example are output from RAINLINK and are stored in file "linkmap_201109102045.dat". This file has been obtained by running RAINLINK with the built-in 2-day sample dataset. These are interpolated values at the middle of the grid cells of the Dutch radar grid ("InterpolationGrid.dat" from https://github.com/overeem11/RAINLINK). In order to plot these, "PlotDataFieldRadarGrid" must be set to "yes" in the configuration script. First, the radar file "RAD_NL25_RAC_MFBS_24H_201805300800_NL.h5" is read, so that the polygons for the plotting can be constructed ("KNMIRadarInputFileName" in the configuration script), where the associated radar values are set to 0 in the Python script. The file "InterpolationGridRowColNrsRadarGrid.dat" ("FileNameRowColNumbersRadarGrid" in the configuration script) is read to obtain the column and row numbers of the radar grid for the CML grid cells. Next, these grid cells are filled with the interpolated CML rainfall values. 

In case path-average rainfall needs to be visualized, the RAINLINK output from file "linkdata_201109102045.dat" is used, and "PlotCMLTimeInterval" and "PlotValuesCMLPathAverage" must be set to "yes" in the configuration script. "PlotDataField" may be set not equal to "yes" to avoid that the CML rainfall field and CML path averages are plotted together. The file "StationMetadataNetwork1.dat" contains locations of the two C-band ground-based weather radars in the Netherlands. These locations can be visualized by setting "PlotGaugeNetwork1" to "yes" in the configuration script. If also "PlotValuesGaugesNetwork1" is set to "yes", the associated rainfall values are visualized (this is meant as an illustration of plotting rainfall values for gauge locations, the radar itself obviously does not observe point values). More information on usage is provided in the configuration and visualization scripts. Note that the first line in a file with locations (e.g., from CMLs, radars, or rain gauges) is not used, because it is assumed that the first line contains the header with variable names. Also note that lines with # at the beginning are skipped (so make sure that the first line does not contain #; this is automatically the case for RAINLINK output with CML locations).

# Plotting for other regions or interpolation grids
For other regions or interpolation grids, "PlotDataFieldRadarGrid" must be made not equal to "yes", and a grid with "fixed changes in degrees", e.g. a 0.02 degree grid in latitude and longitude, must be used (defined by "InterpolationGrid" in the configuration script). Such an interpolation grid can be constructed by following the recipe below, and must be used with RAINLINK to obtain the interpolated data at this grid. Next, it can be used to visualize CML rainfall maps. Here, an example is provided for Sri Lanka, which reproduces Figure 5(e) from https://doi.org/10.1088/1748-9326/ac0fa6 (but that was obtained with a different plotting package), the map "SriLanka.jpg":
```
python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM_SriLanka.py
```
<img src="SriLanka.jpg" alt="drawing" width="500"/>
&copy OpenStreetMap contributors 2025. Distributed under the Open Data Commons Open Database License (ODbL) v1.0. The CML interpolated values are obtained from file "linkmap_201910171615_SriLanka.dat", which is "linkmap_201910171615.dat" in the publicly availabe dataset https://doi.org/10.4121/14166539.v2. A recipe to obtain an interpolation grid for RAINLINK and one for the plotting with "VisualizeCMLsGaugesRadars_CartopyOSM_GM.py" for Sri Lanka is given below:

1. Start R and run this script to obtain the interpolation grid for Sri Lanka:
```
source("MakeInterpolationGridSriLanka.R").
```
3. The output "InterpolationGrid_SriLanka.dat" has to be used by RAINLINK to obtain interpolated rainfall maps (3.5 months of these rainfall maps can be obtained from https://doi.org/10.4121/14166539.v2). Since the CML data are not publicly available, this step cannot be reproduced. Note that this step is not needed to actually obtain an interpolation grid for plotting.
4. After running RAINLINK, store the following output of RAINLINK in the file "CMLLocations_SriLanka_RmeanAvailable_MinMax.dat", by pasting this in the R shell (again this step cannot be reproduced and the file "CMLLocations_SriLanka_RmeanAvailable_MinMax.dat" is not publicly available because of restrictions to share the CML metadata):
```
Rmean <- RainRetrievalMinMaxRSL(Aa=Aa,alpha=alpha,Data=DataOutlierFiltered,kRPowerLawDataH=kRPowerLawDataH,kRPowerLawDataV=kRPowerLawDataV,PmaxCor=Pcor$PmaxCor,PminCor=Pcor$PminCor,Pref=Pref)

DataOutlierFiltered <- OutlierFilterMinMaxRSL(Data=DataPreprocessed,F=WetDry$F,FilterThreshold=FilterThreshold)

cond <- which(Rmean>=0)

dataf <- unique(data.frame(cbind(DataOutlierFiltered$XStart[cond],DataOutlierFiltered$YStart[cond],DataOutlierFiltered$XEnd[cond],DataOutlierFiltered$YEnd[cond])))

write.table(dataf,"CMLLocations_SriLanka_RmeanAvailable_MinMax.dat",row.names=FALSE,col.names=TRUE,quote=FALSE)
```
In this step the unique locations of CMLs are obtained based on those intervals and sub-links with (zero) path-averaged rainfall intensities. This is needed in the next step in order to determine whether the coordinate of a grid cell is close to a start or end coordinate of a CML.

5. Start R and run this script to obtain a grid which can be used for plotting with MapRAINLINK: 
``` 
source("MakeInterpolationGridForPlottingSriLankaMinMax.R")
```
This produces the file "InterpolationGrid_SriLanka_Plus_Indices.dat", which contains the indices for longitude and latitude, i.e., the column and row numbers of an array. Final output is the file "CMLInterpolationGridSriLanka.dat", which is the same interpolation grid as provided at https://doi.org/10.4121/14166539.v2. If "GridPlottingCML" is set to "yes", the last column in "CMLInterpolationGridSriLanka.dat" is set to 1 only if a grid cell is less than 0.05 degrees "distance" from the start and/or end coordinates of one or more CML. Otherwise, the last column is always set to 1. "VisualizeCMLsGaugesRadars_CartopyOSM_GM.py" only plots the interpolated CML rainfall for values of 1 in the last column. In case "GridPlottingCML" is not set to "yes", step 4 can be skipped, and still an interpolation grid file is produced.

6. Simply follow the above steps for deriving your own interpolation grids for interpolation with RAINLINK and/or visualization with MapRAINLINK.

# Warning message
It seems that this warning message is not a problem:

"UserWarning: The input coordinates to pcolormesh are interpreted as cell centers, but are not monotonically increasing or decreasing. This may lead to incorrectly calculated cell edges, in which case, please supply explicit cell edges to pcolormesh."

# Open datasets to use with MapRAINLINK
- Gridded rainfall maps retrieved from CML data from Sri Lanka over a 3.5 month period (https://doi.org/10.4121/14166539.v2).
- The 2-day sample dataset, which is part of RAINLINK, or the ~4-month dataset covering the Netherlands (https://doi.org/10.4121/uuid:323587ea-82b7-4cff-b123-c660424345e5). These datasets need to be processed with RAINLINK to obtain CML path averages or rainfall fields.
- 5 minute precipitation accumulations from climatological gauge-adjusted radar dataset for The Netherlands (1 km) in KNMI HDF5 format (https://dataplatform.knmi.nl/dataset/rad-nl25-rac-mfbs-5min-2-0 or https://dataplatform.knmi.nl/dataset/rad-nl25-rac-mfbs-em-5min-2-0). Can be used as a reference for CML rainfall estimates.
- 1 hour precipitation accumulations from climatological gauge-adjusted radar dataset for The Netherlands (1 km) in KNMI HDF5 format (https://dataplatform.knmi.nl/dataset/rad-nl25-rac-mfbs-01h-2-0). Can be used as a reference for CML rainfall estimates.
- 24 hour precipitation accumulations from climatological gauge-adjusted radar dataset for The Netherlands (1 km) in KNMI HDF5 format (https://dataplatform.knmi.nl/dataset/rad-nl25-rac-mfbs-24h-2-0). Can be used as a reference for CML rainfall estimates. The file "RAD_NL25_RAC_MFBS_24H_201805300800_NL.h5" contains 24-h accumulations for one interval and can be used for testing.

# Visualizing Dutch KNMI radar data
Previous versions of RAINLINK could accumulate Dutch gridded radar data and visualize radar images. With the Python script "AccumulateRadarHDF5KNMIListCount.py", radar data in KNMI-HDF5 format can be accumulated. The files in a certain directory or a list of provided files can be accumulated. Below is an example of accumulating 3 radar images with 5-min precipitation to 15-min precipitation accumulations using Python version 3. The file "RAD_NL25_RAC_RT_202107141500.h5" is only used to construct an output file in case all input files are not valid.
```
python AccumulateRadarHDF5KNMIListCount.py "RAD_NL25_RAC_MFBS_15min_201109102045_NL.h5" "RAD_NL25_RAC_MFBS_5min_201109102035_NL.h5  RAD_NL25_RAC_MFBS_5min_201109102045_NL.h5 RAD_NL25_RAC_MFBS_5min_201109102040_NL.h5" 1 1 RAD_NL25_RAC_RT_202107141500.h5 files RAD_NL25_RAC_MFBS_15min 0
```
Next, a map of 15-min radar precipitation accumulations can be visualized. Note that this is also an example how to overrule the settings in the configuration script:
```
python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py "PlotDataField = 'no'" "PlotDataFieldRadarGrid = 'no'" "PlotKNMIRadar = 'yes'" "KNMIRadarInputFileName = 'RAD_NL25_RAC_MFBS_15min_201109102045_NL.h5'" "PlotCMLTimeInterval = 'no'" "TitlePlot = 'The Netherlands - Radar'" "OutputFileName = 'NetherlandsRadar.jpg'"
```
<img src="NetherlandsRadar.jpg" alt="drawing" width="500"/>

Map made with Natural Earth. Free vector and raster map data &copy naturalearthdata.com. Note that this reproduces Figure 5 (right) from https://doi.org/10.5194/amt-9-2425-2016 (but that was obtained with a different plotting package and different colors and classes). This map can be compared with the CML-based rainfall map given above for the same time interval. The color scale can be changed by modifying "ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py" or by supplying additional arguments ("levels = [0.1,1,2,3,4,5,6]" "LowestValue = 0.1"):
```
python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py ConfigVisualizeCMLsGaugesRadars_CartopyOSM_GM.py "PlotDataField = 'no'" "PlotDataFieldRadarGrid = 'no'" "PlotKNMIRadar = 'yes'" "KNMIRadarInputFileName = 'RAD_NL25_RAC_MFBS_15min_201109102045_NL.h5'" "PlotCMLTimeInterval = 'no'" "TitlePlot = 'The Netherlands - Radar'" "OutputFileName = 'NetherlandsRadarColorScale.jpg'" "levels = [0.1,1,2,3,4,5,6]" "LowestValue = 0.1" "DoColorSetUnder = 'yes'" "ScaleType = 'YellowRed'"
```
<img src="NetherlandsRadarColorScale.jpg" alt="drawing" width="500"/>

Note that when "DoColorSetUnder = 'yes'", a class below the lowest value of the supplied scale and a class above the highest value of the supplied scale are plotted in the legend and in the figure (but for HDF5 radar files in KNMI format, as used here, 0 precipitation values are never plotted, because these become nan values). Otherwise, these classes are not plotted in the legend, although the highest class is plotted in the figure for values exceeding the highest value of the supplied scale.

# Visualizing OPERA radar data
MapRAINLINK can also visualize gridded OPERA radar data in HDF5-ODIM format. Then the file "CoordinatesHDF5ODIMWGS84.dat" is needed. It contains the coordinates of the center of radar grid cells with longitude (first column) and latitude (second column) in degrees (WGS84). More tools for working with OPERA radar data, and the derived climatological dataset EURADCLIM, can be found here: https://github.com/overeem11/EURADCLIM-tools. EURADCLIM is a dataset of 1-h and 24-h precipitation accumulations covering 78% of geographical Europe (https://doi.org/10.21944/ymrk-mr24 & https://doi.org/10.21944/kcd7-3y59; version 2 of the dataset). The accompanying scientific article can be found here: https://doi.org/10.5194/essd-15-1441-2023 (version 1 of the dataset). The EURADCLIM file "RAD_OPERA_24H_RAINFALL_ACCUMULATION_201305311400.h5" contains 24-h accumulations for one interval and can be used for testing.

# Reference
When referring to MapRAINLINK, please check Zenodo and use the citation with doi to the current version (on https://doi.org/10.5281/zenodo.7611398 you can search for the newest version).

# MapRAINLINK
Python script to visualize rain gauge, radar, and commercial microwave link (CML) locations and their rainfall estimates on a map. This Python script is more widely applicable. It can be used to plot lines and/or points (and their values) on a map, not necessarily CML paths and rain gauge locations, and not necessarily rainfall intensities or accumulations. It can also be used to plot other data fields.

This script is meant as a replacement of the visualization by the open-source R package RAINLINK for Dutch radar data and commercial microwave link data. RAINLINK (https://github.com/overeem11/RAINLINK) is a retrieval algorithm for rainfall mapping from microwave links in a cellular communication network. It works with output data from RAINLINK, irrespective of the RAINLINK version. MapRAINLINK can be of general use for opportunistic sensing data, such as those from CML and personal weather stations (PWS).

How to start? Just clone this repository. With the default settings, the locations of CML paths and the interpolated CML rainfall estimates are plotted on a map of the Netherlands for one 15-min time interval. This is based on output from RAINLINK, by running the Dutch sample dataset, which is part of RAINLINK. The map is made by running:
python VisualizeCMLsGaugesRadars_CartopyOSM_GM.py

Open datasets to use with MapRAINLINK:
- Gridded rainfall maps retrieved from CML data from Sri Lanka over a 3.5 month period (https://doi.org/10.4121/14166539.v2).
- The 2-day sample dataset, which is part of RAINLINK, or the ~4-month dataset covering the Netherlands (https://doi.org/10.4121/uuid:323587ea-82b7-4cff-b123-c660424345e5). These datasets need to be processed with RAINLINK to obtain CML path averages or rainfall fields, which can be plotted with MapRAINLINK.
- 5 minute precipitation accumulations from climatological gauge-adjusted radar dataset for The Netherlands (1 km) in KNMI HDF5 format (https://dataplatform.knmi.nl/dataset/rad-nl25-rac-mfbs-5min-2-0). Can be used as a reference for CML rainfall estimates.
- 1 hour precipitation accumulations from climatological gauge-adjusted radar dataset for The Netherlands (1 km) in KNMI HDF5 format (https://dataplatform.knmi.nl/dataset/rad-nl25-rac-mfbs-01h-2-0). Can be used as a reference for CML rainfall estimates.
- 24 hour precipitation accumulations from climatological gauge-adjusted radar dataset for The Netherlands (1 km) in KNMI HDF5 format (https://dataplatform.knmi.nl/dataset/rad-nl25-rac-mfbs-24h-2-0). Can be used as a reference for CML rainfall estimates.

MapRAINLINK can also visualize gridded OPERA radar data in HDF5-ODIM format. More tools for working with OPERA radar data, and the derived climatological dataset EURADCLIM, can be found here: https://github.com/overeem11/EURADCLIM-tools. EURADCLIM is a dataset of 1-h and 24-h precipitation accumulations covering 78% of geographical Europa (https://doi.org/10.21944/7ypj-wn68 & https://doi.org/10.21944/1a54-gg96).

![](Netherlands.jpg)


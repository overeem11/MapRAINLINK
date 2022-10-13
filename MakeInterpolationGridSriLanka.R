# source('MakeInterpolationGridSriLanka.R')
# Make interpolation grid for Sri Lanka, which can be used by RAINLINK for the interpolation.
rm(list=ls(all=TRUE))

# Make interpolation grid for Sri Lanka, 0.02 degrees resolution:
Resolution <- 0.02
# Select geographical area for which to extract interpolation grid:
# http://www.openstreetmap.org/export can be used to select an appropriate area (bounding box). 
StartLon <- 79.66	# Start longitude in WGS84 (degrees).
StartLat <- 5.90	# Start latitude in WGS84 (degrees).
EndLon <- 81.90		# End longitude in WGS84 (degrees).
EndLat <- 9.86		# End latitude in WGS84 (degrees).
LatSeq <- seq(from=StartLat,to=EndLat,by=Resolution)	# Range of latitudes in degrees with step size.
LonSeq <- seq(from=StartLon,to=EndLon,by=Resolution)	# Range of longitudes in degrees with step size.
LatSeqData <- rep(LatSeq,length(LonSeq))
LonSeqData <- c(NA)
for (i in 1:length(LonSeq))
{
	print(i)
	LonSeqData <- c(LonSeqData,rep(LonSeq[i],length(LatSeq)))
}
LonSeqData <- LonSeqData[!is.na(LonSeqData)]
dataf <- data.frame(cbind(LonSeqData,LatSeqData))
write.table(dataf,"InterpolationGrid_SriLanka.dat",row.names=FALSE,col.names=c("X","Y"),quote=FALSE,sep=",")












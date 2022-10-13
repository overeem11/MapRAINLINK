# source('MakeInterpolationGridForPlottingSriLankaMinMax.R')
# Make grid for plotting CML rainfall maps for Sri Lanka for use with Python script "VisualizeCMLsGaugesRadars_CartopyOSM_GM.py".
rm(list=ls(all=TRUE))


# Make indices for longitude & latitude, so that Python script can make an array:
#################################################################################

# First read input file obtained by running R script "MakeInterpolationGridSriLanka.R":
u <- read.table("InterpolationGrid_SriLanka.dat",sep=",",skip=1)
LonUnique <- unique(u[,1])
LatUnique <- unique(u[,2])
IndexLon <- IndexLat <- c(NA)
for (i in 1:length(u[,1]))
{
    IndexLon[i] <- which(LonUnique==u[i,1]) - 1
}
for (i in 1:length(u[,2]))
{
    IndexLat[i] <- which(LatUnique==u[i,2]) - 1
}
dataf <- data.frame(cbind(u[,1],u[,2],IndexLon,IndexLat))
write.table(dataf,"InterpolationGrid_SriLanka_Plus_Indices.dat",sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)


# To create grid for plotting, only for locations within 0.05 degrees (~5 km) of a CML (use "yes"):
GridPlottingCML <- "yes"   # In case this is not equal to "yes", also a grid for plotting is made, but all grid points will be used.
if (GridPlottingCML=="yes")
{
    # Only take into account unique coordinates of CMLs for those CMLs that have data values:
    UniqueLinks <- read.table("CMLLocations_SriLanka_RmeanAvailable_MinMax.dat",skip=1)
    X <- c(UniqueLinks[,1],UniqueLinks[,3])
    Y <- c(UniqueLinks[,2],UniqueLinks[,4])


    # Determine which points of interpolation grid are within 0.05 degrees (~5 km) of a CML:
    LonGrid <- read.table("InterpolationGrid_SriLanka_Plus_Indices.dat",sep=",")[,1]
    LatGrid <- read.table("InterpolationGrid_SriLanka_Plus_Indices.dat",sep=",")[,2]
    MinDist <- c(NA)
    for (i in 1:length(LonGrid))
    {
	    print(i)
    	    DistLonDegrees <- DistLatDegrees <- c(NA)
	    for (z in 1:length(X) )
       	    {
		 DistLonDegrees[z] <- LonGrid[i] - X[z]
		 DistLatDegrees[z] <- LatGrid[i] - Y[z]
	    }
            Dist <- (DistLonDegrees^2 + DistLatDegrees^2)^0.5
	    MinDist[i] <- min(Dist,na.rm=T)
    }
    cond <- which(MinDist >= 0.05)
    # 1 = within 5 km at least 1 coordinate of a link.
    UseGrid <- rep(1,length(LonGrid))
    UseGrid[cond] <- 0
} else {
    LonGrid <- read.table("InterpolationGrid_SriLanka_Plus_Indices.dat",sep=",")[,1]
    UseGrid <- rep(1,length(LonGrid))
}
    
u <- read.table("InterpolationGrid_SriLanka_Plus_Indices.dat",sep=",")
dataf <- data.frame(cbind(u[,1],u[,2],u[,3],u[,4],UseGrid))
write.table("# Longitude[degrees],Latitude[degrees],ColumnNumber,RowNumber,CMLAntennaWithin5km[0=no;1=yes]","CMLInterpolationGridSriLanka.dat",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(dataf,"CMLInterpolationGridSriLanka.dat",sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
# This is the same interpolation grid as provided at https://doi.org/10.4121/14166539.v2





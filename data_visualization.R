library(rgdal) 
library(GISTools)
library(splancs)
library(gstat) # Use gstat's idw routine
library(sp) # Used for the spsample function
library(raster) # Used for clip out
library(tmap) # Install "XML" first, install.packages("XML", type="binary")
setwd("E:/GitHub/Tectonics-EOF")

# Ref: https://mgimond.github.io/Spatial/interpolation-in-r.html
# save.image("0824_work.RData")
load("0824_work.RData")

### Si Visualization

### Taiwan ShapeFile
TW <- readOGR(dsn = "./shp", layer = "Taiwan_county", encoding="utf8") #TWD97
TW <- spTransform(TW, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")) #WGS84
TW.lim <- c(TW@bbox[1,], TW@bbox[2,])


### PC max normalize
df_si_2['PC1_N'] <- df_si_2$PC1 / max(df_si_2$PC1)
df_si_2['PC2_N'] <- df_si_2$PC2 / max(df_si_2$PC2)

### Create SpatialPointsDataFrame with PC
P <- SpatialPointsDataFrame(cbind(df_si_2$lon, df_si_2$lat), df_si_2['PC2_N'], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))

# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(P, "regular", n=500000))
names(grd) <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd) <- TRUE  # Create SpatialPixel object
fullgrid(grd) <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(grd) <- proj4string(P)

# Interpolation using idw
P.idw <- gstat::idw(PC2_N ~ 1, P, newdata=grd, idp=2)

# Interpolation using kirging

# Add X and Y to P
P$X <- coordinates(P)[,1] 
P$Y <- coordinates(P)[,2]
# Define the 1st order polynomial equation
f.1 <- as.formula(PC2_N ~ X + Y) 

# Compute the sample variogram; note that the f.1 trend model is one of the
# parameters passed to variogram(). This tells the function to create the 
# variogram on the de-trended data.
var.smpl <- variogram(f.1, P, cloud = FALSE, cutoff=1000000, width=89900)

# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=14, model="Sph", range=590000, nugget=0))

# Perform the krige interpolation
dat.krg <- krige(f.1, P, grd, dat.fit)

# Convert to raster object then clip to Taiwan
# r <- raster(P.idw)
r <- raster(dat.krg)
r.m  <- mask(r, TW)

# Plot
tm_shape(r.m) + 
  tm_raster(n=10, palette="-RdYlBu", title="Value", breaks=c(-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0), midpoint=NA) + 
  tm_shape(P) + 
    tm_dots(size=0.01, palette="-RdYlBu") + 
  tm_shape(SpatialPoints(cbind(120.54, 22.92), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))) + # ¾_¥¡¦ì¸m
    tm_dots(size=1, col="red", shape=8) + 
  tm_legend(legend.outside=TRUE)

### Ti Visualization

library(ggplot2)

df_ti_2 <- as.data.frame(Ti)[1:2]
date <- seq(as.Date("2016-01-22"), as.Date("2016-02-21"), by="days")

ggplot(df_ti_2, aes(x=date, y=PC2)) + 
  geom_line() +
  theme_minimal()

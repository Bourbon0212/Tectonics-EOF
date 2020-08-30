library(rgdal) 
library(GISTools)
library(splancs)
library(ggplot2)
library(gstat) # Use gstat's idw routine
library(sp) # Used for the spsample function
library(raster) # Used for clip out
library(tmap) # Install "XML" first, install.packages("XML", type="binary")
setwd("E:/GitHub/Tectonics-EOF")

visualize_si <- function(df_si, center, pc) {
  ### Taiwan ShapeFile
  TW <- readOGR(dsn = "./shp", layer = "Taiwan_county", encoding="utf8") #TWD97
  TW <- spTransform(TW, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")) #WGS84
  
  ### Maximum Normalization
  df_si['PC1N'] <- df_si$PC1 / max(df_si$PC1)
  df_si['PC2N'] <- df_si$PC2 / max(df_si$PC2)
  
  ### Process Spatial Data
  P <- SpatialPointsDataFrame(cbind(df_si_2$lon, df_si_2$lat), df_si_2['PC2N'], 
                              proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))
  grd <- as.data.frame(spsample(P, "regular", n=500000)) # Create an empty grid where n is the total number of cells
  names(grd) <- c("X", "Y")
  coordinates(grd) <- c("X", "Y")
  gridded(grd) <- TRUE  # Create SpatialPixel object
  fullgrid(grd) <- TRUE  # Create SpatialGrid object
  proj4string(grd) <- proj4string(P) # Add P's projection to the empty grid
  
  ### Interpolation using Kriging
  P$X <- coordinates(P)[,1] # Add X and Y to P
  P$Y <- coordinates(P)[,2]
  f.1 <- NA # Define the 1st order polynomial equation
  if (pc == "PC1") {
    f.1 <- as.formula(PC1N ~ X + Y) 
  } else if (pc == "PC2") {
    f.1 <- as.formula(PC2N ~ X + Y)
  }
  # Compute the sample variogram
  # Note that the f.1 trend model is one of the parameters passed to variogram()
  # This tells the function to create the variogram on the de-trended data
  var.smpl <- variogram(f.1, P, cloud = FALSE, cutoff=1000000, width=89900)
  # Compute the variogram model by passing the nugget, sill and range values to fit.variogram() via the vgm() function.
  dat.fit <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                            vgm(psill=14, model="Sph", range=590000, nugget=0))
  dat.krg <- krige(f.1, P, grd, dat.fit) # Perform the krige interpolation
  
  ### Interpolation using idw
  # P.idw <- gstat::idw(PC2_N ~ 1, P, newdata=grd, idp=2)
  
  ### Convert to raster and clip to Taiwan
  r <- raster(dat.krg)
  # r <- raster(P.idw)
  r.m  <- mask(r, TW)
  
  ### Plot
  tm_shape(r.m) + 
    tm_raster(n=10, palette="-RdYlBu", title="Value", breaks=c(-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0), midpoint=NA) + 
    tm_shape(P) + 
    tm_dots(size=0.01, palette="-RdYlBu") + 
    tm_shape(SpatialPoints(cbind(center[1], center[2]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))) + # ¾_¥¡¦ì¸m
    tm_dots(size=1, col="red", shape=8) + 
    tm_legend(legend.outside=TRUE)
}

visualize_ti <- function(df_ti, pc) {
  ggplot(df_ti, aes_string(x="date", y=pc)) + 
    geom_line() +
    theme_minimal()
}
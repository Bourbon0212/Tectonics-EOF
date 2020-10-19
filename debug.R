library(readr)
library(tidyverse)
library(zoo)
library(pracma)
library(factoextra)
library(rgdal) 
library(GISTools)
library(splancs)
library(ggplot2)
library(gstat) # Use gstat's idw routine
library(sp) # Used for the spsample function
library(raster) # Used for clip out
library(tmap) # Install "XML" first, install.packages("XML", type="binary")

file_path <- list.files("data_raw") # COR files of all stations
df_sta <- readRDS("df_sta.rds")

### Get the stations of studying area
lat = 24

sta_study <- c()
for (i in 1:length(file_path)) {
  tmp <- read.table(paste0("data_raw/", file_path[i]), col.names = c("time", "lat", "lon", "hgt", "N", "E", "U", "X"))
  
  if (tmp$lat[1] <= lat) { # Based on latitude
    sta_study <- c(sta_study, file_path[i]) 
  } else {
    # print(paste(file_path[i], "is not in the study area."))
  }
}

### Loop loading files
component = "U"

start_path <- paste0("data_raw/", sta_study[1])
df <- read.table(start_path, col.names = c("time", "lat", "lon", "hgt", "N", "E", "U", "X"))[c("time", component)]
for (i in 2:length(sta_study)) {
  # print(sta_study[i])
  tmp <- read.table(paste0("data_raw/", sta_study[i]), col.names = c("time", "lat", "lon", "hgt", "N", "E", "U", "X"))[c("time", component)]
  df <- full_join(df, tmp, by = "time")
}
colnames(df) <- c("time", sta_study) # Rename colname to sta_study
df <- arrange(df, time)

### Time filtering & interpolation

time <- "2010-03-04"
pre <- 20
post <- 50

start <- as.Date(time) - pre; end <- as.Date(time) + post
diff <- as.numeric(end - start)

start_year <- as.numeric(substring(start, 1, 4))
start_julian <- as.numeric(format(as.Date(start), '%j'))
start_time <- round(start_year + ((start_julian - 0.5) / 366), 5)
end_year <- as.numeric(substring(end, 1, 4))
end_julian <- as.numeric(format(as.Date(end), '%j'))
end_time <- round(end_year + ((end_julian - 0.5) / 366), 5)

df_filter <- filter(df, time >= start_time & time <= end_time)
df_select <- df_filter[colSums(!is.na(df_filter)) >= (diff * 5/6)] # Drop stations with too many NA
Cz <- zoo(df_select) # Interpolate with zoo
df_filled <- as.data.frame(na.fill(na.approx(Cz), "extend"))

### Data Detrend
colnames <- df_filled$time
df_filled_t <- t(df_filled[2:length(df_filled)])
colnames(df_filled) <- colnames
df_detrend_t <- detrend(df_filled_t, tt = 'constant')
df_detrend <- t(df_detrend_t)
rownames(df_detrend) <- colnames

### PCA Process

n <- dim(df_detrend)[1] # number of temporal samples
p <- dim(df_detrend)[2] # number of points with data

if (n >= p) {
  # If n>=p we estimate the EOF in the State Space Setting, i.e., using Z'*Z which is p x p.
  R <- t(df_detrend) %*% df_detrend
  data.pca <- prcomp(R)
  pca_eigenvector <- data.pca$rotation
  pca_eigenvalue <- data.pca$sdev ** 2
  
  Si <- as.data.frame(pca_eigenvector) # 丁家Α
  Ti <- as.data.frame(df_detrend %*% pca_eigenvector) # 丁家Α
  # We estimate the expansion coefficients (time series) associated to each EOF
  # for (i in 1:nrow(Ti)) {
    # Ti[,i] = df_detrend %*% pca_eigenvector[,i]
  # }
  
} else if (n < p) {
  # If n<p  we estimate the EOF in the Sample Space Setting, i.e., using Z*Z' which is n x n.
  R <- df_detrend %*% t(df_detrend)
  data.pca <- prcomp(R)
  pca_eigenvector <- data.pca$rotation
  pca_eigenvalue <- data.pca$sdev ** 2
  
  Si <- as.data.frame(pca_eigenvector) # 丁家Α
  Ti <- as.data.frame(t(df_detrend) %*% pca_eigenvector) # 丁家Α
  # Ti1 <- matrix(nrow=p, ncol=n)
  # We estimate the expansion coefficients (time series) associated to each EOF
  # for (i in 1:length(Ti)) {
    # Ti1[,i] = t(df_detrend) %*% pca_eigenvector[,i]
  # }
  # Si = Ti %/% sqrt(pca_eigenvalue)
  # Ti = sqrt(pca_eigenvalue) %*% pca_eigenvector
  
  # Expansion coefficients
  # for (i in 1:nrow(Si)) {
    # Si[,i] = Ti[,i] %/% sqrt(pca_eigenvalue[i])
  # }
  # for (i in 1:nrow(Ti)) {
    # Ti[,i] = sqrt(pca_eigenvalue[i]) %*% pca_eigenvector[,i]
  # }
}


###  PCA
# df_detrend <- df_filled[2:length(df_filled)]
# data.pca <- prcomp(df_detrend)
# fviz_eig(data.pca)
# pca_eigenvector <- data.pca$rotation; D <- data.matrix(df_detrend)
# pca_eigenvalue <- data.pca$sdev ** 2
# Si <- as.data.frame(data.pca$rotation) # 丁家Α
# Ti <- as.data.frame(D %*% pca_eigenvector) # 丁家Α

### Normalization
for (i in 1:length(Si)) {
  Si[i] <- Si[i] / sd(Si[[i]])
}
for (i in 1:length(Ti)) {
  Ti[i] <- Ti[i] * sd(Si[[i]])
}

Ti$date <- seq(as.Date(start), as.Date(end), by="days")

### Binding Si with Locations
Si <- Si[,1:2]
Si <- cbind(newColName = substring(rownames(Si),1,4), Si)
rownames(Si) <- 1:nrow(Si)
colnames(Si) <- c("name", "PC1", "PC2")
df_si <- inner_join(Si, df_sta, by="name")

### Taiwan ShapeFile
TW <- readOGR(dsn = "./shp", layer = "Taiwan_county", encoding="utf8") #TWD97
TW <- spTransform(TW, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")) #WGS84

### Process Spatial Data
P <- SpatialPointsDataFrame(cbind(df_si$lon, df_si$lat), df_si[c('PC1N', 'PC2N')], 
                            proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))
grd <- as.data.frame(spsample(P, "regular", n=50000)) # Create an empty grid where n is the total number of cells
names(grd) <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd) <- TRUE  # Create SpatialPixel object
fullgrid(grd) <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(P) # Add P's projection to the empty grid

### Interpolation using Kriging
P$X <- coordinates(P)[,1] # Add X and Y to P
P$Y <- coordinates(P)[,2]
P <- sp::remove.duplicates(P) # Remove duplicate points
f.1 <- as.formula(PC1N ~ X + Y) 
# Compute the sample variogram
# Note that the f.1 trend model is one of the parameters passed to variogram()
# This tells the function to create the variogram on the de-trended data
var.smpl <- variogram(f.1, P, cloud = FALSE, cutoff=1000000, width=89900)
# Compute the variogram model by passing the nugget, sill and range values to fit.variogram() via the vgm() function.
dat.fit <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                         vgm(psill=14, model="Sph", range=590000, nugget=0))
dat.krg <- krige(f.1, P, grd, dat.fit) # Perform the krige interpolation
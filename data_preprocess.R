library(readr)
library(tidyverse)
library(zoo)
setwd("E:/GitHub/Tectonics-EOF")

data_process <- function(time, pre, post, lat, component) {
  file_path <- list.files("data_raw") # COR files of all stations
  
  ### Get the stations of studying area
  sta_study <- c()
  for (i in 1:length(file_path)) {
    tmp <- read.table(paste0("data_raw/", file_path[i]), col.names = c("time", "lat", "lon", "hgt", "E", "N", "U", "X"))
    
    if (tmp$lat[1] <= lat) { # Based on latitude
      sta_study <- c(sta_study, file_path[i]) 
    } else {
      # print(paste(file_path[i], "is not in the study area."))
    }
  }
  
  ### Loop loading files
  start_path <- paste0("data_raw/", sta_study[1])
  df <- read.table(start_path, col.names = c("time", "lat", "lon", "hgt", "E", "N", "U", "X"))[c("time", component)]
  for (i in 2:length(sta_study)) {
    # print(sta_study[i])
    tmp <- read.table(paste0("data_raw/", sta_study[i]), col.names = c("time", "lat", "lon", "hgt", "E", "N", "U", "X"))[c("time", component)]
    df <- full_join(df, tmp, by = "time")
  }
  colnames(df) <- c("time", sta_study) # Rename colname to sta_study
  df <- arrange(df, time)
  
  ### Time filtering & interpolation
  
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
  
  ###  PCA
  data.pca <- prcomp(df_filled[2:length(df_filled)])
  Si <- as.data.frame(data.pca$rotation) # 空間模式
  pca_eigenvector <- data.pca$rotation; D = data.matrix(df_filled[,2:length(df_filled)])
  Ti <- as.data.frame(D %*% pca_eigenvector) # 時間模式
  Ti$date <- seq(as.Date(start), as.Date(end), by="days")
  
  return(list(Si = Si, Ti = Ti))
}

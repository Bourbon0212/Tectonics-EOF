setwd("E:/GitHub/Tectonics-EOF")

### Stations List with Locations
file_path <- list.files("data_raw") # COR files of all stations

sta_name <- c()
sta_lat <- c()
sta_lon <- c()

for (i in 1:length(file_path)) {
  tmp <- read.table(paste0("data_raw/", file_path[i]), col.names = c("time", "lat", "lon", "hgt", "E", "N", "U", "X"))
  print(file_path[i])
  sta_name[i] <- substring(file_path[i],1,4)
  sta_lat[i] <- round(tmp$lat[1], 5)
  sta_lon[i] <- round(tmp$lon[1], 5)
}

df_sta <- data.frame(name = sta_name, lat = sta_lat, lon = sta_lon)

### Tidying, Interpolations & EOF
library(readr);library(tidyverse);library(zoo)
if(!exists("data_process", mode="function")) source("data_preprocess.R")

start <- "2016-01-22"
end <- "2016-02-21"

cc <- data_process(start, end, 24, "E")

### Binding Si with Locations
Si <- cc$Si[,1:2]
Si <- cbind(newColName = substring(rownames(Si),1,4), Si)
rownames(Si) <- 1:nrow(Si)
colnames(Si) <- c("name", "PC1", "PC2")
df_si <- inner_join(Si, df_sta, by="name")
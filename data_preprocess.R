library(readr)
library(tidyverse)
setwd("E:/GitHub/Tectonics-EOF")

### 20160206 美濃地震，Julian Date 時間格式 (天數 - 0.5)/366

file_path <- list.files("data_raw")

### Get the STAs of studying area

sta_study <- c()

for (i in 1:length(file_path)) {
  tmp <- read.table(paste0("data_raw/", file_path[i]), col.names = c("time", "lat", "lon", "hgt", "E", "N", "U", "X"))
  
  if (tmp$lat[1] <= 24) { # Based on latitude
    sta_study <- c(sta_study, file_path[i]) 
  } else {
    print(paste(file_path[i], "is not in the study area."))
  }
}

### Loop loading files

start_path <- paste0("data_raw/", sta_study[1])
df <- read.table(start_path, col.names = c("time", "lat", "lon", "hgt", "E", "N", "U", "X"))[c("time","E")]

for (i in 2:length(sta_study)) {
  print(sta_study[i])
  tmp <- read.table(paste0("data_raw/", sta_study[i]), col.names = c("time", "lat", "lon", "hgt", "E", "N", "U", "X"))[c("time","E")]
  df <- full_join(df, tmp, by = "time")
}

colnames(df) <- c("time", sta_study) # Rename colname to sta_study
df <- arrange(df, time)

### Filter the data between 20160122 - 20160221

start <- round(2016 + ((22 - 0.5) / 366), 5)
end <- round(2016 + ((52 - 0.5) / 366), 5)

df_filter <- filter(df, time >= start & time <= end)
df_select <- df_filter[colSums(!is.na(df_filter)) >= 25] # Drop stations with too many NA

### Deal with Missing Values

Cz <- zoo(df_select) # Interpolate with zoo
df_filled <- as.data.frame(na.fill(na.approx(Cz), "extend"))

###  PCA
data.pca <- prcomp(df_filled[2:237])
pca_eigenvector <- as.data.frame(data.pca$rotation[,1:2]) # Si (空間模式)
pca_eigenvalue <- as.data.frame(data.pca$x[,1:2])

Si = data.pca$rotation; D = data.matrix(df_filled[,2:237])
Ti = D %*% Si # Si (時間模式)

# save.image("0820_work.RData")
load("0820_work.RData")

library(readr)
library(tidyverse)
library(zoo)
library(pracma)
library(factoextra)
library(ggplot2)
setwd("E:/GitHub/Tectonics-EOF")

time <- "2010-03-04" # 2010.171
pre <- 33
post <- 59

# "PAOL" (23.10863, 120.7029)
# "TAYN" (23.15935, 120.7642)

sta_study <- c('AKND','ALIS','BANP','BDES','BKBL','C002','CAOT','CHEN','CHGO','CHIA',
               'CKSV','CLON','CTOU','CWEN','DASI','DAWU','DCHU','DNAN','DNFU','DSIN',
               'ERLN','ERPN','FALI','FENP','FKDO','FLNM','FUGN','FUNY','GAIS','GFES',
               'GGDS','GS05','GS06','GS07','GS17','GS18','GS25','GS26','GS27','GS28',
               'GS29','GS30','GS31','GS33','GS34','GS35','GS41','GS42','GS43','GS44',
               'GS45','GS46','GUFU','GUKN','HENC','HLIU','HNES','HRGN','HUAL','HUWE',
               'ICHU','JHCI','JLUT','JPEI','JPIN','JSUI','JULI','KAFN','KASH','KASU',
               'KAWN','KDNM','KNKO','KUAN','KULN','KZN1','LANY','LAOL','LGUE','LIKN',
               'LIUC','LNCH','LNJS','LONT','LUDA','MAJA','MITO','MLON','MOTN','NCKU',
               'NDHU','NEMN','NHSI','NJOU','PAKU','PAOL','PEIN','PING','PINT',        # 'PEIM' 
               'PKGM','PLIM','PTUN','PUSN','RENI','RFES','S011','S103','S104','S105',
               'S106','S167','S170','SAND','SANL','SCES','SGAN','SHAN','SHUL','SHWA',
               'SINL','SINY','SOFN','SSUN','SUAN','SUN1','T103','TAPE','TAPO','TAPU',
               'TASI','TATA','TAYN','TKJS','TMAL','TMAM','TSLN','TTSH','TUNS','VR02', # 'VS03'
               'W021','WANS','WARO','WDAN','WULU','WUST','YENL','YSAN','YSLL',
               'YULI','YUSN','ZEND'                                                   # 'ZWEN'
)

### Loop loading files
start_path <- paste0("data_raw/", sta_study[1], ".COR")
df <- read.table(start_path, col.names = c("time", "lat", "lon", "hgt", "N", "E", "U", "X"))[c("time", "U")]
for (i in 2:length(sta_study)) {
  # print(sta_study[i])
  tmp <- read.table(paste0("data_raw/", sta_study[i], ".COR"), col.names = c("time", "lat", "lon", "hgt", "N", "E", "U", "X"))[c("time", "U")]
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

### Data Detrend
colnames <- df_filled$time
df_filled_t <- t(df_filled[2:length(df_filled)])
df_detrend_t <- detrend(df_filled_t, tt = 'linear')
df_detrend <- as.data.frame(t(df_detrend_t))
rownames(df_detrend) <- colnames
df_detrend <- cbind(time = as.numeric(rownames(df_detrend)), df_detrend)
rownames(df_detrend) <- NULL

###  PCA
# df_detrend <- df_filled[2:length(df_filled)]
data.pca <- prcomp(df_detrend)

### Plot time series
ggplot() +
  geom_line(data = df_filled, aes(x=time, y=PAOL), col='red') +
  geom_line(data = df_detrend, aes(x=time, y=PAOL), col='blue') +
  geom_vline(xintercept=2010.171, linetype="dotted", col='green')

filter(df, time >= 2008) %>%
ggplot() + 
  geom_line(aes(x=time, y=PAOL)) +
  geom_vline(xintercept=2010.171, linetype="dotted", col='green')

ggplot() +
  geom_line(data = df_filled, aes(x=time, y=TAYN), col='red') +
  geom_line(data = df_detrend, aes(x=time, y=TAYN), col='blue') +
  geom_vline(xintercept=2010.171, linetype="dotted", col='green')

filter(df, time >= 2008) %>%
  ggplot() + 
  geom_line(aes(x=time, y=TAYN)) +
  geom_vline(xintercept=2010.171, linetype="dotted", col='green')
  

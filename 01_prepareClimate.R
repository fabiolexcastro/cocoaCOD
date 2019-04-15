
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, gtools, tidyverse, stringr, velox, sf, foreach, doSNOW, parallel)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Functions to use --------------------------------------------------------
extMask <- function(vr){
  # vr <- 'prec'
  fl <- grep(vr, fls, value = TRUE)
  st <- stack(fl)
  ct <- raster::crop(st, adm_dss)
  ct <- raster::mask(ct, adm_dss)
  ct <- unstack(ct)
  Map('writeRaster', x = ct, filename = paste0('../data/raster/clm/crn/', basename(fl), '.asc'), overwrite = TRUE)
  return(ct)
}
extMaskFtr <- function(gc, yr){
  # gc <- m50[1]
  # yr <- yrs[[2]]
  print(paste0(yr, ' ', gc))  
  fl <- grep(yr, fls, value = TRUE) %>%
    grep(gc, ., value = TRUE)
  st <- stack(fl) 
  ct <- raster::crop(st, adm_dss)
  ct <- raster::mask(ct, adm_dss)
  ct <- unstack(ct)
  dir.create(paste0('../data/raster/clm/ftr/', yr, '/', gc))
  pt <- paste0('../data/raster/clm/ftr/', yr, '/', gc)
  Map('writeRaster', x = ct, filename = paste0(pt, '/', basename(fl), '.asc'), overwrite = TRUE)
  print('Done!')
}
reviewGCM <- function(yr){
  # yr <- '2020_2049'
  x <- paste0('../data/raster/clm/ftr/', yr) %>%
    list.files()
  y <- gcm
  z <- setdiff(y, x)
  return(z)
}

# Load data ---------------------------------------------------------------
lbr <- raster::getData(name = 'GADM', country = 'LBR', level = 0)
cod <- raster::getData(name = 'GADM', country = 'LBR', level = 0)
srr <- raster::getData(name = 'GADM', country = 'SLE', level = 0)

adm <- rbind(lbr, cod, srr)
writeOGR(obj = adm, dsn = '../data/shp/base', layer = 'limite_countries', driver = 'ESRI Shapefile')
rm(lbr, cod, srr)
adm@data$gid <- 1
adm_dss <- aggregate(adm, 'gid')

# Extract by mask - Current -----------------------------------------------
pth <- '//dapadfs/data_cluster_4/observed/gridded_products/worldclim/Global_30s'
vrs <- c(paste0('prec_', 1:12, '$'), paste0('tmin_', 1:12, '$'), paste0('tmean_', 1:12, '$'), paste0('tmax_', 1:12, '$'))
fls <- list.files(pth, full.names = TRUE) %>%
  grep(paste0(vrs, collapse = '|'), ., value = TRUE) %>%
  mixedsort()
vrs <- c('prec', 'tmin', 'tmean', 'tmax')
crn <- Map('extMask', vr = vrs)

# Extract by mask - Future ------------------------------------------------
pth <- '//dapadfs/data_cluster_2/gcm/cmip5/downscaled/rcp60/global_30s'
gcm <- list.files(pth)
yrs <- c('2020_2049', '2040_2069')
fls <- paste0(pth, '/', gcm)
f30 <- paste0(fls, '/', 'r1i1p1', '/', yrs[[1]]) %>%
  list.files(full.names = TRUE) %>%
  grep(paste0(vrs, collapse = '|'), ., value = TRUE) %>%
  mixedsort()
f50 <- paste0(fls, '/', 'r1i1p1', '/', yrs[[2]]) %>%
  list.files(full.names = TRUE) %>%
  grep(paste0(vrs, collapse = '|'), ., value = TRUE) %>%
  mixedsort()
fls <- c(f30, f50)

cl <- makeCluster(19)
registerDoSNOW(cl)

foreach(i = 1:length(gcm), .packages = c('raster', 'rgdal', 'gtools', 'foreach', 'sp', 'stringr', 'tidyverse'), .verbose = TRUE) %dopar% {
  foreach(y = 1:length(yrs)) %do% {
    print(gcm[i])
    extMaskFtr(gc = gcm[i], yr = yrs[y])  
  }
}

# Review GCM's ------------------------------------------------------------
miss <- Map('reviewGCM', yr = yrs)



# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse, gtools, dismo, foreach, parallel, doSNOW)

# Initial setup  ----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)
source('functionsR.r')

# Load data ---------------------------------------------------------------
msk <- shapefile('../data/shp/base/limite_countries.shp')
yrs <- list.files('../data/raster/clm/ftr')
gcm <- list.files('../data/raster/clm/ftr/2020_2049')
xtr <- list.files('../data/raster/clm/crn', full.names = TRUE, pattern = 'et_sol') %>% 
  mixedsort() %>%
  stack()
msk_lyr <- raster('../data/raster/clm/crn/bio1.asc') * 0 + 1

# Calculating the etp variables by each gcm -------------------------------
cl <- makeCluster(18)
registerDoSNOW(cl)
foreach(i = 1:length(gcm), .packages = c('raster', 'rgdal', 'dplyr', 'gtools', 'foreach', 'sp', 'stringr', 'dismo'), .verbose = TRUE) %dopar% {
  calcETP(gc = gcm[i], yr = yrs[2])
}
stopCluster(cl)



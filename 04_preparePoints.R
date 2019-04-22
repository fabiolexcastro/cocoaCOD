

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, gtools, tidyverse, stringr, velox, sf, foreach, doSNOW, parallel, rgbif)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
gbl <- read_csv('../data/tbl/cocoa_global.csv')
pnt <- shapefile('../data/tbl/_codCocoa/IturiCocoaFields_geo.shp')
adm <- shapefile('../data/shp/base/limite_countries.shp')
msk <- raster('../data/raster/clm/crn/prec_1.asc') * 0 + 1 

# Download cocoa data from GBIF -------------------------------------------
occ <- occ_data(scientificName = 'Theobroma cacao',
                limit = 20000,
                hasCoordinate = TRUE,
                hasGeospatialIssue = FALSE)[[2]]
occ <- occ %>% 
  filter(country == 'Congo, Democratic Republic of the') %>%
  dplyr::select(name, decimalLongitude, decimalLatitude)

# Filter data -------------------------------------------------------------
lb1 <- gbl %>% filter(Country %in% c('Liberia', 'Democratic Republic of the Congo', 'Sierra Leone'))

# Join the databases ------------------------------------------------------
lb1 <- lb1 %>% 
  mutate(id = 1:nrow(.)) %>% 
  dplyr::select(id, Longitude, Latitude) %>% 
  setNames(c('id', 'lon', 'lat'))
occ <- occ %>%
  dplyr::select(-name) %>%
  mutate(id = 1:nrow(.)) %>%
  dplyr::select(id, decimalLongitude, decimalLatitude) %>%
  setNames(c('id', 'lon', 'lat'))
occ <- rbind(lb1, occ)
pnt <- coordinates(pnt) %>%
  as_tibble() %>%
  mutate(id = 1:nrow(.)) %>% 
  setNames(c('lon', 'lat', 'id')) %>%
  dplyr::select(id, lon, lat)
pnt <- rbind(occ, pnt)
saveRDS(object = pnt, file = '../data/rds/pnt.rds')
coordinates(pnt)

rm(occ, lb1, gbl)

# Extract climate values --------------------------------------------------
stk <- list.files('../data/raster/clm/crn', full.names = TRUE, pattern = '.asc$') %>%
  mixedsort() %>%
  stack()
crn_vls <- raster::extract(stk, pnt[,2:3]) %>%
  cbind(., pnt[,2:3]) %>% 
  as_tibble()
crn_vls <- crn_vls[complete.cases(crn_vls),]














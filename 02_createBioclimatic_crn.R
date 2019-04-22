
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse, gtools, dismo)

# Initial setup -----------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Functions to use --------------------------------------------------------
myExtract <- function(fle){
  print(fle)
  lyr <- raster(fle)
  cut <- raster::crop(lyr, msk)
  cut <- raster::mask(cut, msk)
  print('To write')
  writeRaster(x = cut, filename = paste0('../_data/_tif/_climate/_crn/', gsub('.tif', '.asc', basename(fle))))
}
cumTemp <- function(x) {
  
  p <- matrix(nrow = 1, ncol = 4)
  colnames(p) <- paste('bio', 21:24, sep = '')
  
  w <- x[25:36] ### tmax
  y <- x[13:24] ### tmean
  x <- x[1:12]  ### Prec-PET
  z <- x
  
  ### if the values are NA the bios are NA
  if(all(is.na(x))) {
    p[,'bio21'] <- NA
    p[,'bio22'] <- NA
    p[,'bio23'] <- NA
    p[,'bio24'] <- NA
  } else {
    
    ## cumulative deficit to determine dry season (=Bio22)
    
    # print('Bio 22...')
    
    x <- z
    lng <- length(x)
    x <- c(x, x[1:12])
    x[x>0] <- NA
    cumdef <- matrix(ncol = 12, nrow = lng)
    for (i in 1:12) {
      cumdef[, i] <- x[i:(lng + i - 1)]
    }
    p[,'bio22'] <- min(c(0,apply(cumdef, MARGIN = 1, FUN = cumsum)),na.rm=T)
    
    ## cumulative surplus to determine growing season
    x <- z
    lng <- length(x)
    x <- c(z, z[1:12])
    x[x<0] <- NA
    cumplus <- matrix(ncol = 12, nrow = lng)
    
    for (i in 1:12) {
      
      cumplus[, i] <- x[i:(lng + i - 1)]
      
    }
    
    ### If there is no dry season
    ### the length becomes 0
    ### the growing season temp is the mean of monthly mean temp
    ### the dry season max temp is the max temp of the driest month 
    
    if(p[,'bio22']==0){
      
      p[,'bio21'] <- 0
      p[,'bio23'] <- mean(y)
      p[,'bio24'] <- w[which.min(z)]
      
    } else {
      
      ### the mean temperatures for all possible seasons
      y <- c(y, y[1:12])
      n <- matrix(ncol = 12, nrow = lng)
      for (i in 1:12) {
        
        n[, i] <- y[i:(lng + i - 1)]
        
      }
      
      meantemp <- apply(n, MARGIN = 1, FUN = cumsum)
      
      ### the max temperatures for all possible seasons
      w <- c(w, w[1:12])
      n <- matrix(ncol = 12, nrow = lng)
      
      for (i in 1:12) {
        
        n[, i] <- w[i:(lng + i - 1)]
        
      }
      maxtemp <- apply(n, MARGIN = 1, FUN = cumsum)
      
      ### Consecutive months with Prec<PET (=bio21)
      x <- z
      x <- c(x, x[1:12])
      x[x>0] <- NA
      x[x<0] <- 1
      o <- matrix(ncol = 12, nrow = lng)
      
      for (i in 1:12) {
        
        o[, i] <- x[i:(lng + i - 1)]
        
      }
      
      con_months <- max(apply(o,1,cumsum),na.rm=T)
      p[,'bio21'] <- con_months
      
      ### if the dry season is 12 months the growing season mean is the mean of the wettest month
      
      if(con_months==12){
        
        p[,'bio23'] <- y[which.max(z)]
        
      } else { 
        
        ### The meantemp of the wettest season
        p[,'bio23'] <- meantemp[which.max(apply(cumplus, MARGIN = 1, FUN = cumsum))]/(12-con_months)
        
      }
      ### The mean maxtemp of the driest season
      
      p[,'bio24'] <- maxtemp[which.min(apply(cumdef, MARGIN = 1, FUN = cumsum))]/con_months    
      
    }
    
  }
  
  return(p)
  
}
etpvars <- function(x){
  p <- matrix(nrow = 1, ncol = 9)
  colnames(p) = paste("bio", 25:33, sep = "")
  
  tavg <- x[25:36] ### Temp
  prec <- x[13:24] ### PREC
  pet <- x[1:12]  ### PET
  
  ### if the values are NA the bios are NA
  if(all(is.na(x))) {
    return(p)
  } else {
    
    window <- function(x)  {
      lng <- length(x)
      x <- c(x,  x[1:3])
      m <- matrix(ncol=3, nrow=lng)
      for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
      apply(m, MARGIN=1, FUN=sum)
    }
    
    ### BIO_25: Annual PET
    p[,1] <- sum(pet)
    ### BIO_26: PET seasonality (Coefficient of Variation)
    p[,2] <- cv(pet)
    ### BIO_27: MAX PET
    p[,3] <- max(pet)
    ### BIO_28: Min PET
    p[,4] <- min(pet)
    ### BIO_29: Range of PET (PETmax-PETmin)
    p[,5] <- p[,3]-p[,4]
    
    wet <- window(prec)
    hot <- window(tavg)/3
    pet2 <- c(pet,pet[1:2])
    
    ### BIO_30: PET of wettest quarter
    p[,6] <- sum(pet2[c(which.max(wet):(which.max(wet)+2))])
    ### BIO_31:	PET of driest quarter
    p[,7] <- sum(pet2[c(which.min(wet):(which.min(wet)+2))])
    ### BIO_32:	PET of warmest quarter
    p[,8] <- sum(pet2[c(which.max(hot):(which.max(hot)+2))])
    ### BIO_33:	PET of coldest quarter
    p[,9] <- sum(pet2[c(which.min(hot):(which.min(hot)+2))])
  }
  round(p,digits=2)
  return(p)
}
cumDry <- function(x) {
  p <- matrix(nrow = 1, ncol = 3)
  colnames(p) = paste("bio", c("20_40mm","20_50mm","20_100mm"), sep = "")
  
  if(all(is.na(x))) {
    p[,"bio20_40mm"] <- NA
    p[,"bio20_50mm"] <- NA
    p[,"bio20_100mm"] <- NA
  } else {
    ##40mm
    y <- x
    lng <- length(y)
    y <- c(y, y[1:12])
    y[y<40] <- 1
    y[y>=40] <- NA
    m <- matrix(ncol = 12, nrow = lng)
    for (i in 1:12) {
      m[, i] <- y[i:(lng + i - 1)]
    }
    cumdry <-  max(c(0,apply(m, MARGIN = 1, FUN = cumsum)),na.rm=T)
    p[,"bio20_40mm"] <- cumdry
    
    ##50mm
    y <- x
    lng <- length(y)
    y <- c(y, y[1:12])
    y[y<50] <- 1
    y[y>=50] <- NA
    m <- matrix(ncol = 12, nrow = lng)
    for (i in 1:12) {
      m[, i] <- y[i:(lng + i - 1)]
    }
    cumdry <-  max(c(0,apply(m, MARGIN = 1, FUN = cumsum)),na.rm=T)
    p[,"bio20_50mm"] <- cumdry
    
    ##100mm
    y <- x
    lng <- length(y)
    y <- c(y, y[1:12])
    y[y<100] <- 1
    y[y>=100] <- NA
    m <- matrix(ncol = 12, nrow = lng)
    for (i in 1:12) {
      m[, i] <- y[i:(lng + i - 1)]
    }
    cumdry <-  max(c(0,apply(m, MARGIN = 1, FUN = cumsum)),na.rm=T)
    p[,"bio20_100mm"] <- cumdry
  }
  return(p)
}

# Load data ---------------------------------------------------------------
msk <- shapefile('../data/shp/base/limite_countries.shp')

# Climate data
fls <- list.files('../data/raster/clm/crn', full.names = TRUE, patter = '.asc') %>%
  mixedsort()
var <- c('prec', 'tmin', 'tmean', 'tmax')
for(i in 1:length(var)){
  eval(parse(text = paste0(var[i], ' <- grep(var[', i, '], fls, value = T) %>% stack()')))
}
msk_lyr <- prec[[1]] * 0 + 1

# Solar radiation
xtr <- list.files('Z:/_colombiaETP/_data/_tif/ET_SolRad', full.names = TRUE) %>%
  grep('/et_', ., value = TRUE) %>%
  mixedsort() %>% 
  stack() %>% 
  raster::crop(., msk) %>%
  raster::mask(., msk)
xtr <- resample(x = xtr, y = msk_lyr, method = 'ngb')
xtr <- xtr*c(31,29,31,30,31,30,31,31,30,31,30,31)
nms <- names(xtr)
xtr <- unstack(xtr)
Map('writeRaster', x = xtr, paste0('../data/raster/clm/crn/', nms, '.asc'))
xtr <- stack(xtr)

# Calculating the ETP variables -------------------------------------------
etp <- 0.0013 * 0.408 * xtr * (tmean + 17) * (tmax - tmean - 0.0123 * prec) ^ 0.76
etp <- unstack(etp)
Map('writeRaster', x = etp, filename = paste0('../data/raster/clm/crn/etp_', 1:12, '.asc'), overwrite = TRUE)
etp <- stack(etp)

# To calculate bioclimatic variables --------------------------------------
biostack <- biovars(prec,tmin,tmax)
nms <- names(biostack)
biostack <- unstack(biostack)
Map('writeRaster', x = biostack, filename = paste0('../data/raster/clm/crn/', nms, '.asc'), overwrite = TRUE)

# Bioclimatic 20 Max cons months less than mm
prc_mtx <- as.matrix(prec)
gcmbiofolder <- '../data/raster/clm/crn/'
bio20ies <-  t(apply(prc_mtx, 1, cumDry))
mm <- c(40,50,100)
for(i in 1:3){ 
  bio20 <- msk_lyr
  values(bio20) <- bio20ies[,i]
  writeRaster(bio20,paste(gcmbiofolder,"/","bio20_",mm[i],"mm",sep=""),
              format = "ascii",
              overwrite = T, 
              NAflag=-9999)
} 

# Deficit stack
dfc_stk <- prec - etp
DefAndTemp <- cbind(as.matrix(dfc_stk),as.matrix(tmean),as.matrix(tmax))
biovalues <-  t(apply(DefAndTemp,1,cumTemp))

msk_lyr <- etp[[1]] * 0 + 1
biovalues <-  t(apply(DefAndTemp,1,cumTemp))

### BIO_21: Consecutive Months with less Prec than PET
bio21 <- msk_lyr
values(bio21) <- biovalues[,1]
writeRaster(bio21, paste('../data/raster/clm/crn/', 'bio21', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

### BIO_22: Sum of water deficit during dry season
bio22 <- msk_lyr
values(bio22) <- biovalues[,2]
writeRaster(bio22, paste('../data/raster/clm/crn/', 'bio22', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

### BIO_23: Mean temperature during growing season
bio23 <- msk_lyr
values(bio23) <- biovalues[,3]
writeRaster(bio23, paste('../data/raster/clm/crn/', 'bio23', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

### BIO_24: Max dry season temperature
bio24 <- msk_lyr
values(bio24) <- biovalues[,4]
writeRaster(bio24, paste('../data/raster/clm/crn/', 'bio24', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

ETPAndPrec <- cbind(as.matrix(etp),as.matrix(prec),as.matrix(tmean))
etpbios <-  t(apply(ETPAndPrec,1,etpvars))

lapply(1:nrow(etpbios), function(){
  
  
  
})

# Bio 25 Annual PET
bio25 <- msk_lyr
values(bio25) <- etpbios[,1]
plot(bio25)
writeRaster(bio25, paste('../data/raster/clm/crn/', 'bio25', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 26 PET seasonality (Coefficient of Variation)
bio26 <- msk_lyr
values(bio26) <- etpbios[,2]
plot(bio26)
writeRaster(bio26, paste('../data/raster/clm/crn/', 'bio26', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 27 MAX PET
bio27 <- msk_lyr
values(bio27) <- etpbios[,3]
writeRaster(bio27, paste('../data/raster/clm/crn/', 'bio27', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 28 Min PET
bio28 <- msk_lyr
values(bio28) <- etpbios[,4]
writeRaster(bio28, paste('../data/raster/clm/crn/', 'bio28', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 29 Range of PET (PETmax-PETmin)
bio29 <- msk_lyr
values(bio29) <- etpbios[,5]
writeRaster(bio29, paste('../data/raster/clm/crn/', 'bio29', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 30 PET of wettest quarter
bio30 <- msk_lyr
values(bio30) <- etpbios[,6]
writeRaster(bio30, paste('../data/raster/clm/crn/', 'bio30', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 31:	PET of driest quarter
bio31 <- msk_lyr
values(bio31) <- etpbios[,7]
writeRaster(bio31, paste('../data/raster/clm/crn/', 'bio31', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 32:	PET of warmest quarter
bio32 <- msk_lyr
values(bio32) <- etpbios[,8]
writeRaster(bio32, paste('../data/raster/clm/crn/', 'bio32', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)

# Bio 33:	PET of coldest quarter
bio33 <- msk_lyr
values(bio33) <- etpbios[,9]
writeRaster(bio33, paste('../data/raster/clm/crn/', 'bio33', sep = ''), format = 'ascii', overwrite = T, NAflag = -9999)












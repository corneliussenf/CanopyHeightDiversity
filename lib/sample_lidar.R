
### Libraries

library(raster)
library(tidyverse)
library(rgeos)
library(data.table)
library(fasterize)
library(sf)

### Set temp directory for raster analysis

rasterOptions(datatype = "FLT4S", 
              progress = "", 
              tmptime = 6)

### Loop through landscapes

for (site in c("berchtesgaden", "boehmerwald")) {
  
  ### Load data

  disturbance_onset <- raster(paste0("data/disturbance/disturbance_mmu5pixels_", site, ".dat"))
  disturbance_onset[!disturbance_onset %in% c(1985:2016)] <- 0
  forestmask <- raster(paste0("data/forest/forestmask_1985_", site, ".dat"))
  forestmask[forestmask == 255] <- NA

  if (site == "boehmerwald") {
    dist_after_lidar <- 2012
  } else if (site == "berchtesgaden") {
    dist_after_lidar <- 2009
  }
  
  disturbance_onset[disturbance_onset > dist_after_lidar] <- 0
  
  lidar_height <- raster(paste0("data/lidar/", site, "/lidar_chm_", site, "_utm33.tif"))
  lidar_elev <- raster(paste0("data/lidar/", site, "/lidar_dtm_", site, "_utm33.tif"))

  stratum <- shapefile(paste0("data/sites/", site, "_no_intervention.shp"))

  ### Create sample-raster

  sampleraster <- crop(mask(disturbance_onset > 0, stratum), stratum)
  idraster <- sampleraster
  values(idraster) <- 1:ncell(sampleraster)
  sampleraster_points <- rasterToPoints(sampleraster, spatial = TRUE)
  sampleraster_points$disturbance <- sampleraster_points$layer
  sampleraster_points$disturbance_onset <- raster::extract(disturbance_onset, sampleraster_points)
  sampleraster_points$disturbance_onset <- ifelse(sampleraster_points$disturbance_onset %in% c(1985:2016), sampleraster_points$disturbance_onset, 0)
  sampleraster_points$layer <- NULL
  sampleraster_points$cell <- raster::extract(idraster, sampleraster_points)
  sampleraster_points$xcoord <- sampleraster_points@coords[, 1]
  sampleraster_points$ycoord <- sampleraster_points@coords[, 2]
  sampleraster_points_buffer <- gBuffer(sampleraster_points, byid = TRUE,
                                        id = sampleraster_points$cell, width = 15,
                                        capStyle = "SQUARE")
  
  sampleraster_1m <- fasterize(st_as_sf(sampleraster_points_buffer), raster = lidar_height, field = "cell")

  forestraster <- crop(mask(forestmask, stratum), stratum)
  forest_points <- as.data.frame(rasterToPoints(forestraster, fun = function(x) !is.na(x)))
  names(forest_points) <- c("xcoord", "ycoord", "forest")
  sampleraster_points_buffer <- merge(sampleraster_points_buffer, forest_points, by = c("xcoord", "ycoord"))

  ### Extract Lidar CHM

  lidar_height_sample <- data.table(chm = values(lidar_height),
                                    cellnumber = values(sampleraster_1m)) %>%
    filter(., !is.na(cellnumber))

  ### Functions
  
  # Structural diversity metrics
  get_diversity_metrics <- function (cellnumbers) {
    
    ind <- which(sampleraster_points_buffer$cell %in% cellnumbers)
    
    # if (length(ind) > 0) {
    if (mean(cellnumbers %in% sampleraster_points_buffer$cell) > 0) {
      
      lidar_height_sample_list <- lidar_height_sample %>%
        filter(cellnumber %in% cellnumbers) %>%
        split(.$cellnumber) %>%
        map(~ as.double(.$chm))
      
      height_abundances <- lidar_height_sample_list %>%
        map2(.y = paste0("site_", 1:length(lidar_height_sample_list)),
             ~ data.frame(height = .) %>%
               mutate(height = ifelse(height < 2, NA, height)) %>%
               mutate(., height_class = tryCatch(cut(height, seq(0, 100, 1), include.lowest = TRUE),
                                                 error = function(e) return(NA))) %>%
               group_by(., height_class) %>%
               summarize(., pi = ifelse(length(height[!is.na(height)]) == 0, NA, length(height[!is.na(height)]))) %>%
               ungroup(.) %>%
               dplyr::rename(!!.y := pi)) %>%
        plyr::join_all(., by = "height_class", type = "full") %>%
        filter(!is.na(height_class)) %>%
        column_to_rownames("height_class") %>%
        mutate_all(function(x) ifelse(is.na(x), 0, x))
      
      if (nrow(height_abundances) > 0) {
      
        alpha0 = vegetarian::d(t(height_abundances), lev = "alpha", q = 0)
        beta0 = vegetarian::d(t(height_abundances), lev = "beta", q = 0)
        gamma0 = vegetarian::d(t(height_abundances), lev = "gamma", q = 0)
        
        alpha1 = vegetarian::d(t(height_abundances), lev = "alpha", q = 1)
        beta1 = vegetarian::d(t(height_abundances), lev = "beta", q = 1)
        gamma1 = vegetarian::d(t(height_abundances), lev = "gamma", q = 1)
        
        alpha2 = vegetarian::d(t(height_abundances), lev = "alpha", q = 2)
        beta2 = vegetarian::d(t(height_abundances), lev = "beta", q = 2)
        gamma2 = vegetarian::d(t(height_abundances), lev = "gamma", q = 2)
        
        return(data.frame(alpha0 = alpha0, 
                          alpha1 = alpha1, 
                          alpha2 = alpha2, 
                          beta0 = beta0, 
                          beta1 = beta1, 
                          beta2 = beta2, 
                          gamma0 = gamma0, 
                          gamma1 = gamma1, 
                          gamma2 = gamma2))
        
      } else {
        
        return(data.frame(alpha0 = 0, 
                          alpha1 = 1, 
                          alpha2 = NA, 
                          beta0 = NA, 
                          beta1 = 1, 
                          beta2 = NA, 
                          gamma0 = 0, 
                          gamma1 = 1, 
                          gamma2 = NA))
        
      }
      
    } else {
      
      return(data.frame(alpha0 = NA, 
                        alpha1 = NA, 
                        alpha2 = NA, 
                        beta0 = NA, 
                        beta1 = NA, 
                        beta2 = NA, 
                        gamma0 = NA, 
                        gamma1 = NA, 
                        gamma2 = NA))
      
    }
    
  }
  
  # Topographic metrics
  get_topography_metrics <- function (cellnumbers, size) {
    
    if (mean(cellnumbers %in% sampleraster_points_buffer$cell) > 0) {
      
      shp <- subset(sampleraster_points_buffer, sampleraster_points_buffer$cell %in% cellnumbers)
      lidar_elev_tmp <- crop(lidar_elev, shp)
      tri <- terrain(lidar_elev_tmp, opt = "TRI")
      slope <- terrain(lidar_elev_tmp, opt = "slope", unit = "degrees")
      
      varcof <- function(x, na.rm) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
      
      topo <- data.frame(elevation = values(lidar_elev_tmp),
                         tri = values(tri),
                         slope = values(slope)) %>%
        summarize_all(.funs = funs(mean, varcof, .args = list(na.rm = TRUE)))
      
      return(topo)
      
    } else {
      
      return(data.frame(elevation_mean = 0, 
                        tri_mean = 0,
                        slope_mean = 0,
                        elevation_varcof = 0,
                        tri_varcof = 0,
                        slope_varcof = 0))
      
    }
    
  }
  
  # Disturbance metrics
  
  get_disturbance_metrics <- function (cellnumbers, size) {
    
    ind <- which(sampleraster_points_buffer@data$cell %in% cellnumbers)
    dist_prop <- sum(sampleraster_points_buffer@data$disturbance[ind], na.rm = TRUE) / 
      sum(sampleraster_points_buffer@data$forest[ind] == 1, na.rm = TRUE)
    dist_onset <- sampleraster_points_buffer@data$disturbance_onset[ind]
    dist_onset_median <- median(dist_onset[dist_onset != 0 & !is.na(dist_onset) & dist_onset != 1985])
    
    shp <- subset(sampleraster_points_buffer, sampleraster_points_buffer$cell %in% cellnumbers)
    
    sampleraster_tmp <- crop(sampleraster, shp)
    sampleraster_tmp[sampleraster_tmp == 0] <- NA
    
    if (!is.nan(dist_prop) & dist_prop > 0) {
      
      return(data.frame(disturbance_rate = dist_prop, 
                        disturbance_onset = dist_onset_median))
      
    } else {
      
      return(data.frame(disturbance_rate = 0, 
                        disturbance_onset = 0))
    }
    
    
  }
  
  # Forest metrics
  
  get_forest_metrics <- function (cellnumbers, size) {
    
    ind <- which(sampleraster_points_buffer@data$cell %in% cellnumbers)
    for_prop <- mean(sampleraster_points_buffer@data$forest[ind] == 1, na.rm = TRUE)
    
    shp <- subset(sampleraster_points_buffer, sampleraster_points_buffer$cell %in% cellnumbers)
    
    forestmask_tmp <- crop(forestmask, shp)
    forestmask_tmp[forestmask_tmp == 2] <- NA
    
    if (!is.nan(for_prop) & for_prop > 0) {
      
      
      return(data.frame(forest_rate = for_prop))
      
    } else {
      
      return(data.frame(forest_rate = 0))
    }
    
    
  }
  
  ### Sample sub-landscapes and calculate structural diversity metrics + predictors
  
  n <- 150
  sizes <- c(10, 15, 20, 25, 30, 35)
  diversity <- vector("list", length(sizes))
  j <- 0
  
  for (size in sizes) {
    
    j <- j + 1
    
    neighborhood_matrix <- matrix(rep(1, size ^ 2), ncol = size, nrow = size)
    neighborhood_matrix[ceiling(size / 2), ceiling(size / 2)] <- 0
    
    ii <- sample(unique(sampleraster_points_buffer@data$cell), n)
    
    div <- vector("list", length = length(ii))
    dist <- vector("list", length = length(ii))
    topo <- vector("list", length = length(ii))
    forest <- vector("list", length = length(ii))
    
    k <- 0
    
    for (i in ii) {
      
      k <- k + 1
      
      print(paste0("Size: ", size, "; Sublandscape: ", k))
      
      cellnumbers <- adjacent(sampleraster, 
                              i, 
                              directions = neighborhood_matrix, 
                              pairs = FALSE, 
                              include = TRUE, 
                              id = TRUE)
      
      div[[k]] <- get_diversity_metrics(cellnumbers)
      dist[[k]] <- get_disturbance_metrics(cellnumbers)
      topo[[k]] <- get_topography_metrics(cellnumbers, size)
      forest[[k]] <- get_forest_metrics(cellnumbers, size)
      
    }
    
    div <- do.call("rbind", div)
    dist <- do.call("rbind", dist)
    topo <- do.call("rbind", topo)
    forest <- do.call("rbind", forest)
    
    windowsize <- rep(size, length(ii))
    
    diversity[[j]] <- cbind(div, dist, topo, forest, windowsize)
    
  }
  
  diversity <- bind_rows(diversity)
  diversity$sublandscape <- 1:n
  assign(paste0("diversity_", site), diversity)
  
  save(list = paste0("diversity_", site), file = paste0("data/diversity_", site, ".RData"))
  
}


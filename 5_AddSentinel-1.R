# ============================================================================
# SENTINEL-1 VARIABLE EXTRACTION
# ============================================================================
# Description:
#   This script extracts Sentinel-1 radar variables and texture metrics for
#   harvester grid cells. The output is a grid file with the extracted SAR
#   variables.
#
# Author: Muhammad Lukman Yazid (Luk)
# Thesis: Timber Volume Estimation Using Harvester Data and Multisource Remote 
#          Sensing: Integrating Airborne Laser Scanning, Sentinel Imagery, 
#          Topographic, and Stand Data
# ============================================================================

# 1. CLEAN ENVIRONMENT =======================================================

rm(list = ls())
gc()

# Load required libraries
library(terra)         # Raster handling
library(sf)            # Vector/spatial data
library(dplyr)         # Data manipulation
library(glcm)          # Texture metrics
library(exactextractr) # Raster extraction
library(lidR)          # LiDAR
library(future)        # Parallel processing

# Set up parallel processing
n_cores <- parallel::detectCores()
plan(multisession, workers = n_cores)

# 2. LOAD DATA ==============================================================

# Relative paths (replace with your project structure)
grid_path <- "data/coreGrid_ALS_DTM_S1_S2_S_NW.gpkg"
sar_file <- "data/S1.tif"
output_path <- "results/coreGrid_ALS_DTM_SAR_NW2.gpkg"

# Load harvester grid shapefile
grid <- st_read(grid_path)

# Remove duplicates and unnecessary rows (if applicable)
grid <- grid[!duplicated(grid$gridIndex), ]
grid <- grid[-c(108:116)] # Optional: remove unwanted rows

# Load SAR raster
sar_stack <- rast(sar_file)

# Rename bands (if required)
names(sar_stack) <- c("VV", "VH", "VV_VH_ratio", "VV_minus_VH",
                      "VV_VH_norm_diff", "VV_texture_glcm_contrast",
                      "VV_texture_glcm_entropy", "VV_texture_glcm_homogeneity")[1:nlyr(sar_stack)]

# 3. EXTRACT SAR METRICS ====================================================

# Extract mean values for each grid cell
sar_metrics <- exact_extract(sar_stack, grid, 'mean')

# 4. CALCULATE TEXTURE METRICS ==============================================

if ("VV" %in% names(sar_stack)) {
  # Convert VV band to raster for GLCM calculations
  vv_raster <- raster::raster(sar_stack[["VV"]])
  
  # Compute GLCM texture features
  vv_texture <- glcm::glcm(vv_raster, window = c(3,3),
                           statistics = c("contrast", "entropy", "homogeneity"))
  
  names(vv_texture) <- paste0("VV_texture_", names(vv_texture))
  
  # Convert texture output to terra format
  vv_texture_terra <- terra::rast(vv_texture)
  
  # Extract mean texture metrics for each grid cell
  texture_metrics <- exact_extract(vv_texture_terra, grid, 'mean')
  
  # Combine SAR and texture metrics
  all_metrics <- bind_cols(as.data.frame(sar_metrics), as.data.frame(texture_metrics))
} else {
  all_metrics <- as.data.frame(sar_metrics)
}

# 5. JOIN METRICS TO GRID ===================================================

grid_with_metrics <- bind_cols(grid, all_metrics)

# 6. SAVE OUTPUT =============================================================

st_write(grid_with_metrics, output_path, delete_dsn = TRUE)

# Optional: summary of a key metric
if ("VV_VH_ratio" %in% names(grid_with_metrics)) {
  print(summary(grid_with_metrics$VV_VH_ratio))
}

# Optional: view final dataset
# View(grid_with_metrics)

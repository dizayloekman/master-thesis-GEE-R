# ============================================================================
# SENTINEL-2 VARIABLE EXTRACTION
# ============================================================================
# Description:
#   This script extracts Sentinel-2 band reflectivity and vegetation indices
#   for harvester grid cells. The output is a grid file with the extracted
#   Sentinel-2 variables.
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
library(terra)         # Raster data handling
library(sf)            # Vector/spatial data
library(raster)        # Raster support
library(caret)         # Modeling framework
library(ggplot2)       # Visualization
library(tmap)          # Mapping
library(future)        # Parallel processing
library(dplyr)         # Data manipulation
library(exactextractr) # Raster extraction
library(lidR)          # LiDAR data handling

# Set up parallel processing
n_cores <- parallel::detectCores()
plan(multisession, workers = n_cores)

# 2. LOAD DATA ==============================================================

# Relative paths (replace with your project structure)
shapefile_path <- "data/coreGrid_ALS_DTM_SAR_NW2.gpkg"
sentinel_path <- "data/S2.tif"
output_file <- "results/GridHarvesterSentinel.gpkg"

# Load Sentinel-2 raster
sentinel2 <- stack(sentinel_path)  # Multi-band raster

# Load harvester grid shapefile
grid <- st_read(shapefile_path)

# Ensure CRS match
if (!st_crs(grid) == crs(sentinel2)) {
  sentinel2 <- projectRaster(sentinel2, crs = st_crs(grid)$proj4string)
}

# 3. PREPARE SENTINEL-2 VARIABLES ===========================================

# Clip Sentinel-2 to grid extent
sentinel_clipped <- crop(sentinel2, extent(grid))

# Band reflectivity variables (Bands 2-12)
band_vars <- sentinel_clipped[[2:12]]
names(band_vars) <- paste0("B", 2:12)

# Calculate vegetation indices
rvi       <- sentinel_clipped[[8]] / sentinel_clipped[[4]]
rvire     <- sentinel_clipped[[8]] / sentinel_clipped[[5]]
dvi       <- sentinel_clipped[[8]] - sentinel_clipped[[4]]
ndvi      <- (sentinel_clipped[[8]] - sentinel_clipped[[4]]) / 
  (sentinel_clipped[[8]] + sentinel_clipped[[4]])
ndvire1   <- (sentinel_clipped[[8]] - sentinel_clipped[[5]]) / 
  (sentinel_clipped[[8]] + sentinel_clipped[[5]])
ndvire2   <- (sentinel_clipped[[8]] - sentinel_clipped[[6]]) / 
  (sentinel_clipped[[8]] + sentinel_clipped[[6]])
mndvi     <- (sentinel_clipped[[8]] - sentinel_clipped[[4]]) / 
  (sentinel_clipped[[8]] + sentinel_clipped[[4]] - 2 * sentinel_clipped[[2]])
mndvire   <- (sentinel_clipped[[8]] - sentinel_clipped[[5]]) / 
  (sentinel_clipped[[8]] + sentinel_clipped[[5]] - 2 * sentinel_clipped[[2]])
nd11      <- (sentinel_clipped[[8]] - sentinel_clipped[[11]]) / 
  (sentinel_clipped[[8]] + sentinel_clipped[[11]])

# Stack vegetation indices
veg_indices <- stack(rvi, rvire, dvi, ndvi, ndvire1, ndvire2, mndvi, mndvire, nd11)
names(veg_indices) <- c("RVI", "RVIre", "DVI", "NDVI", "NDVIre1", 
                        "NDVIre2", "mNDVI", "mNDVIre", "ND11")

# 4. EXTRACT VARIABLES FOR GRID CELLS ======================================

# Mean values for each grid cell
extracted_band_vars <- extract(band_vars, grid, fun = mean, df = TRUE, na.rm = TRUE)
extracted_veg_indices <- extract(veg_indices, grid, fun = mean, df = TRUE, na.rm = TRUE)

# Combine extracted data
extracted_data <- cbind(grid, extracted_band_vars, extracted_veg_indices)

# Remove rows with missing values in key band (e.g., B2)
extracted_data_clean <- extracted_data[!is.na(extracted_data$B2), ]

# 5. SAVE OUTPUT ============================================================

# CSV for reference
write.csv(extracted_data_clean, "results/GridHarvesterSentinel_noNA.csv", row.names = FALSE)

# Save as GeoPackage
st_write(extracted_data_clean, output_file, delete_dsn = TRUE)

# Optional: View output
# View(extracted_data_clean)

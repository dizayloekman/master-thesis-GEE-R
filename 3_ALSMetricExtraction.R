# ============================================================================
# LASMETRIC EXTRACTION
# ============================================================================
# Description:
#   This script extracts lasmetrics from ALS point clouds for
#   harvester grid cells.
#
# Author: Muhammad Lukman Yazid (Luk)
# Thesis: Timber Volume Estimation Using Harvester Data and Multisource Remote 
#          Sensing: Integrating Airborne Laser Scanning, Sentinel Imagery, 
#          Topographic, and Stand Data
# ============================================================================

# 1. CLEAN ENVIRONMENT =======================================================

rm(list = ls())
gc()
options(scipen = 999)  # Disable scientific notation

# Load required libraries
suppressPackageStartupMessages({
  library(sf)           # Spatial vector data
  library(lidR)         # ALS point cloud processing
  library(stringr)      # String operations
  library(terra)        # Raster processing
  library(dplyr)        # Data manipulation
  library(future)       # Parallel processing
  library(lidRmetrics)  # LiDAR metrics
  library(tidyverse)    # Data wrangling
})

# 2. SET UP PARALLEL PROCESSING =============================================

n_cores <- parallel::detectCores()
plan(multisession, workers = n_cores)

# 3. DEFINE PATHS ============================================================

las_folder      <- "data/ALS_normalized/"      # Folder with .las/.laz files
grid_shapefile  <- "data/coreCellsHarvesterALS.gpkg"  # Harvester core grid
output_shapefile <- "results/coreGridMetrics.gpkg"
output_rdata    <- "results/coreGridMetrics.rds"

# 4. LOAD HARVESTER GRID ====================================================

grid <- st_read(grid_shapefile)
grid <- grid %>% mutate(fid = row_number())  # Add unique ID

# 5. LOAD AND PROCESS ALS FILES =============================================

# List all LAS/LAZ files
laz_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, recursive = TRUE)

# Create LAScatalog
ctg <- readLAScatalog(laz_files)
opt_chunk_size(ctg) <- 2000  # Optional: adjust for memory/performance

# Extract standard LiDAR metrics for each grid cell
ctg_metrics <- plot_metrics(ctg, .stdmetrics_z, grid)

message("??? Lasmetric extraction completed for all LAS files.")

# 6. FILTER CORE GRID ========================================================

# Calculate grid cell area in hectares
grid$area <- as.numeric(st_area(grid)) / 10000

# Filter grid cells with area >= 0.0256 ha
core_grid <- grid %>% filter(area >= 0.0256)

# Save core grid metrics
st_write(core_grid, output_shapefile, append = FALSE)
saveRDS(core_grid, output_rdata)

# 7. ADD ALS PROJECT METADATA =================================================

# Example project shapefiles (replace with your relative paths)
shapefiles <- list(
  Valdres2013 = "data/ALS_projects/Valdres2013.shp",
  NordreLand  = "data/ALS_projects/NordreLand2016.shp"
  # Add other projects as needed
)

# Load and harmonize CRS
shapefile_list <- lapply(shapefiles, st_read, quiet = TRUE)
for (name in names(shapefile_list)) {
  shapefile_list[[name]] <- shapefile_list[[name]] %>%
    st_transform(st_crs(core_grid)) %>%
    mutate(
      ALSproject = name,
      ALSyear = as.integer(str_extract(shapefiles[[name]], "\\d{4}"))
    )
}

# Combine all projects
combined_projects <- bind_rows(shapefile_list)

# Spatial join: assign ALS project info to grid
grid_with_metadata <- st_join(core_grid, combined_projects, join = st_intersects, left = TRUE)

# Keep only relevant columns
grid_with_metadata <- grid_with_metadata %>%
  select(-matches("SOSI_ID|OBJTYPE|KARTID|KARTSERIE|Shape_Leng|Shape_Area")) 

# Save final dataset
saveRDS(grid_with_metadata, "results/coreGridMetrics_withALS.rds")
st_write(grid_with_metadata, "results/coreGridMetrics_withALS.gpkg", append = FALSE)

message("??? Core grid metrics with ALS project metadata saved successfully.")

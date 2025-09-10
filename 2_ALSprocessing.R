# ============================================================================
# ALS DATA PROCESSING: REMOVAL OF OUTLIERS, NOISE, AND NORMALIZATION
# ============================================================================
# Description:
#   This script processes raw ALS point clouds (.las/.laz) by removing
#   duplicates, classifying and removing noise, normalizing heights, and
#   setting a standard CRS for all files.
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
suppressPackageStartupMessages({
  library(sf)           # Spatial vector data
  library(lidR)         # ALS processing
  library(stringr)      # String manipulation
  library(terra)        # Raster operations
  library(dplyr)        # Data manipulation
  library(lidRmetrics)  # ALS metrics
  library(tools)        # File path utilities
  library(tidyverse)    # Data wrangling
})

# 2. SET UP PARALLEL PROCESSING =============================================

n_cores <- parallel::detectCores()
plan(multisession, workers = n_cores)
set_lidr_threads(n_cores)

# 3. DEFINE PATHS ============================================================

input_dir  <- "data/raw_ALS/"         # Input LAS/LAZ files
output_dir <- "data/normalized_ALS/"  # Output normalized files

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 4. DEFINE TARGET CRS ======================================================

# UTM zone 32N (EPSG:25832)
crs_target <- st_crs(25832)

# 5. PROCESS EACH LAS FILE ==================================================

las_files <- list.files(input_dir, pattern = "\\.(las|laz)$", full.names = TRUE)

for (las_file in las_files) {
  
  # Read LAS file
  las <- readLAS(las_file)
  
  # Skip invalid or empty files
  if (is.null(las) || npoints(las) == 0) {
    message(paste("Skipping empty or invalid LAS file:", las_file))
    next
  }
  
  # Remove duplicated points
  las <- filter_duplicates(las)
  
  # Classify noise using Statistical Outlier Removal (SOR)
  las <- classify_noise(las, algorithm = sor(k = 10, m = 3, quantile = FALSE))
  
  # Remove noise points
  las <- filter_poi(las, Classification != LASNOISE)
  
  # Normalize heights
  las_normalized <- normalize_height(las, knnidw())
  
  # Transform to target CRS
  las_normalized <- st_transform(las_normalized, crs_target)
  
  # Define output path
  output_file <- file.path(output_dir,
                           paste0(file_path_sans_ext(basename(las_file)), "_norm.laz"))
  
  # Save normalized LAS file
  writeLAS(las_normalized, output_file)
  
  message(paste("Processed LAS file saved to:", output_file))
}

message("??? ALS normalization process completed for all LAS files.")

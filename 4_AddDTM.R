# ============================================================================
# DTM VARIABLE EXTRACTION
# ============================================================================
# Description:
#   This script extracts DTM-based terrain variables (TPI and Qfit metrics)
#   for harvester grid cells. The output is a grid file with added DTM metrics.
#
# Author: Muhammad Lukman Yazid (Luk)
# Modified from Lennart Noordermeer's code
# Thesis: Timber Volume Estimation Using Harvester Data and Multisource Remote 
#          Sensing: Integrating Airborne Laser Scanning, Sentinel Imagery, 
#          Topographic, and Stand Data
# ============================================================================

# 1. CLEAN ENVIRONMENT =======================================================

rm(list = ls())
gc()

# Load required libraries
suppressPackageStartupMessages({
  library(sf)           # Vector/spatial data
  library(terra)        # Raster data handling
  library(dplyr)        # Data manipulation
  library(exactextractr)# Raster extraction
  library(MultiscaleDTM)# Terrain metrics
  library(future)       # Parallel processing
  library(stringr)      # String operations
})

# Set up parallel processing
n_cores <- parallel::detectCores()
plan(multisession, workers = n_cores)

# 2. DEFINE PATHS ============================================================

# Relative paths (replace with your project structure)
grid_path   <- "data/coreGrid_ALS_NW.gpkg"
dtm_folder  <- "data/DTM/"        # Folder containing DTM .tif files
output_path <- "results/coreGrid_ALS_DTM_NW.gpkg"

# 3. LOAD GRID ===============================================================

grid <- st_read(grid_path, quiet = TRUE)
grid$fid <- 1:nrow(grid)  # Unique ID for merging metrics later

# List DTM files
dtm_files <- list.files(dtm_folder, pattern = "\\.tif$", full.names = TRUE)

# 4. DEFINE AGGREGATION FACTORS =============================================

agg_factors <- c(5, 10, 20)          # Aggregation factors for multi-scale metrics
suffixes    <- as.character(agg_factors)

# Initialize list to store results
all_results <- vector("list", length(dtm_files))
names(all_results) <- basename(dtm_files)

# 5. PROCESS DTM FILES =======================================================

for (dtm_file in dtm_files) {
  cat("Processing:", basename(dtm_file), "\n")
  
  # Load DTM and reproject to grid CRS
  dtm <- rast(dtm_file)
  dtm <- project(dtm, st_crs(grid)$wkt)
  
  # Identify overlapping grid cells
  grid_sub <- grid[st_intersects(grid, st_as_sfc(st_bbox(dtm)), sparse = FALSE), ]
  if (nrow(grid_sub) == 0) next
  
  # Initialize dataframe to collect metrics
  metrics_df <- data.frame(fid = grid_sub$fid)
  
  # Loop over aggregation scales
  for (i in seq_along(agg_factors)) {
    factor <- agg_factors[i]
    suffix <- suffixes[i]
    
    cat(" Aggregation factor:", factor, "\n")
    
    # Aggregate DTM
    dtm_agg <- terra::aggregate(dtm, fact = factor, fun = mean)
    
    # Compute TPI
    tpi <- TPI(dtm_agg, shape = "rectangle", na.rm = TRUE)
    
    # Compute Qfit metrics
    qmetrics <- Qfit(dtm_agg, w = c(3, 3), unit = "degrees",
                     metrics = c("elev", "qslope", "qaspect", "qeastness", "qnorthness",
                                 "profc", "planc", "twistc", "meanc", "maxc", "minc", "features"),
                     na.rm = TRUE)
    
    # Extract Qfit metrics
    for (metric in names(qmetrics)) {
      cat("  Extracting:", metric, "\n")
      values <- exact_extract(qmetrics[[metric]], grid_sub, 'mean')
      metrics_df[[paste0(metric, "_", suffix)]] <- values
    }
    
    # Extract TPI
    metrics_df[[paste0("TPI_", suffix)]] <- exact_extract(tpi, grid_sub, 'mean')
  }
  
  # Merge metrics with subset of grid
  result_grid <- left_join(grid_sub, metrics_df, by = "fid")
  all_results[[basename(dtm_file)]] <- result_grid
}

# 6. COMBINE RESULTS ========================================================

final_result <- do.call(rbind, all_results)
final_result$fid <- 1:nrow(final_result)  # Ensure unique IDs

# 7. SAVE OUTPUT ============================================================

st_write(final_result, output_path, delete_dsn = TRUE)
cat("??? DTM variable extraction complete.\n")

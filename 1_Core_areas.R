# ============================================================================
# CORE GRID HARVESTER POLYGON SELECTION
# ============================================================================
# Description:
#   This script selects core polygons and stems for different inventory methods,
#   ensuring that the same areas are used for fair comparison. Exported ITC core
#   segments also include ALS metrics.
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

options(digits = 4)
options(scipen = 999)

# 2. LOAD REQUIRED LIBRARIES =================================================

suppressPackageStartupMessages({
  library(sf)         # Spatial vector data
  library(stringr)    # String manipulation
  library(sp)         # Spatial functions
  library(plyr)       # Data manipulation
  library(dplyr)      # Data manipulation
  library(data.table) # Efficient data handling
})

# 3. DEFINE PATHS ============================================================

grid_file <- "data/gridData_compiled.gpkg"   # GeoPackage of polygons
pts_file  <- "data/stemData_compiled.gpkg"   # GeoPackage of stem points
output_pdf <- "results/Core_Area/1_core_area.pdf"
output_gpkg <- "results/Core_Area/coreCells.gpkg"

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_gpkg), recursive = TRUE, showWarnings = FALSE)

# 4. READ ALL LAYERS ========================================================

grid_layers <- st_layers(grid_file)$name
pts_layers  <- st_layers(pts_file)$name

grid <- lapply(grid_layers, function(layer) st_read(grid_file, layer = layer, quiet = TRUE))
pts  <- lapply(pts_layers,  function(layer) st_read(pts_file, layer = layer, quiet = TRUE))

# 5. FUNCTION TO NORMALIZE COLUMN NAMES ====================================

normalize_names <- function(df) {
  names(df) <- iconv(names(df), to = "ASCII//TRANSLIT")
  df <- df[, !is.na(names(df))]  # Remove columns with NA names
  return(df)
}

# 6. PDF OUTPUT FOR PLOTS ====================================================

pdf(file = output_pdf)

n_grids <- length(grid)

for (i in seq_len(n_grids)) {
  
  sq <- grid[[i]]
  
  # Identify non-border polygons using a small buffer
  buf <- st_buffer(sq, 1) %>% st_union() %>% st_buffer(-1.1)
  idx <- lengths(st_within(sq, buf)) > 0  # TRUE if polygon is not on the border
  
  # Keep only stands with sufficient area (>= 2 da)
  if (sum(idx) >= 8) {
    
    sq <- sq[idx, ]  # Subset core cells
    
    # Plot polygons and stems
    plot(sq$geom, col = "green", main = unique(sq$operationDir))
    plot(pts[[i]]$geom, col = "#00000080", pch = 20, add = TRUE)
    
    # Export or append core ABA cells
    if (!file.exists(output_gpkg)) {
      st_write(sq, output_gpkg)
    } else {
      existing <- st_read(output_gpkg)
      existing <- existing[existing$operationDir != unique(sq$operationDir), ]
      
      sq <- normalize_names(sq)
      
      # Remove count columns
      existing <- existing[, !grepl("count", names(existing), ignore.case = TRUE)]
      sq       <- sq[, !grepl("count", names(sq), ignore.case = TRUE)]
      
      combined <- rbind(existing, sq)
      st_write(combined, output_gpkg, delete_dsn = TRUE)
    }
    
  } else {
    message(paste("Stand", i, "removed because area < 2 da"))
  }
  
  message(paste("Completed stand", i, "out of", n_grids))
}

dev.off()
message("??? Core grid selection completed. PDF and GPKG exported.")

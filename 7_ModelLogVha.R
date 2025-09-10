# ============================================================================
# VOLUME MODELING USING RANDOM FOREST
# ============================================================================
# Description:
#   This script estimates log wood volume using Random Forest and multiple
#   sets of remote sensing features. It performs feature selection, model 
#   training, evaluation, and visualization of results.
#
# Author: Muhammad Lukman Yazid (Luk)
# Thesis: Timber Volume Estimation Using Harvester Data and Multisource Remote 
#          Sensing: Integrating Airborne Laser Scanning, Sentinel Imagery, 
#          Topographic, and Stand Data
# ============================================================================

# 1. SETUP ENVIRONMENT =======================================================

# Clear console and workspace
cat("\014")  # Clear console
rm(list = ls())  # Remove all objects
gc()  # Free memory

# Load required packages
library(randomForestSRC) # Random Forest modeling
library(randomForest)
library(caret)           # Modeling framework
library(dplyr)           # Data manipulation
library(ggplot2)         # Plotting
library(tidyr)           # Data reshaping
library(readr)           # CSV reading/writing
library(sf)              # Spatial data handling
library(pdp)             # Partial dependence plots
library(future)          # Parallel processing
library(viridis)
library(RColorBrewer)

# Set up parallel processing
n_cores <- parallel::detectCores()
plan(multisession, workers = n_cores)

# 2. LOAD AND PREPARE DATA ===================================================

# Load cleaned dataset (ensure sensitive paths/data are removed)
data <- st_read("data/clean_data.gpkg")  # Replace with relative path

# Define target variable
prediction_vars <- c("logV_ha")

# Define feature groups
ALS_vars <- c("zmax","zmin","zmean","zvar","zsd","zcv","zskew","zkurt",
              paste0("zq", seq(5,95,5)),
              "pzabovemean","pzabove2","pzabove5",
              "ziqr","zMADmean","zMADmedian","CRR","zentropy",
              paste0("zpcum",1:9))

terrain_vars <- c("elev_5","qslope_5","qaspect_5","qeastness_5","qnorthness_5",
                  "profc_5","planc_5","twistc_5","meanc_5","maxc_5","minc_5",
                  "features_5","TPI_5","elev_10","qslope_10","qaspect_10",
                  "qeastness_10","qnorthness_10","profc_10","planc_10",
                  "twistc_10","meanc_10","maxc_10","minc_10","features_10",
                  "TPI_10","elev_20","qslope_20","qaspect_20","qeastness_20",
                  "qnorthness_20","profc_20","planc_20","twistc_20","meanc_20",
                  "maxc_20","minc_20","features_20","TPI_20")

s1_vars <- c("mean.VV","mean.VH","mean.VV_VH_ratio","mean.VV_VH_norm_diff",
             "mean.VV_texture_glcm_contrast","mean.VV_texture_glcm_entropy",
             "mean.VV_texture_glcm_homogeneity")

s2_vars <- c("B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12",
             "RVI","RVIre","DVI","NDVI","NDVIre1","NDVIre2","mNDVI",
             "mNDVIre","ND11")

stand_vars <- c("BONITET","HOGSTKLASSE","ALDER")

# Define feature combinations
feature_sets <- list(
  All = c(ALS_vars, s1_vars, s2_vars, terrain_vars, stand_vars),
  ALS_Terrain_Stand = c(ALS_vars, terrain_vars, stand_vars),
  ALS_S1_Terrain_Stand = c(ALS_vars, s1_vars, terrain_vars, stand_vars),
  ALS_S2_Terrain_Stand = c(ALS_vars, s2_vars, terrain_vars, stand_vars),
  S1_S2_Terrain_Stand = c(s1_vars, s2_vars, terrain_vars, stand_vars),
  S1_Terrain_Stand = c(s1_vars, terrain_vars, stand_vars),
  S2_Terrain_Stand = c(s2_vars, terrain_vars, stand_vars),
  ALS = ALS_vars,
  S1 = s1_vars,
  S2 = s2_vars
)

# Verify all required variables are present
all_vars <- c(prediction_vars, ALS_vars, s1_vars, s2_vars, terrain_vars, stand_vars)
missing_vars <- setdiff(all_vars, names(data))
if(length(missing_vars) > 0) {
  stop("Missing variables: ", paste(missing_vars, collapse = ", "))
}

# Subset data to required variables and convert to data.frame
data <- data[, all_vars]
data_sf <- data  # keep spatial copy
data <- as.data.frame(data)
data$geom <- NULL  # remove geometry column

# Convert stand variables to numeric
data$BONITET <- as.numeric(data$BONITET)
data$HOGSTKLASSE <- as.numeric(data$HOGSTKLASSE)

# 3. TRAIN-TEST SPLIT =======================================================

set.seed(123)
train_idx <- createDataPartition(data$logV_ha, p = 0.8, list = FALSE)
train_data <- data[train_idx, ]
test_data <- data[-train_idx, ]

# 4. RANDOM FOREST TUNING ===================================================

tune_result <- tune.rfsrc(
  formula = logV_ha ~ .,
  data = train_data,
  mtryStart = floor(sqrt(ncol(data) - 1)),
  nodesizeTry = c(1:9, seq(10,100,5)),
  ntreeTry = 100,
  stepFactor = 1.25,
  improve = 1e-3,
  strikeout = 3,
  maxIter = 25,
  trace = TRUE,
  doBest = TRUE
)

optimal_mtry <- tune_result$optimal[["mtry"]]
optimal_nodesize <- tune_result$optimal[["nodesize"]]
print(tune_result)

# 5. MODELING, FEATURE SELECTION, AND RESULTS ===============================

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("plots/vimp_top30", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/vimp_selected", recursive = TRUE, showWarnings = FALSE)
dir.create("models/final", recursive = TRUE, showWarnings = FALSE)

# Initialize containers
importance_tables <- list()
results <- list()
selected_var_summary <- data.frame()
failed_models <- list()

# Loop through response variables and feature sets
for (response_var in prediction_vars) {
  cat("\n====== Modeling for:", response_var, "======\n")
  
  for (fs_name in names(feature_sets)) {
    predictors <- feature_sets[[fs_name]]
    formula_all <- as.formula(paste(response_var, "~", paste(predictors, collapse = "+")))
    
    cat("Fitting model with:", fs_name, "\n")
    
    # Fit Random Forest with error handling
    rf_model <- tryCatch({
      randomForest(formula_all, data = train_data,
                   mtry = optimal_mtry,
                   nodesize = optimal_nodesize,
                   importance = TRUE, ntree = 500)
    }, error = function(e) {
      warning(paste("RandomForest failed for", response_var, fs_name, ":", e$message))
      failed_models[[length(failed_models)+1]] <<- data.frame(Response=response_var,
                                                              FeatureSet=fs_name,
                                                              Stage="RandomForest",
                                                              Error=e$message)
      return(NULL)
    })
    
    if (is.null(rf_model)) next
    
    # Compute variable importance
    imp <- tryCatch({
      importance(rf_model, type = 1)
    }, error = function(e1) {
      warning("Type 1 importance failed, trying type 2")
      tryCatch(importance(rf_model, type = 2),
               error = function(e2) { warning("Both importance types failed"); return(NULL) })
    })
    
    if (is.null(imp)) next
    
    imp_col <- ifelse("%IncMSE" %in% colnames(imp), "%IncMSE", colnames(imp)[1])
    imp_df <- data.frame(Variable=rownames(imp), Importance=imp[, imp_col]) %>%
      arrange(desc(Importance))
    
    importance_tables[[paste(response_var, fs_name, sep="_")]] <- imp_df
    write_csv(imp_df, paste0("results/imp_raw_", response_var, "_", fs_name, ".csv"))
    
    # Plot top 30 importance variables
    top_imp_df <- imp_df %>% slice_head(n=30)
    top_plot <- ggplot(top_imp_df, aes(x=reorder(Variable, Importance), y=Importance)) +
      geom_bar(stat="identity") +
      scale_fill_viridis(option="D", direction=-1) +
      coord_flip() +
      theme_minimal() +
      labs(title=paste("Top 30 Variable Importance\n", response_var, "-", fs_name),
           x="Variable", y="% Increase in MSE")
    ggsave(paste0("plots/vimp_top30/vimp_", response_var, "_", fs_name, ".png"),
           top_plot, width=7, height=6)
    
    # Select important variables
    threshold <- 0.2 * max(imp_df$Importance, na.rm=TRUE)
    important_vars <- imp_df %>% filter(Importance >= threshold) %>% pull(Variable)
    if (length(important_vars) < 2) {
      important_vars <- imp_df %>% slice_max(order_by=Importance, n=min(5,nrow(imp_df))) %>% pull(Variable)
    }
    
    selected_var_summary <- bind_rows(selected_var_summary,
                                      data.frame(Response=response_var,
                                                 FeatureSet=fs_name,
                                                 Selected_N=length(important_vars)))
    
    selected_imp_df <- imp_df %>% filter(Variable %in% important_vars)
    write_csv(selected_imp_df, paste0("results/imp_selected_", response_var, "_", fs_name, ".csv"))
    
    # Plot selected variables
    selected_plot <- ggplot(selected_imp_df, aes(x=reorder(Variable, Importance), y=Importance)) +
      geom_bar(stat="identity") +
      scale_fill_viridis(option="C", direction=-1) +
      coord_flip() +
      theme_minimal() +
      labs(title=paste("Selected Variables\n", response_var, "-", fs_name),
           x="Variable", y="% Increase in MSE")
    ggsave(paste0("plots/vimp_selected/vimp_selected_", response_var, "_", fs_name, ".png"),
           selected_plot, width=6, height=5)
    
    # Train final model
    formula_final <- as.formula(paste(response_var, "~", paste(important_vars, collapse="+")))
    final_model <- tryCatch({
      ctrl <- trainControl(method="cv", number=5, savePredictions="final", returnResamp="all")
      train(formula_final, data=test_data, method="rf", trControl=ctrl,
            importance=TRUE, ntree=500)
    }, error=function(e){
      warning(paste("Final model failed:", response_var, fs_name))
      failed_models[[length(failed_models)+1]] <<- data.frame(Response=response_var,
                                                              FeatureSet=fs_name,
                                                              Stage="Final_Model",
                                                              Error=e$message)
      return(NULL)
    })
    
    if (is.null(final_model)) next
    saveRDS(final_model, file=paste0("models/final/model_", response_var, "_", fs_name, ".rds"))
    
    # Performance metrics
    pred_df <- final_model$pred
    pred_df <- pred_df[pred_df$mtry == final_model$bestTune$mtry, ]
    pred <- pred_df$pred
    obs <- pred_df$obs
    bias <- mean(pred - obs)
    rBias <- 100 * bias / mean(obs)
    rmse <- sqrt(mean((pred - obs)^2))
    rrmse <- 100 * rmse / mean(obs)
    r <- cor(pred, obs)
    mef <- 1 - sum((obs - pred)^2)/sum((obs - mean(obs))^2)
    
    results[[paste(response_var, fs_name, sep="_")]] <- data.frame(
      Response=response_var, Model=fs_name, R=r, RMSE=rmse,
      rRMSE=rrmse, Bias=bias, rBias=rBias, MEF=mef
    )
  }
}

# Save performance and selected variables
all_perf <- bind_rows(results)
write_csv(all_perf, "results/performance_summary.csv")
write_csv(selected_var_summary, "results/selected_variables_summary.csv")

if(length(failed_models) > 0){
  failed_df <- bind_rows(failed_models)
  write_csv(failed_df, "results/log_failed_models.csv")
}

# ============================================================================  
# End of Script  
# ============================================================================
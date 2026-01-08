#############################################################################
# This script uses WB SpaDES modules code to produce lichen biomass maps from 
# data simulated in 2018 with SpaDES biomass_core.
#############################################################################
source("https://raw.githubusercontent.com/pedrogit/rUtils/refs/heads/main/rutils.R")
# source("G:/Home/working/projects/base_withSimInit/modules/rUtils/rutils.R")

currentDir <- get_script_dir()

source("https://raw.githubusercontent.com/pedrogit/WB_HartJohnstoneForestClasses/refs/heads/main//R/WB_HartJohnstoneForestClasses.r")
# source(file.path(currentDir, "../base_withSimInit/modules/WB_HartJohnstoneForestClasses/R/WB_HartJohnstoneForestClasses.r"))

library(data.table)
library(terra)

# Read cohortData and pixelGroupMap from the data folder
currentYear <- "2031"

cohortDataFileName <- paste0("cohortData_year", currentYear, ".rds")
cohortDatafilePath <- file.path(currentDir, "data", "CanESM2_run1", cohortDataFileName)
cohortData <- readRDS(cohortDatafilePath)

pixelGroupMapFile <- paste0("pixelGroupMap_year", currentYear, ".rds")
pixelGroupMapfilePath <- file.path(currentDir, "data", "CanESM2_run1", pixelGroupMapFile)
pixelGroupMap <- readRDS(pixelGroupMapfilePath)
pixelGroupMap <- terra::rast(pixelGroupMap)

# Determine the output folder and create it
outputFolder <- paste0(currentDir, "/output/CanESM2_run1")
dir.create(outputFolder, recursive = TRUE)

# Write pixelGroupMap to disk for display in QGIS
outputRasterFileName <- paste0("pixelGroupMap_year", currentYear, ".tif")
outputRasterfilePath <- file.path(outputFolder, outputRasterFileName)
writeRaster(pixelGroupMap, outputRasterfilePath, overwrite=TRUE)

# Generate the forest classification map
WB_HJForestClassesMap <- classifyStand(
  cohortData = cohortData, 
  pixelGroupMap = pixelGroupMap,
  classificationTablePath = file.path(outputFolder, paste0("classificationTable_year", currentYear, ".csv"))
)

# Display a frequency table
freq_table(WB_HJForestClassesMap)

# Write WB_HJForestClassesMap to disk for display in QGIS
outputRasterFileName <- paste0("WB_HJForestClassesMap_year", currentYear, ".tif")
outputRasterfilePath <- file.path(outputFolder, outputRasterFileName)
writeRaster(WB_HJForestClassesMap, outputRasterfilePath, datatype = "INT1U", overwrite=TRUE)

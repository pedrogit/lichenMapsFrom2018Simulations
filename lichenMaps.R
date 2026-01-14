################################################################################
# This script uses WB SpaDES modules code to produce lichen biomass maps from 
# data simulated in 2018 with SpaDES biomass_core.
################################################################################
library(data.table)
library(terra)
library(reproducible)
library(whitebox)

source("https://raw.githubusercontent.com/pedrogit/rUtils/refs/heads/main/rutils.R")
# source(file.path(scriptDir, "../base_withSimInit/modules/rUtils/rutils.R")
scriptDir <- get_script_dir()

source("https://raw.githubusercontent.com/pedrogit/WB_HartJohnstoneForestClasses/refs/heads/main/R/WB_HartJohnstoneForestClasses.r")
# source("https://raw.githubusercontent.com/pedrogit/WB_VegBasedDrainage/refs/heads/main/R/WB_VegBasedDrainage.r")

# source(file.path(scriptDir, "../base_withSimInit/modules/WB_HartJohnstoneForestClasses/R/WB_HartJohnstoneForestClasses.r"))
source(file.path(scriptDir, "../base_withSimInit/modules/WB_VegBasedDrainage/R/WB_VegBasedDrainage.r"))

# Set the output folder and create it
outputFolder <- file.path(scriptDir, "output")
dir.create(outputFolder, recursive = TRUE)

# Set the cache folder and create it
cacheFolder <- file.path(scriptDir, "cache")
dir.create(cacheFolder, recursive = TRUE)

# Set the data folder
dataFolder <- file.path(scriptDir, "data")

# Define the set of runs to process
runs <- c("CanESM2_run1")

# Define the years to process
years <- c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")

################################################################################
# Beginning of the main loop to process each run
################################################################################
currentRun <- runs[1]

# Set the run output folder and create it
currentRunFolder <- file.path(outputFolder, currentRun)
dir.create(currentRunFolder, recursive = TRUE)

################################################################################
# Beginning of the sub loop to process each year within each run
################################################################################
currentYear <- years[1]

################################################################################
# 1 - Read pixelGroupMap from the data folder and set it as the base raster
################################################################################
pixelGroupMapFile <- paste0("pixelGroupMap_year", currentYear, ".rds")
pixelGroupMapfilePath <- file.path(dataFolder, currentRun, pixelGroupMapFile)
pixelGroupMap <- readRDS(pixelGroupMapfilePath)
pixelGroupMap <- terra::rast(pixelGroupMap)

# Write pixelGroupMap to disk for display in QGIS
outputRasterFileName <- paste0("pixelGroupMap_year", currentYear, ".tif")
outputRasterfilePath <- file.path(currentRunFolder, outputRasterFileName)
writeRaster(pixelGroupMap, outputRasterfilePath, overwrite=TRUE)

# Define this raster as the base raster for the rest of the analysis
baseRast <- rast(pixelGroupMap)
baseExtent <- ext(baseRast)
baseCRS <- crs(baseRast)
baseExtentPoly <- vect(baseExtent, crs = baseCRS)

################################################################################
# Generate initial and static data 
################################################################################
if (currentYear == "2011") {
  ##############################################################################
  # 1 - General an initial no drainage forest classification.
  ##############################################################################
  # Read cohortData and pixelGroupMap from the data folder
  cohortDataFileName <- paste0("cohortData_year", currentYear, ".rds")
  cohortDatafilePath <- file.path(dataFolder, currentRun, cohortDataFileName)
  cohortData <- readRDS(cohortDatafilePath)

  # Generate the forest classification map
  WB_HJForestClassesMap <- classifyStand(
    cohortData = cohortData, 
    pixelGroupMap = pixelGroupMap,
    classificationTablePath = file.path(currentRunFolder, paste0("classificationTable_year", currentYear, ".csv"))
  )

  # Display a frequency table
  freq_table(WB_HJForestClassesMap)

  # Write WB_HJForestClassesMap to disk for display in QGIS
  outputRasterFileName <- paste0("WB_HJForestClassesMap_year", currentYear, ".tif")
  outputRasterfilePath <- file.path(currentRunFolder, outputRasterFileName)
  writeRaster(WB_HJForestClassesMap, outputRasterfilePath, datatype = "INT1U", overwrite=TRUE)

  ################################################################################
  # 2 - Fit the WB_Drainage model and generate a drainage map.
  ################################################################################
  
  # Read and clean plot data from the module repository
  plotdataFilePath <- "https://raw.githubusercontent.com/pedrogit/WB_VegBasedDrainage/refs/heads/main/data/plotData.csv"
  drainagePlotPoints <- getAndcleanPlotData(plotdataFilePath, baseCRS)
  
  ################################################################################
  # Compute the joined area covered by the plot data AND the pixelGroupMap raster
  # (actually WB_HartJohnstoneForestClassesMap here) in order to prepInputs() 
  # covariates to this area to fit the model. Once the model is fitted, we crop 
  # the covariate back to the pixelGroupMap area
  ################################################################################
  
  # Make a buffer around the plot data
  plotPoints100KmBuffers <- aggregate(buffer(drainagePlotPoints, width = 100000))  # 1000 m buffer (if CRS is in meters)
  
  # Merge it with the base extent polygon
  plotAndPixelGroupArea <- aggregate(rbind(plotPoints100KmBuffers, baseExtentPoly))
  
  # Create a raster with the merged area so projectTo does not crop it
  # see https://github.com/PredictiveEcology/reproducible/issues/431
  plotAndPixelGroupAreaRast <- terra::extend(baseRast, plotAndPixelGroupArea)
  
  ################################################################################
  # Download and postProcess the Medium Resolution Digital Elevation Model (MRDEM)
  ################################################################################
  message("##############################################################################")   
  message("Downloading/cropping/reprojecting/resampling/masking medium resolution ")
  message("MRDEM dem (80GB) to union of studyarea and a 100km buffer around buffered plot points...")
  
  plotAndPixelGroupAreaDemPath <- file.path(cacheFolder, "plotAndPixelGroupAreaDem.tif")
  MRDEMMap <- Cache(
    prepInputs,
    url = "https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dtm.tif",
    targetFile = "mrdem-30-dtm.tif",
    destinationPath = cacheFolder,
    fun = terra::rast,
    cropTo = plotAndPixelGroupArea,
    projectTo = plotAndPixelGroupAreaRast,
    align_only = TRUE,
    maskTo = plotAndPixelGroupArea,
    method = "bilinear",
    writeTo = plotAndPixelGroupAreaDemPath,
    userTags = "MRDEMMap",
    cachePath = cacheFolder,
    overwrite = TRUE
  )
  ################################################################################
  # Generate the TWI map
  ################################################################################
  message("##############################################################################")   
  message("Generating twi from MRDEMMap...")   
  twi <- generateTWIMap(
    dem = MRDEMMap,
    dem_path = plotAndPixelGroupAreaDemPath,
    dem_filled_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_filled.tif"),
    slope_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_slope.tif"),
    flow_acc_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_flowAccum.tif"),
    final_twi_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_TWI.tif"),
    cachePath = cacheFolder
  )
  
  ################################################################################
  # Generate the Downslope Distance map
  ################################################################################
  message("##############################################################################")   
  message("Generating downslope_dist from MRDEMMap...")   
  downslope_dist <- generateDownslopeDistMap(
    dem = MRDEMMap,
    dem_path = plotAndPixelGroupAreaDemPath,
    dem_breach_filled_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_breachFilledDep.tif"),
    bf_flow_acc_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_breachFilledFlowAccum.tif"),
    streams_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_streams.tif"),
    downslope_dist_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_downslopeDist.tif"),
    cachePath = cacheFolder
  )
  
  ##############################################################################
  # Generate an aspect map from the MRDEM if it is not supplied
  ##############################################################################
  message("##############################################################################")   
  message("Generating aspect from MRDEMMap...")
  aspectPath <- file.path(cacheFolder, "plotAndPixelGroupAreaDem_aspect.tif")
  aspect <- Cache(
    cacheableWhiteboxFct,
    cacheable_input = MRDEMMap,
    fun_name = "wbt_aspect",
    dem = plotAndPixelGroupAreaDemPath,
    output = aspectPath,
    userTags = "plotAndPixelGroupAreaDem_aspect.tif",
    cachePath = cacheFolder
  )
  names(aspect) <- "aspect"

  ##############################################################################
  # Download and patch CANSIS soil maps with SoilGrids data
  ##############################################################################
  message("##############################################################################")   
  message("Download and patch CANSIS soil maps with SoilGrids data...")
  CANSISMapToProcess <- c("Clay", "Sand", "Silt", "BD") # BD is bulk_density
  equivSoilGridsMaps <- c("clay", "sand", "silt", "bdod") # bdod is bulk_density
  soilMaps <- list()
  sapply(CANSISMapToProcess, function(mapName){
    message("##############################################################################")
    varMapName <- paste0("WB_VBD_", mapName, "Map") # e.g. WB_VBD_clayMap
    soilMaps[[varMapName]] <- getAndPatchCANSISSoilMap(
      mapName = mapName,
      plotAndPixelGroupArea = plotAndPixelGroupArea,
      plotAndPixelGroupAreaRast = plotAndPixelGroupAreaRast,
      equivSoilGridsMaps = equivSoilGridsMaps,
      CANSISMapToProcess = CANSISMapToProcess,
      destinationPath = cacheFolder,
      cachePath = cacheFolder
    )
  })

  ##############################################################################
  # Download the EcoProvince map
  ##############################################################################
  message("##############################################################################")   
  message("Download the EcoProvince map...")
  ecoProvVect <- Cache(
    prepInputs,
    url = extractURL("EcoProvincesMap", sim),
    targetFile = "NA_CEC_Eco_Level3.shp",
    destinationPath = getPaths()$cache,
    projectTo = plotAndPixelGroupAreaRast,
    cropTo = plotAndPixelGroupArea,
    writeTo = file.path(getPaths()$cache, "NA_CEC_Eco_Level3_postProcessed.shp"),
    fun = terra::vect,
    userTags = "NA_CEC_Eco_Level3_postProcessed.shp",
    cachePath = cacheFolder,
    overwrite = TRUE
  )
  ecoProvVect <- ecoProvVect[, c("NA_L3NAME")]
  ecoProvVect$NA_L3NAME <- as.factor(ecoProvVect$NA_L3NAME)
  names(ecoProvVect)[names(ecoProvVect) == "NA_L3NAME"] <- "ecoprov"
  ecoprov <- rasterize(ecoProvVect, plotAndPixelGroupAreaRast, field = "ecoprov")
}
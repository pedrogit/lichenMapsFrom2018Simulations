################################################################################
# This script uses WB SpaDES modules code to produce lichen biomass maps from 
# data simulated in 2018 with SpaDES biomass_core.
################################################################################
library(data.table)
library(terra)
library(sf)
library(reproducible)
library(whitebox)
library(caret)

# Source the utility functions
source("https://raw.githubusercontent.com/pedrogit/rUtils/refs/heads/main/rutils.R")
# source(file.path(scriptDir, "../base_withSimInit/modules/rUtils/rutils.R")
scriptDir <- get_script_dir()

# Source the WB_HartJohnstoneForestClasses functions 
source("https://raw.githubusercontent.com/pedrogit/WB_HartJohnstoneForestClasses/refs/heads/main/R/WB_HartJohnstoneForestClasses.r")
# source(file.path(scriptDir, "../base_withSimInit/modules/WB_HartJohnstoneForestClasses/R/WB_HartJohnstoneForestClasses.r"))

# Source the WB_VegBasedDrainage functions 
source("https://raw.githubusercontent.com/pedrogit/WB_VegBasedDrainage/refs/heads/main/R/WB_VegBasedDrainage.r")
# source(file.path(scriptDir, "../base_withSimInit/modules/WB_VegBasedDrainage/R/WB_VegBasedDrainage.r"))

# Source the WB_NonForestedVegClasses functions 
source("https://raw.githubusercontent.com/pedrogit/WB_NonForestedVegClasses/refs/heads/main/R/WB_NonForestedVegClasses.r")
# source(file.path(scriptDir, "../base_withSimInit/modules/WB_NonForestedVegClasses/R/WB_NonForestedVegClasses.r"))

# Source the WB_LichenBiomass functions 
source("https://raw.githubusercontent.com/pedrogit/WB_LichenBiomass/refs/heads/main/R/WB_LichenBiomass.r")
# source(file.path(scriptDir, "../base_withSimInit/modules/WB_LichenBiomass/R/WB_LichenBiomass.r"))

################################################################################

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
baseRast <- ifel(is.na(pixelGroupMap), NA, 1)
baseExtent <- ext(baseRast)
baseCRS <- crs(baseRast)
baseExtentPoly <- vect(baseExtent, crs = baseCRS)

################################################################################
# Generate initial and static data that do not change over the years
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
  twi <- Cache(
    generateTWIMap(
      dem = MRDEMMap,
      dem_path = plotAndPixelGroupAreaDemPath,
      dem_filled_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_filled.tif"),
      slope_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_slope.tif"),
      flow_acc_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_flowAccum.tif"),
      final_twi_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_TWI.tif"),
      cachePath = cacheFolder
    ),
    cachePath = cacheFolder,
    userTags = "twi"
  )
  
  ################################################################################
  # Generate the Downslope Distance map
  ################################################################################
  message("##############################################################################")   
  message("Generating downslope_dist from MRDEMMap...")   
  downslope_dist <- Cache(
    generateDownslopeDistMap(
      dem = MRDEMMap,
      dem_path = plotAndPixelGroupAreaDemPath,
      dem_breach_filled_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_breachFilledDep.tif"),
      bf_flow_acc_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_breachFilledFlowAccum.tif"),
      streams_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_streams.tif"),
      downslope_dist_path = file.path(cacheFolder, "plotAndPixelGroupAreaDem_downslopeDist.tif"),
      cachePath = cacheFolder
    ),
    cachePath = cacheFolder,
    userTags = "downslope_dist"
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

  for (mapName in CANSISMapToProcess) {
    message("##############################################################################")
    message("Processing ", mapName, "...")
    assign(tolower(mapName), Cache(
        getAndPatchCANSISSoilMap(
        mapName = mapName,
        plotAndPixelGroupArea = plotAndPixelGroupArea,
        plotAndPixelGroupAreaRast = plotAndPixelGroupAreaRast,
        equivSoilGridsMaps = equivSoilGridsMaps,
        CANSISMapToProcess = CANSISMapToProcess,
        destinationPath = cacheFolder,
        userTags = tolower(mapName),
        cachePath = cacheFolder
      ),
      userTags = tolower(mapName),
      cachePath = cacheFolder
    ))
  }

  ##############################################################################
  # Download the EcoProvince map
  ##############################################################################
  message("##############################################################################")   
  message("Download the EcoProvince map...")
  ecoProvVect <- Cache(
    prepInputs,
    url ="https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip",
    targetFile = "NA_CEC_Eco_Level3.shp",
    destinationPath = cacheFolder,
    projectTo = plotAndPixelGroupAreaRast,
    cropTo = plotAndPixelGroupArea,
    writeTo = file.path(cacheFolder, "NA_CEC_Eco_Level3_postProcessed.shp"),
    fun = terra::vect,
    userTags = "NA_CEC_Eco_Level3_postProcessed.shp",
    cachePath = cacheFolder,
    overwrite = TRUE
  )
  ecoProvVect <- ecoProvVect[, c("NA_L3NAME")]
  ecoProvVect$NA_L3NAME <- as.factor(ecoProvVect$NA_L3NAME)
  names(ecoProvVect)[names(ecoProvVect) == "NA_L3NAME"] <- "ecoprov"
  ecoprov <- rasterize(ecoProvVect, plotAndPixelGroupAreaRast, field = "ecoprov")
  
  ##############################################################################
  # Fit a drainage model
  ##############################################################################
  message("##############################################################################")   
  # List the covariates from which to extract values
  # element's names are the names of the column to create in the plotPoints dataframe (e.g clay)
  # element values are maps to extract values from (e.g. twi)
  covariateMapList = list(
    twi = twi,
    downslope_dist = downslope_dist,
    aspect = aspect,
    clay = clay,
    sand = sand,
    silt = silt,
    bd = bd,
    ecoprov = ecoprov
  )
  
  # Fit the model
  drainageModel <- fit_WB_VegBasedDrainageModel(
    plotPoints = drainagePlotPoints,
    covariateMapList = covariateMapList,
    pixelDist = 2
  )
  
  # Crop covariate maps back to groupPixelMap now that the model is fitted and 
  # we don't need to extract covariate values at plot points anymore.
  # From now on, the module will, at each iteration step, use the model to predict 
  # drainage from the list of covariate maps.
  message("##############################################################################")   
  message("Now that the model is fitted using maps covering the plot data, we can crop")
  message("the covariates back to the study area (baseRast)...")
  for (mapName in names(covariateMapList)) {
    message("Cropping ", mapName, " back to study area...")   
    assign(
      mapName,
      Cache(
        postProcessTo,
        covariateMapList[[mapName]],
        cropTo = baseRast,
        userTags = mapName,
        cachePath = cacheFolder,
        overwrite = TRUE
      )
    )
  }

  ##############################################################################
  # Generate the non forested area map
  ##############################################################################
  message("##############################################################################")   
  message("Generate the non forested area map...")   
  # We do not have the equivalent of rasterToMatch which is a raster
  # representation of the study area without NA holes. We therefore have to 
  # produce it from pixelGroupMap.
  
  message("------------------------------------------------------------------------------")   
  message("Generate the study area and rasterize it...")   
  message("1 - Vectorizing the base raster (pixelGroupMap)...")
  baseRastPoly <- Cache(
    terra::as.polygons,
    x = baseRast,
    userTags = "baseRastPoly",
    cachePath = cacheFolder
  )
  
  message("2 - Removing holes...")
  # display_ring_and_holes(baseRastPoly, "terra::fillHoles")
  baseRastPolyWithoutHoles <- Cache(
    terra::fillHoles,
    x = baseRastPoly,
    userTags = "baseRastPolyWithoutHoles",
    cachePath = cacheFolder
  )
  
  message("3 - Convert to sf object...")
  # display_ring_and_holes(baseRastPolyWithoutHoles, "sf::st_as_sf")
  baseRastPolyWithoutHolesSf <- Cache(
    sf::st_as_sf,
    x = baseRastPolyWithoutHoles, 
    userTags = "baseRastPolyWithoutHolesSf",
    cachePath = cacheFolder
  )
  
  message("4 - Generate concave hull...")
  # display_ring_and_holes(baseRastPolyWithoutHolesSf, "sf::st_concave_hull")
  concaveHull <- Cache(
    sf::st_concave_hull,
    x = baseRastPolyWithoutHolesSf,
    ratio = 0.01,
    allow_holes = FALSE,
    userTags = "concaveHull",
    cachePath = cacheFolder
  )
  st_write(concaveHull, file.path(currentRunFolder, "concaveHull.shp"), delete_layer=TRUE)
  
  message("5 - Rasterizing to an equivalent raster...")
  # display_ring_and_holes(concaveHull, "terra::rasterize")
  rasterToMatch <- Cache(
    terra::rasterize,
    x = concaveHull,
    y = baseRast,
    userTags = "rasterToMatch",
    cachePath = cacheFolder
  )
  writeRaster(rasterToMatch, file.path(currentRunFolder, "rasterToMatch.tif"), overwrite=TRUE)
}


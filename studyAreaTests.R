library(terra)
library(sf)
library(nngeo)

source("https://raw.githubusercontent.com/pedrogit/rUtils/refs/heads/main/rutils.R")
# source(file.path(scriptDir, "../base_withSimInit/modules/rUtils/rutils.R")
scriptDir <- get_script_dir()

# Set the output folder and create it
outputFolder <- file.path(scriptDir, "output")
dir.create(outputFolder, recursive = TRUE)

# Set the data folder
dataFolder <- file.path(scriptDir, "data")
################################################################################
# 1 - Read pixelGroupMap from the data folder and set it as the base raster
################################################################################
pixelGroupMapFile <- paste0("pixelGroupMap_year2011.rds")
pixelGroupMapfilePath <- file.path(dataFolder, "CanESM2_run1", pixelGroupMapFile)
pixelGroupMap <- readRDS(pixelGroupMapfilePath)
pixelGroupMap <- terra::rast(pixelGroupMap)

################################################################################
basePath <- file.path(outputFolder, "temp")

baseRast <- ifel(is.na(pixelGroupMap), NA, 1)
# plot(baseRast)
# writeRaster(baseRast, file.path(basePath, "baseRast.tif"), overwrite=TRUE)

# baseRastNoHoles <- fillHoles(baseRast, nearest=TRUE)
# plot(baseRastNoHoles)
# writeRaster(baseRastNoHoles, file.path(basePath, "baseRastNoHoles.tif"), overwrite=TRUE)

baseRastPoly <- as.polygons(baseRast)
# plot(baseRastPoly)
# writeVector(baseRastPoly, file.path(basePath, "baseRastPoly.shp"), overwrite=TRUE)

baseRastPolyNoHoles <- fillHoles(baseRastPoly)
# plot(baseRastPolyNoHoles)
# writeVector(baseRastPolyNoHoles, file.path(basePath, "baseRastPolyNoHoles.shp"), overwrite=TRUE)

################################################################################
total_exterior_rings <- function(sf_obj) {
  x <- sf::st_as_sf(sf_obj)
  sum(sapply(st_geometry(x), function(g) {
    if (inherits(g, "POLYGON")) {
      1
    } else if (inherits(g, "MULTIPOLYGON")) {
      length(g)  # one exterior ring per polygon part
    } else {
      0
    }
  }))
}

total_holes <- function(sf_obj) {
  x <- sf::st_as_sf(sf_obj)
  sum(sapply(st_geometry(x), function(g) {
    if (inherits(g, "POLYGON")) {
      max(0, length(g) - 1)  # interior rings
    } else if (inherits(g, "MULTIPOLYGON")) {
      sum(sapply(g, function(poly) max(0, length(poly) - 1)))
    } else {
      0
    }
  }))
}

display_ring_and_holes <- function(x, fct = ""){
  message(fct, "() with ", nrow(x), " (multi)polygons, ", total_exterior_rings(x), " exterior rings, ", total_holes(x), " holes...")
}

################################################################################
if (FALSE) {
# ConcaveHull - Default params = convexHull
concaveHullRatio <- hull(baseRastPolyNoHoles, type="concave_ratio", allowHoles=FALSE)
# plot(concaveHullRatio)
# writeVector(concaveHullRatio, file.path(basePath, "concaveHullRatio.shp"), overwrite=TRUE)

# ConcaveHull - params = 0.5
concaveHullRatio_0_5 <- hull(baseRastPolyNoHoles, type="concave_ratio", param=0.5, allowHoles=FALSE)
# plot(concaveHullRatio_0_5)
# writeVector(concaveHullRatio_0_5, file.path(basePath, "concaveHullRatio_0_5.shp"), overwrite=TRUE)

# ConcaveHull - params = 0.2
concaveHullRatio_0_2 <- hull(baseRastPolyNoHoles, type="concave_ratio", param=0.2, allowHoles=FALSE)
# plot(concaveHullRatio_0_2)
# writeVector(concaveHullRatio_0_2, file.path(basePath, "concaveHullRatio_0_2.shp"), overwrite=TRUE)

# ConcaveHull - params = 0.05 - Best
concaveHullRatio_0_05 <- hull(baseRastPolyNoHoles, type="concave_ratio", param=0.05, allowHoles=FALSE)
# plot(concaveHullRatio_0_05)
# writeVector(concaveHullRatio_0_05, file.path(basePath, "concaveHullRatio_0_05.shp"), overwrite=TRUE)

################################################################################
# ConcaveHull - Default params = convexHull
concaveHullLength <- hull(baseRastPolyNoHoles, type="concave_length", allowHoles=FALSE)
# plot(concaveHullLength)
# writeVector(concaveHullLength, file.path(basePath, "concaveHullLength.shp"), overwrite=TRUE)

# ConcaveHull - params = 0.5
concaveHullLength_0_5 <- hull(baseRastPolyNoHoles, type="concave_length", param=0.5, allowHoles=FALSE)
# plot(concaveHullLength_0_5)
# writeVector(concaveHullLength_0_5, file.path(basePath, "concaveHullLength_0_5.shp"), overwrite=TRUE)

# ConcaveHull - params = 0.2
concaveHullLength_0_2 <- hull(baseRastPolyNoHoles, type="concave_length", param=0.2, allowHoles=FALSE)
# plot(concaveHullLength_0_2)
# writeVector(concaveHullLength_0_2, file.path(basePath, "concaveHullLength_0_2.shp"), overwrite=TRUE)

# ConcaveHull - params = 0.05
concaveHullLength_0_05 <- hull(baseRastPolyNoHoles, type="concave_length", param=0.05, allowHoles=FALSE)
# plot(concaveHullLength_0_05)
# writeVector(concaveHullLength_0_05, file.path(basePath, "concaveHullLength_0_05.shp"), overwrite=TRUE)

# ConcaveHull - params = 0.02
concaveHullLength_0_02 <- hull(baseRastPolyNoHoles, type="concave_length", param=0.02, allowHoles=FALSE)
# plot(concaveHullLength_0_02)
# writeVector(concaveHullLength_0_02, file.path(basePath, "concaveHullLength_0_02.shp"), overwrite=TRUE)

system.time({
# ConcaveHull - params = 0.01
concaveHullLength_0_01 <- hull(baseRastPolyNoHoles, type="concave_length", param=0.01, allowHoles=FALSE)
})
# plot(concaveHullLength_0_01)
# writeVector(concaveHullLength_0_01, file.path(basePath, "concaveHullLength_0_01.shp"), overwrite=TRUE)



################################################################################
baseRastPolySf <- sf::st_as_sf(baseRastPoly)
display_ring_and_holes(baseRastPolySf)
}

baseRastPolyNoHolesSf <- sf::st_as_sf(baseRastPolyNoHoles)
display_ring_and_holes(baseRastPolyNoHolesSf)

if (FALSE) {
# SF ConcaveHull - Default params = convexHull
concaveHullSf_0_5 <- sf::st_concave_hull(baseRastPolyNoHolesSf, ratio=0.5, allow_holes=FALSE)
# plot(concaveHullSf_0_5)
st_write(concaveHullSf_0_5, file.path(basePath, "concaveHullSf_0_5.shp"), delete_layer=TRUE)

# ConcaveHull - params = 0.2
concaveHullSf_0_2 <- sf::st_concave_hull(baseRastPolyNoHolesSf, ratio=0.2, allow_holes=FALSE)
# plot(concaveHullSf_0_2)
st_write(concaveHullSf_0_2, file.path(basePath, "concaveHullSf_0_2.shp"), delete_layer=TRUE)

# ConcaveHull - params = 0.05
concaveHullSf_0_05 <- sf::st_concave_hull(baseRastPolyNoHolesSf, ratio=0.05, allow_holes=FALSE)
# plot(concaveHullSf_0_05)
st_write(concaveHullSf_0_05, file.path(basePath, "concaveHullSf_0_05.shp"), delete_layer=TRUE)

# ConcaveHull - params = 0.02
concaveHullSf_0_02 <- sf::st_concave_hull(baseRastPolyNoHolesSf, ratio=0.02, allow_holes=FALSE)
# plot(concaveHullSf_0_02)
st_write(concaveHullSf_0_02, file.path(basePath, "concaveHullSf_0_02.shp"), delete_layer=TRUE)

system.time({
  # ConcaveHull - params = 0.01 - Best 177 seconds
  concaveHullSfWithHoles_0_01 <- sf::st_concave_hull(baseRastPolySf, ratio=0.01, allow_holes=FALSE)
})
st_write(concaveHullSfWithHoles_0_01, file.path(basePath, "concaveHullSfWithHoles_0_01.shp"), delete_layer=TRUE)
}

system.time({
  # ConcaveHull - params = 0.01 - Best - 53 seconds
  concaveHullSfWithoutHoles_0_01 <- sf::st_concave_hull(baseRastPolyNoHolesSf, ratio=0.01, allow_holes=FALSE)
})
st_write(concaveHullSfWithoutHoles_0_01, file.path(basePath, "concaveHullSfWithoutHoles_0_01.shp"), delete_layer=TRUE)


# plot(concaveHullSf_0_01)
st_write(concaveHullSf_0_01, file.path(basePath, "concaveHullSf_0_01.shp"), delete_layer=TRUE)

################################################################################
if (FALSE) {
addRemoveBuffer <- function(x, width = NULL){
  if (is.null(width)){
    message("width is NULL. Returning x...")
    return(x)
  }
  if (width == 0){
    message("width is 0. Returning x...")
    return(x)
  }
  x <- sf::st_as_sf(x)
  xcrs <- sf::st_crs(x)
  
  # x <- sf::st_geometry(x)
  
  display_ring_and_holes(x, "st_exterior_ring")
  # message("st_remove_holes()...")
  # x <- nngeo::st_remove_holes(x)
  # x <- sf::st_exterior_ring(x)
  x <- sf::st_cast(sf::st_exterior_ring(x), "POLYGON")
  
  display_ring_and_holes(x, "st_buffer")
  x <- sf::st_buffer(x, width)
  
  display_ring_and_holes(x, "st_union")
  x <- sf::st_union(x)
  # message("st_remove_holes()...")
  # x <- nngeo::st_remove_holes(x)
  
  display_ring_and_holes(x, "-st_buffer")
  x <- sf::st_buffer(x, -width)
  
  display_ring_and_holes(x, "st_exterior_ring")
  x <- sf::st_cast(sf::st_exterior_ring(x), "POLYGON")
  
  sf::st_crs(x) <- xcrs
  return(x)
}

# ChatGPT version (takes forever)
addRemoveBufferGPT <- function(sf_obj, buffer_dist = 10) {
  sf_obj <- sf::st_as_sf(sf_obj)
  
  # 1. Remove empty geometries (helps speed and avoids errors)
  sf_obj <- sf_obj[!st_is_empty(sf_obj), ]
  if (nrow(sf_obj) == 0) return(NULL)
  
  # 2. Apply positive buffer
  display_ring_and_holes(sf_obj, "st_buffer")
  buffered <- st_buffer(sf_obj, buffer_dist)
  
  # 3. Union all geometries into a single geometry
  display_ring_and_holes(buffered, "st_union")
  unioned <- st_union(buffered)
  
  # 4. Negative buffer to shrink back, removing small holes automatically
  display_ring_and_holes(unioned, "-st_buffer")
  hull <- st_buffer(unioned, -buffer_dist)
  
  # 5. Ensure no holes: st_make_valid + st_cast to POLYGON
  display_ring_and_holes(hull, "st_make_valid")
  hull <- st_make_valid(hull)
  
  display_ring_and_holes(hull, "st_cast")
  hull <- st_cast(hull, "POLYGON")
  
  # Return as sf object
  st_sf(geometry = hull)
}

system.time({ # 330 seconds
  y30 <- addRemoveBuffer(y, res(baseRast)[1] * 30)
})

system.time({
  y30GPT <- addRemoveBufferGPT(sf::st_as_sf(y), res(baseRast)[1] * 30)
})

# plot(y30)
st_write(sf::st_as_sf(y30), file.path(basePath, "y30.shp"), delete_layer=TRUE)
}

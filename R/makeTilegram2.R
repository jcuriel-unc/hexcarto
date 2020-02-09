#' Create a hexagonal map 
#' 
#' This function takes regular shapefiles with a population field and firt converts the maps to a hexagonal grid, with data preserved. Note that 
#' this function is heavily based off of the original, though is an improvement.   
#' @param shp The spatial object to be read in, providing the coordinates and the dbf 
#' @param cellsize The optional argument to change cellsize, Default is null, and probably best to leave empty 
#' @return The transformed spatial object and data frame, though polygons now in hexagon format. 
#' @export
#' @examples 
#' mi_shp <- readOGR(eguia_path, "mi_county_shp")
#' test_hexagon <- hexcarto2b(mi_shp)
#' 
makeTilegram2 <- function (sp, cellsize = NULL){
  list.of.packages <- c("rgdal","rgeos","BAMMtools","GISTools","moments","plotrix","sp","clue","devtools","roxygen2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  if(!"makeTilegram" %in% installed.packages()[,"Package"])
    devtools::install_git("https://gitlab.com/lajh87/makeTilegram")
  #####Here is the make tile gram path 
  sp <- sp::spTransform(sp, sp::CRS("+init=EPSG:32663"))
  tiles <- hex_tiles(sp)
  tiles <- tiles[sp, ]
  pts <- rgeos::gCentroid(sp, byid = T)
  pts <- sp::SpatialPointsDataFrame(pts, data.frame(sp@data, pt_id = row.names(pts),
                                                    stringsAsFactors = F))##still works here; data frame present
  tileCentroids <- rgeos::gCentroid(tiles, T)
  tileCentroids <- sp::SpatialPointsDataFrame(tileCentroids, 
                                              data.frame(id = row.names(tileCentroids), stringsAsFactors = F))
  distance <- rgeos::gDistance(tileCentroids, pts, byid = T)
  tile_pref <- t(apply(distance, 1, function(x) rank(x, ties.method = "random")))
  solved <- clue::solve_LSAP(tile_pref, maximum = FALSE)
  solved_cols <- as.numeric(solved)
  newDat <- data.frame(tile_region = row.names(tile_pref), 
                       id = as.numeric(colnames(tile_pref)[solved_cols]), stringsAsFactors = F)
  newDat <- cbind(newDat, pts@data)
  newTiles <- tiles
  newTiles@data <- plyr::join(newTiles@data, newDat, by = "id")
  newTiles <- newTiles[!is.na(newTiles$tile_region), ]
  return(newTiles)
}
#' Create a hexagonal cartogram map 
#' 
#' This function takes shapefiles with a population field and first converts the maps to a hexagonal grid, and then weights areas with 
#' higher populations to appear larger 
#' @param shp The spatial object to be read in, providing the coordinates and the dbf 
#' @param pop_field The name of the population field within the data set that will later be transformed to provide the weights 
#' @param choro_field The name of the optional field with the values that will later go onto provide the values for the choropleth map
#' @param colval The name of the optional field that if true, will provide colors by polygon for the map. Note that if this field and the 
#' choro/jenks fields are empty, then the map will be made blue. The colval field should take the form of a value that can be recognized as 
#' as color (i.e. blue). Can also be manipulated before hand by the user so that it reflects some non-jenks categorization of values.
#' @param jenks The True/False field that if true, will produce a grayscale natural breaks choropleth, as provided with the choro_field.
#' @param label_field The name of the field with polygon labels, which if filled, will overlay the polygons with names. 
#' @param quant_carto_breaks The True/False field as to how the user wants the population weights. If true, then the weights will be transformed
#' such that the weights will be broken by the values quantile (i.e. pops in fifth pct will have size weight of 0.05). False, and the values 
#' will be the normalized population values, which might result in oddities if the population is not normally distributed. 
#' @return The dataframe and plot of the map will be provided, with the following fields: 
#'     \itemize{
#'     \item xcor = The x coordinate for the centroid of the hexagon
#'     \item ycor = The y coordinate for the centroid of the hexagon
#'     \item pop_field = The populations of the polygonal units provided  
#'     \item color_field = The color values for the polygon 
#'     \item norm_pop = The normalized population values on a 0 - 1 scale
#'     \item label = If user provided, the labels of the polygons 
#'     \item norm_pop2 = The quantile apportioned weights for the map   
#'     
#' }
#' @export
#' @examples 
#' mi_shp <- readOGR(eguia_path, "mi_county_shp")
#' test_hexagon <- hexcarto2b(mi_shp,pop_field="VAP2010",choro_field="gop_vote_s",jenks = TRUE,quant_carto_breaks = TRUE)

#' 

hexcarto2b <- function(shp,pop_field,choro_field,jenks=c(TRUE,FALSE), colval,label_field, quant_carto_breaks=c(TRUE,FALSE)){
  list.of.packages <- c("rgdal","rgeos","BAMMtools","GISTools","moments","plotrix","sp","clue","devtools","roxygen2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
    library(rlang)
  for(pkg in (list.of.packages)){
    eval(bquote(library(.(pkg))))
  }
  if(!"makeTilegram" %in% installed.packages()[,"Package"])
    devtools::install_git("https://gitlab.com/lajh87/makeTilegram")
  library(makeTilegram)
  #####Here is the make tile gram path 
  sp <- sp::spTransform(shp, sp::CRS("+proj=longlat +datum=WGS84 +EPSG:4326"))
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
  ###this is the end of tilegram code 
  #return(newTiles)
  shp_cents <- gCentroid(newTiles, byid = TRUE)
  cents <- SpatialPointsDataFrame(coords=shp_cents, data=newTiles@data, 
                                  proj4string=CRS("+proj=longlat +datum=WGS84 +EPSG:4326"))
  cents$xcor <- cents@coords[,1]
  cents$ycor <- cents@coords[,2]
  pop_position <- match(pop_field, names(cents))
  cents$pop_field <- cents@data[,pop_position]
  if(is.character(cents$pop_field)==TRUE){
    cents$pop_field <- as.numeric(cents$pop_field)
  }else if(is.factor(cents$pop_field)==TRUE){
    cents$pop_field <- as.numeric(as.character(cents$pop_field))
    warning("Pop. field was factor and converted to numeric.")
  }
  temp_skew <- skewness(cents$pop_field,na.rm=T)
  temp_kurtosis <- kurtosis(cents$pop_field,na.rm=T)
  print(paste0("Pop. field skewness is: ", temp_skew))
  print(paste0("Pop. field kurtosis is: ", temp_kurtosis))
  if(missing(choro_field)==FALSE){
    choro_pos <- match(choro_field, names(cents) )
    cents$choro_field <- as.numeric(cents@data[,choro_pos])
  }else{
    print("No field for choropleth values provided.")
  }
  try(col_pos <- match(colval, names(cents)))
  try(cents$color_field <- cents@data[,col_pos])
  
  if(missing(jenks)==TRUE | jenks==FALSE){
    print("Skipping the jenks creation process.")
  }else if(missing(jenks)==FALSE | jenks==TRUE ){
    if(missing(choro_field)==T){
      stop("No choropleth field provided for purposes of creating jenks breaks.")
    }
    vals_jenks <- getJenksBreaks(cents$choro_field, 5)
    grays <- grey.colors(6)
    cents$color_field <- NA
    cents$color_field[cents$choro_field < vals_jenks[1] ] <- grays[6]
    cents$color_field[cents$choro_field >= vals_jenks[1] & cents$choro_field < vals_jenks[2]  ] <- grays[5]
    cents$color_field[cents$choro_field >= vals_jenks[2] & cents$choro_field < vals_jenks[3]  ] <- grays[4]
    cents$color_field[cents$choro_field >= vals_jenks[3] & cents$choro_field < vals_jenks[4]  ] <- grays[3]
    cents$color_field[cents$choro_field >= vals_jenks[4] & cents$choro_field < vals_jenks[5]  ] <- grays[2]
    cents$color_field[cents$choro_field >= vals_jenks[5]  ] <- grays[1]
  }else if(missing(jenks)==TRUE | jenks==FALSE | missing(col_val)==FALSE){
    if(missing(col_val)==T){
      print("Color field missing, so will make map blue.")
      cents$color_field <- "blue"
    }
  }
  data <- data.frame(cents$xcor,cents$ycor,cents$pop_field, cents$color_field)
  colnames(data) <- c("xcor","ycor","pop_field","color_field")
  if(missing(choro_field)==FALSE){
    data <- cbind(data, cents$choro_field)
    colnames(data)[5] <- "choro_field"
  }else{
    data$choro_field <- 0
  }
  norm1 <- as.data.frame(apply(as.data.frame(data$pop_field), 2, function(x) (x - min(x))/(max(x)-min(x))))
  print("Pop field normalized.")
  data <- cbind(data, norm1)
  colnames(data)[6] <- "norm_pop"
  if(missing(label_field)==FALSE){
    label_pos <- match(label_field, names(cents))
    data <- cbind(data, cents@data[,label_pos])
    colnames(data)[7] <- "label"
  }
  data$norm_pop[data$norm_pop==0] <- 0.0001
  
  pop_quants <- quantile(data$norm_pop, seq(0,1,by=0.05))
  data$norm_pop2 <- 0.05
  data$norm_pop2[data$norm_pop >= pop_quants[2] & data$norm_pop < pop_quants[5] ] <- 0.1
  data$norm_pop2[data$norm_pop >= pop_quants[5] & data$norm_pop < pop_quants[7] ] <- 0.15
  data$norm_pop2[data$norm_pop >= pop_quants[7] & data$norm_pop < pop_quants[9] ] <- 0.2
  data$norm_pop2[data$norm_pop >= pop_quants[9] & data$norm_pop < pop_quants[11] ] <- 0.25
  data$norm_pop2[data$norm_pop >= pop_quants[11] & data$norm_pop < pop_quants[13] ] <- 0.3
  data$norm_pop2[data$norm_pop >= pop_quants[13] & data$norm_pop < pop_quants[15] ] <- 0.35
  data$norm_pop2[data$norm_pop >= pop_quants[15] & data$norm_pop < pop_quants[17] ] <- 0.4
  data$norm_pop2[data$norm_pop >= pop_quants[17] & data$norm_pop < pop_quants[19] ] <- 0.45
  data$norm_pop2[data$norm_pop >= pop_quants[19] ] <- 0.5
  ###finally plotting here 
  plot(min(data$xcor),min(data$ycor), type="n", xlim=c(min(data$xcor),max(data$xcor)+.5), ylim=c(min(data$ycor),max(data$ycor)+.5),
       frame.plot=F, xaxt="n", yaxt="n", xlab="", ylab="")
  if(missing(quant_carto_breaks)==TRUE){
    for(i in 1:nrow(data)){
      hexagon(data$xcor[i],data$ycor[i],col=data$color_field[i], unitcell=data$norm_pop2[i],border="white")
    }
  }else if(quant_carto_breaks==T){
    for(i in 1:nrow(data)){
      hexagon(data$xcor[i],data$ycor[i],col=data$color_field[i], unitcell=data$norm_pop2[i],border="white")
    }
  }else if(quant_carto_breaks==F){
    for(i in 1:nrow(data)){
      hexagon(data$xcor[i],data$ycor[i],col=data$color_field[i], unitcell=data$norm_pop[i],border="white")
    }
  }
  #apply(data, 1, function(zone) hexagon(zone[1],zone[2],col=zone[4], unitcell=zone[6],border="white"))
  if(missing(label_field)==FALSE){
    text(cents$xcor+0.1,cents$ycor+0.1,labels=cents$NAME, cex=0.4)
  }
  return(data)
  #note: while it seems that everything mostly works, it is definitely the case that the unit cell is too small. Therefore, I'll create
  # 5 sizes so as to prevent odd maps 
}

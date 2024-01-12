## ----initialize, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
### needed libraries
library(cmocean)
library(mapdata)
library(plotdap)
library(rerddap)
library(rerddapXtracto)


## ----install, eval = FALSE----------------------------------------------------
#  install.packages("ncdf4", dependencies = TRUE)
#  install.packages("parsedate", dependencies = TRUE)
#  install.packages("plotdap", dependencies = TRUE)
#  install.packages("rerddap", dependencies = TRUE)
#  install.packages("sp", dependencies = TRUE)

## ----install_package, eval = FALSE--------------------------------------------
#  install.packages("rerddapXtracto", dependencies = TRUE)

## ----installGit, eval = FALSE-------------------------------------------------
#  install.packages("devtools")
#  remotes::install_github("rmendels/rerddapXtracto", subdir = 'development')

## ----eval = FALSE-------------------------------------------------------------
#  library("rerddapXtracto")

## ----eval = FALSE-------------------------------------------------------------
#  library("gganimate")
#  library("ggplot2")
#  library("plotdap")

## ----info, echo = TRUE, eval = FALSE------------------------------------------
#  require("rerddap")
#  ## base URL does not need to given because it is the default one
#  dataInfo <- info('erdMBchla1day')
#  dataInfo

## ----dataInfo, echo = TRUE, eval = FALSE--------------------------------------
#  <ERDDAP info> erdMBchla1day
#   Base URL: https://upwell.pfeg.noaa.gov/erddap/
#   Dimensions (range):
#       time: (2006-01-01T12:00:00Z, 2020-10-27T12:00:00Z)
#       altitude: (0.0, 0.0)
#       latitude: (-45.0, 65.0)
#       longitude: (120.0, 320.0)
#   Variables:
#       chlorophyll:
#           Units: mg m-3
#  

## ----getMarlinChl, eval = FALSE, echo = TRUE----------------------------------
#  require("rerddap")
#  require("rerddapXtracto")
#  
#  # First we will copy the Marlintag38606 data into a variable
#  # called tagData  so that subsequent code will be more generic.
#  
#  tagData <- Marlintag38606
#  xpos <- tagData$lon
#  ypos <- tagData$lat
#  tpos <- tagData$date
#  zpos <- rep(0., length(xpos))
#  swchlInfo <- rerddap::info('erdSWchla8day')
#  swchl1 <- rxtracto(swchlInfo, parameter = 'chlorophyll',
#                     xcoord = xpos, ycoord = ypos, tcoord = tpos, zcoord = zpos,
#                     xlen = .2, ylen = .2, progress_bar = TRUE)

## ----meantrackPlot, eval = FALSE, echo = TRUE---------------------------------
#  require("ggplot2")
#  require("plotdap")
#  
#  myPlot <- plotTrack(swchl1, xpos, ypos, tpos, plotColor = 'algae')
#  myPlot

## ----animateTrack, echo = TRUE, eval = FALSE----------------------------------
#  myPlot <- plotTrack(swchl1, xpos, ypos, tpos, plotColor = 'algae',
#                      animate = TRUE, cumulative = TRUE)
#  

## ----tidy_track, echo = TRUE, eval = FALSE------------------------------------
#  
#  swchl1_tidy <- as.data.frame(swchl1)
#  

## ----topotag, eval = FALSE, echo = TRUE---------------------------------------
#  require("ggplot2")
#  require("plotdap")
#  require("rerddap")
#  require("rerddapXtracto")
#  ylim <- c(15, 30)
#  xlim <- c(-160, -105)
#  topoInfo <- rerddap::info('etopo360')
#  topo <- rxtracto(topoInfo, parameter = 'altitude', xcoord = xpos, ycoord = ypos, xlen = .1, ylen = .1)
#  myFunc = function(x) -x
#  topoPlot <- plotTrack(topo, xpos, ypos, NA, plotColor = 'dense',
#                        name = 'Depth', myFunc = myFunc)
#  topoPlot

## ----extract3D, echo = TRUE, eval = FALSE-------------------------------------
#  require("rerddap")
#  urlBase <- "https://erddap.marine.ie/erddap/"
#  parameter <- "Sea_water_temperature"
#  dataInfo <- rerddap::info("IMI_CONN_3D", url = urlBase)
#  #get the actual last 3 times,  and extract from data frame
#  dataInfo1 <- read.csv("https://erddap.marine.ie/erddap/griddap/IMI_CONN_3D.csv0?time[last-2:1:last]", stringsAsFactors = FALSE, header = FALSE, row.names = NULL)
#  sstTimes <- dataInfo1[[1]]
#  sstLats <- c(53.505758092414446, 53.509303546859805, 53.51284900130517)
#  sstLons <- c(-10.25975390624996, -10.247847656249961, -10.23594140624996)
#  sstDepths <- c(2, 6, 10)
#  sstTrack <- rxtracto(dataInfo, parameter = parameter, xcoord = sstLons, ycoord = sstLats, tcoord = sstTimes, zcoord = sstDepths, xlen = .05, ylen = .05, zlen = 0., zName = 'altitude')
#  str(sstTrack)

## ----extract3D_struct, echo = TRUE, eval = FALSE------------------------------
#  List of 13
#   $ mean Sea_water_temperature  : num [1:3] 11.7 11.7 11.6
#   $ stdev Sea_water_temperature : num [1:3] 0.0407 0.0555 0.0859
#   $ n                           : int [1:3] 493 491 484
#   $ satellite date              : chr [1:3] "2020-10-31T22:00:00Z" "2020-10-31T23:00:00Z" "2020-11-01T00:00:00Z"
#   $ requested lon min           : num [1:3] -10.3 -10.3 -10.3
#   $ requested lon max           : num [1:3] -10.2 -10.2 -10.2
#   $ requested lat min           : num [1:3] 53.5 53.5 53.5
#   $ requested lat max           : num [1:3] 53.5 53.5 53.5
#   $ requested z min             : num [1:3] 2 6 10
#   $ requested z max             : num [1:3] 2 6 10
#   $ requested date              : chr [1:3] "2020-10-31T22:00:00Z" "2020-10-31T23:00:00Z" "2020-11-01T00:00:00Z"
#   $ median Sea_water_temperature: num [1:3] 11.7 11.7 11.6
#   $ mad Sea_water_temperature   : num [1:3] 0.0208 0.0361 0.0986
#   - attr(*, "row.names")= chr [1:3] "1" "2" "3"
#   - attr(*, "class")= chr [1:2] "list" "rxtractoTrack"
#  

## ----dateline_track, echo = TRUE,  eval = FALSE-------------------------------
#  dataInfo <- rerddap::info('jplMURSST41mday')
#  parameter <- 'sst'
#  xcoord <- c(179.7, 179.8, 179.9, 180., 180.1, 180.2, 180.3, 180.4)
#  ycoord <- c(40, 40, 40, 40, 40, 40, 40, 40)
#  tcoord <- c('2018-03-16', '2018-03-16', '2018-03-16','2018-03-16','2018-03-16',
#              '2018-03-16','2018-03-16','2018-03-16')
#  xlen <- .05
#  ylen <- .05
#  extract <- rxtracto(dataInfo, parameter = parameter, xcoord = xcoord,
#                      ycoord = ycoord, tcoord = tcoord,
#                      xlen = xlen, ylen = ylen)
#  str(extract)

## ----dateline_track_struct, echo = TRUE,  eval = FALSE------------------------
#  List of 13
#   $ mean sst         : num [1:8] 11.1 11.1 11.1 11.2 11.1 ...
#   $ stdev sst        : num [1:8] 0.01192 0.00602 0.01025 0.00876 0.01446 ...
#   $ n                : int [1:8] 30 30 35 25 30 30 30 35
#   $ satellite date   : chr [1:8] "2018-03-16T00:00:00Z" "2018-03-16T00:00:00Z" "2018-03-16T00:00:00Z" "2018-03-16T00:00:00Z" ...
#   $ requested lon min: num [1:8] 180 180 180 180 180 ...
#   $ requested lon max: num [1:8] 180 180 180 180 180 ...
#   $ requested lat min: num [1:8] 40 40 40 40 40 ...
#   $ requested lat max: num [1:8] 40 40 40 40 40 ...
#   $ requested z min  : logi [1:8] NA NA NA NA NA NA ...
#   $ requested z max  : logi [1:8] NA NA NA NA NA NA ...
#   $ requested date   : chr [1:8] "2018-03-16" "2018-03-16" "2018-03-16" "2018-03-16" ...
#   $ median sst       : num [1:8] 11.1 11.1 11.1 11.2 11.1 ...
#   $ mad sst          : num [1:8] 0.01416 0.0052 0.01149 0.00887 0.01744 ...
#   - attr(*, "row.names")= chr [1:8] "1" "2" "3" "4" ...
#   - attr(*, "class")= chr [1:2] "list" "rxtractoTrack"
#  

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  interp <- c('interpolation method', 'number of neighbors')

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  interp <- c('Mean',  '16')

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  my_interp <- c('Mean',  '16')
#  #
#  #  use the coastwatch ERDDAP as it is a more recent version
#  #
#  swchlInfo <- rerddap::info('erdSWchla8day',
#               url = 'https://coastwatch.pfeg.noaa.gov/erddap/')
#  swchl1 <- rxtracto(swchlInfo, parameter = 'chlorophyll',
#                     xcoord = xpos, ycoord = ypos, tcoord = tpos, zcoord = zpos,
#                     interp = my_interp, progress_bar = TRUE)
#  

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  PB_Argos_subset <- PB_Argos[1:100, ]
#  # datetime is not in format for ERDDAP
#  PB_Argos_subset$DateTime <- lubridate::as_datetime(PB_Argos_subset$DateTime)
#  head(PB_Argos_subset, 4)
#  
#  # A tibble: 4 × 4
#    DateTime            QualClass   Lat   Lon
#    <dttm>              <chr>     <dbl> <dbl>
#  1 2009-04-20 17:01:39 B          70.4 -132.
#  2 2009-04-20 17:23:00 A          71.0 -131.
#  3 2009-04-20 18:12:15 A          71.0 -131.
#  4 2009-04-20 20:43:17 A          70.9 -131.
#  

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  url=  'https://polarwatch.noaa.gov/erddap/'
#  datasetid = 'nsidcG02202v4nh1day'
#  dataInfo <- rerddap::info(datasetid, url)
#  proj_crs_code_index <- which(dataInfo$alldata$NC_GLOBAL$attribute_name == "proj_crs_code" )
#  proj_crs_code <- dataInfo$alldata$NC_GLOBAL$value[proj_crs_code_index]
#  proj_crs_code
#  [1] "EPSG:3411"

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  temp_df <- data.frame(Lat = PB_Argos_subset$Lat, Lon = PB_Argos_subset$Lon)
#  # transform PB_Argos_subset to sf object
#  # EPSG:4326 is basic lat-lon coordinates
#  temp_df <- sf::st_as_sf(temp_df, coords = c("Lon", "Lat"), crs = 4326)
#  # project data
#  temp_df <- sf::st_transform(temp_df, crs = 3411)
#  # get projection coordinates
#  coordinates <- sf::st_coordinates(temp_df)

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  xcoord <- coordinates[, 1]
#  ycoord <- coordinates[, 2]
#  # R datetimes are passed as number,  require ISO character string
#  tcoord <- as.character(PB_Argos_subset$DateTime)
#  parameter <- 'cdr_seaice_conc'
#  extract <- rxtracto(dataInfo,
#                      xName="xgrid",
#                      yName="ygrid",
#                      tName="time",
#                      parameter=parameter,
#                      xcoord = xcoord,
#                      ycoord = ycoord,
#                      tcoord = tcoord
#                      )
#  str(extract)
#  List of 13
#   $ mean cdr_seaice_conc  : num [1:100] 1 1 1 1 1 1 1 1 1 1 ...
#   $ stdev cdr_seaice_conc : num [1:100] NA NA NA NA NA NA NA NA NA NA ...
#   $ n                     : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
#   $ satellite date        : chr [1:100] "2009-04-21T00:00:00Z" "2009-04-21T00:00:00Z" "2009-04-21T00:00:00Z" "2009-04-21T00:00:00Z" ...
#   $ requested x min       : num [1:100] -2143958 -2077605 -2077124 -2082576 -2036149 ...
#   $ requested x max       : num [1:100] -2143958 -2077605 -2077124 -2082576 -2036149 ...
#   $ requested y min       : num [1:100] -118477 -131404 -131955 -130368 -120720 ...
#   $ requested y max       : num [1:100] -118477 -131404 -131955 -130368 -120720 ...
#   $ requested z min       : logi [1:100] NA NA NA NA NA NA ...
#   $ requested z max       : logi [1:100] NA NA NA NA NA NA ...
#   $ requested date        : chr [1:100] "2009-04-20 17:01:39" "2009-04-20 17:23:00" "2009-04-20 18:12:15" "2009-04-20 20:43:17" ...
#   $ median cdr_seaice_conc: num [1:100] 1 1 1 1 1 1 1 1 1 1 ...
#   $ mad cdr_seaice_conc   : num [1:100] 0 0 0 0 0 0 0 0 0 0 ...
#   - attr(*, "row.names")= chr [1:100] "1" "2" "3" "4" ...
#   - attr(*, "class")= chr [1:2] "list" "rxtractoTrack"
#  

## ----echo = TRUE,  eval = FALSE-----------------------------------------------
#  xgrid <- 'requested x min'
#  ygrid <- 'requested y min'
#  temp_df <- data.frame(xgrid = extract[[xgrid]], ygrid = extract[[ygrid]])
#  temp_df <- sf::st_as_sf(temp_df, coords = c("xgrid", "ygrid"), crs = 3411)
#  temp_df <- sf::st_transform(temp_df, crs = 4326)
#  coordinates <- sf::st_coordinates(temp_df)
#  

## ----VIIRSchla, echo = TRUE, eval = FALSE-------------------------------------
#  require("rerddap")
#  require("rerddapXtracto")
#  
#  xpos <- c(-125, -120)
#  ypos <- c(39, 36)
#  tpos <- c("last", "last")
#  tpos <- c("2017-04-15", "2017-04-15")
#  VIIRSInfo <- rerddap::info('erdVH3chlamday')
#  VIIRS <- rxtracto_3D(VIIRSInfo, parameter = 'chla', xcoord = xpos, ycoord = ypos, tcoord = tpos)

## ----VIIRSLogPlot, echo = TRUE, eval = FALSE----------------------------------
#  require("ggplot2")
#  require("plotdap")
#  myFunc <- function(x) log(x)
#  chlalogPlot <- plotBBox(VIIRS, plotColor = 'algae', myFunc = myFunc)
#  chlalogPlot

## ----viirs_tidy, echo = TRUE,  eval = FALSE-----------------------------------
#  VIIRS_tidy <- tidy_grid(VIIRS)

## ----viirs_tidy_str, echo = TRUE,  eval = FALSE-------------------------------
#  'data.frame':	8833 obs. of  4 variables:
#   $ time     : POSIXlt, format: "2017-04-15" "2017-04-15" ...
#   $ latitude : num [1:8833(1d)] 36 36.1 36.1 36.1 36.2 ...
#   $ longitude: num [1:8833(1d)] -125 -125 -125 -125 -125 ...
#   $ chla     : num  0.194 0.193 0.186 0.18 0.163 ...
#  

## ----dateline_3D, echo = TRUE,  eval = FALSE----------------------------------
#  dataInfo <- rerddap::info('jplMURSST41mday')
#  parameter <- 'sst'
#  xcoord <- c(175, 185)
#  ycoord <- c(40, 50)
#  tcoord <- c('2019-03-16', '2019-03-16')
#  mur_dateline <- rxtracto_3D(dataInfo, parameter, xcoord = xcoord, ycoord = ycoord,
#                         tcoord = tcoord)

## ----world2hires, echo = TRUE, eval = FALSE-----------------------------------
#  xlim <- c(170, 190)
#  ylim <- c(40, 55)
#  remove <- c("UK:Great Britain", "France", "Spain", "Algeria", "Mali", "Burkina Faso", "Ghana", "Togo")
#  w <- map("world2Hires", xlim = xlim, ylim = ylim, fill = TRUE, plot = FALSE)
#  w <- map("mapdata::world2Hires", regions = w$names[!(w$names %in% remove)], plot = FALSE, fill = TRUE, ylim = ylim, xlim = xlim)
#  

## ----world2hires_map, echo = TRUE, eval = FALSE-------------------------------
#  mapFrame <- function(longitude, latitude, my_data) {
#    my_data_name <- names(my_data)
#    temp_data <- drop(my_data[[1]])
#    dims <- dim(temp_data)
#    temp_data <- array(temp_data, dims[1] * dims[2])
#    my_frame <- expand.grid(x = longitude, y = latitude)
#    my_frame[my_data_name] <- temp_data
#    return(my_frame)
#  }
#  mur_frame <- mapFrame(mur_dateline$longitude, mur_dateline$latitude, mur_dateline['sst'])
#  mycolor <- cmocean::cmocean('thermal')(256)
#    myplot <- ggplot(data = mur_frame, aes(x = x, y = y, fill = sst)) +
#    geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +     geom_tile(interpolate = FALSE) +
#      scale_fill_gradientn(colours = mycolor, na.value = NA) +
#      theme_bw() + ylab("latitude") + xlab("longitude") +
#      coord_fixed(1.3, xlim = xlim, ylim = ylim)
#  myplot

## ----mbnmsChla, echo = TRUE, eval = FALSE-------------------------------------
#  require("rerddapXtracto")
#  dataInfo <- rerddap::info('erdVH3chlamday')
#  parameter = 'chla'
#  tpos <- c("2014-09-01", "2014-10-01")
#  #tpos <-as.Date(tpos)
#  xpos <- mbnms$Longitude
#  ypos <- mbnms$Latitude
#  sanctchl <- rxtractogon(dataInfo, parameter = parameter, xcoord = xpos, ycoord = ypos,  tcoord = tpos)
#  str(sanctchl)

## ----mbnmsChla_struct, echo = TRUE, eval = FALSE------------------------------
#  List of 6
#   $ chla       : num [1:50, 1:57, 1:2] NA NA NA NA NA NA NA NA NA NA ...
#   $ datasetname: chr "erdVH3chlamday"
#   $ longitude  : num [1:50(1d)] -123 -123 -123 -123 -123 ...
#   $ latitude   : num [1:57(1d)] 35.6 35.6 35.6 35.7 35.7 ...
#   $ altitude   : logi NA
#   $ time       : POSIXlt[1:2], format: "2014-09-15" "2014-10-15"
#   - attr(*, "class")= chr [1:2] "list" "rxtracto3D"
#  

## ----mbnmsChlaPlot, echo = TRUE,  eval = FALSE--------------------------------
#  require("ggplot2")
#  require("plotdap")
#  myFunc <- function(x) log(x)
#  sanctchl1 <- sanctchl
#  sanctchl1$chla <- sanctchl1$chla[, , 2]
#  sanctchl1$time <- sanctchl1$time[2]
#  sanctchlPlot <- plotBBox(sanctchl1, plotColor = 'algae', myFunc = myFunc)
#  sanctchlPlot

## ----animate, eval = FALSE----------------------------------------------------
#  require("gganimate")
#  #> Loading required package: gganimate
#  require("ggplot2")
#  require("plotdap")
#  myFunc <- function(x) log(x)
#  sanctchlPlot <- plotBBox(sanctchl, plotColor = 'algae', myFunc = myFunc, time = identity, animate = TRUE)

## ----mbnmsBathy, echo = TRUE,  eval = FALSE-----------------------------------
#  require("rerddap")
#  dataInfo <- rerddap::info('etopo180')
#  xpos <- mbnms$Longitude
#  ypos <- mbnms$Latitude
#  bathy <- rxtractogon(dataInfo, parameter = 'altitude', xcoord = xpos, ycoord = ypos)
#  str(bathy)

## ----mbnmsBathy_struct, echo = TRUE,  eval = FALSE----------------------------
#  List of 6
#   $ depth      : num [1:123, 1:141, 1] NA NA NA NA NA NA NA NA NA NA ...
#   $ datasetname: chr "etopo180"
#   $ longitude  : num [1:123(1d)] -123 -123 -123 -123 -123 ...
#   $ latitude   : num [1:141(1d)] 35.5 35.6 35.6 35.6 35.6 ...
#   $ altitude   : logi NA
#   $ time       : logi NA
#   - attr(*, "class")= chr [1:2] "list" "rxtracto3D"
#  

## ----mbnmsBathyPlot, echo = TRUE, eval = FALSE--------------------------------
#  require("ggplot2")
#  require("mapdata")
#  myFunc = function(x) -x
#  bathyPlot <- suppressMessages((plotBBox(bathy, plotColor = 'dense', myFunc = myFunc, name = 'Depth')))
#  bathyPlot

## ----soda70, echo = TRUE,  eval = FALSE---------------------------------------
#  require("rerddap")
#  dataInfo <- rerddap::info('erdSoda331oceanmday')
#  xpos <- c(185.25, 240.25)
#  ypos <- c(20.25, 60.25)
#  zpos <- c(76.80285, 76.80285)
#  tpos <- c('2010-12-15', '2010-12-15')
#  soda70 <- rxtracto_3D(dataInfo, parameter = 'temp', xcoord = xpos, ycoord = ypos, tcoord = tpos, zcoord = zpos, zName = 'depth')

## ----soda70Plot, echo = TRUE, eval = FALSE------------------------------------
#  require("ggplot2")
#  require("plotdap")
#  sodaPlot <- plotBBox(soda70, plotColor = 'thermal', name = 'temp_at_70m', maxpixels = 30000)
#  sodaPlot
#  

## ----NAtlSSS, eval = FALSE, echo = TRUE---------------------------------------
#  require("rerddap")
#  urlBase <- "https://erddap.marine.ie/erddap/"
#  parameter <- "sea_surface_salinity"
#  sssTimes <- c("last", "last")
#  sssLats <- c(48.00625, 57.50625)
#  sssLons <- c(-17.99375, -1.00625)
#  dataInfo <- rerddap::info("IMI_NEATL", url = urlBase)
#  NAtlSSS <- rxtracto_3D(dataInfo, parameter = parameter, xcoord = sssLons, ycoord = sssLats, tcoord = sssTimes)
#  

## ----NAtlSSSplot, eval = FALSE, echo = TRUE-----------------------------------
#  require("ggplot2")
#  require("plotdap")
#  NAtlSSSPlot <- plotBBox(NAtlSSS, plotColor = 'haline', name = "salinity", maxpixels = 30000)
#  NAtlSSSPlot

## ----NAtlSSSplot1, eval = FALSE, echo = TRUE----------------------------------
#  require("ggplot2")
#  require("plotdap")
#  haline = cmocean::cmocean('haline')(256)
#  add_ggplot(NAtlSSSPlot, scale_colour_gradientn(colours = haline, na.value = NA, limits = c(32, 36)), scale_fill_gradientn(colours = haline, na.value = NA, limits = c(32, 36)))

## ----nearGrid, eval = FALSE---------------------------------------------------
#  latitude[which.min(abs(latitude - ypos[1]))]  # minimum latitude
#  latitude[which.min(abs(latitude - ypos[2]))]  # maximum latitude
#  longitude[which.min(abs(longitude- xpos[1]))] # minimum longitude
#  longitude[which.min(abs(longitude - xpos[2]))] # maximum longitude
#  isotime[which.min(abs(time - tpos[1]))] # minimum time
#  isotime[which.min(abs(time - tpos[2]))] # maximum time


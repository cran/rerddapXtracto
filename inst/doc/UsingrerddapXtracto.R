## ----initialize, echo = FALSE--------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE----------------------------------------------
### needed libraries
library(mapdata)
library(plotdap)
library(rerddap)
library(rerddapXtracto)


## ----install, eval = FALSE-----------------------------------------------
#  install.packages("ncdf4", dependencies = TRUE)
#  install.packages("parsedate", dependencies = TRUE)
#  install.packages("plotdap", dependencies = TRUE)
#  install.packages("rerddap", dependencies = TRUE)
#  install.packages("sp", dependencies = TRUE)

## ----install_package, eval = FALSE---------------------------------------
#  install.packages("rerddapXtracto", dependencies = TRUE)

## ----installGit, eval = FALSE--------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("rmendels/rerddapXtracto")

## ---- eval = FALSE-------------------------------------------------------
#  library("rerddapXtracto")

## ---- eval = FALSE-------------------------------------------------------
#  library("gganimate")
#  library("ggplot2")
#  library("plotdap")

## ----info----------------------------------------------------------------
require("rerddap")
## base URL does not need to given because it is the default one
dataInfo <- info('erdMBchla1day')
dataInfo

## ----getMarlinChl, eval = FALSE, echo = TRUE-----------------------------
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
#  swchl1 <- rxtracto(swchlInfo, parameter = 'chlorophyll', xcoord = xpos, ycoord = ypos, tcoord = tpos, zcoord = zpos, xlen = .2, ylen = .2)

## ----meantrackPlot, eval = FALSE, echo = TRUE----------------------------
#  require("ggplot2")
#  require("plotdap")
#  
#  myPlot <- plotTrack(swchl1, xpos, ypos, tpos, plotColor = 'chlorophyll')
#  myPlot

## ----animateTrack, echo = TRUE, eval = FALSE-----------------------------
#  myPlot <- plotTrack(swchl1, xpos, ypos, tpos, plotColor = 'chlorophyll',
#                      animate = TRUE, cumulative = TRUE)
#  

## ----topotag, eval = FALSE, echo = TRUE----------------------------------
#  require("ggplot2")
#  require("plotdap")
#  require("rerddap")
#  require("rerddapXtracto")
#  ylim <- c(15, 30)
#  xlim <- c(-160, -105)
#  topoInfo <- rerddap::info('etopo360')
#  topo <- rxtracto(topoInfo, parameter = 'altitude', xcoord = xpos, ycoord = ypos, xlen = .1, ylen = .1)
#  myFunc = function(x) -x
#  topoPlot <- plotTrack(topo, xpos, ypos, NA, plotColor = 'density', name = 'Depth', myFunc = myFunc)
#  topoPlot

## ----extract3D-----------------------------------------------------------
require("rerddap")
urlBase <- "https://erddap.marine.ie/erddap/"
parameter <- "Sea_water_temperature"
dataInfo <- rerddap::info("IMI_CONN_3D", url = urlBase)
#get the actual last 3 times,  and extract from data frame
dataInfo1 <- read.csv("https://erddap.marine.ie/erddap/griddap/IMI_CONN_3D.csv0?time[last-2:1:last]",stringsAsFactors = FALSE, header = FALSE, row.names = NULL)
sstTimes <- dataInfo1[[1]]
sstLats <- c(53.505758092414446, 53.509303546859805, 53.51284900130517)
sstLons <- c(-10.25975390624996, -10.247847656249961, -10.23594140624996)
sstDepths <- c(2, 6, 10)
sstTrack <- rxtracto(dataInfo, parameter = parameter, xcoord = sstLons, ycoord = sstLats, tcoord = sstTimes, zcoord = sstDepths, xlen = .05, ylen = .05, zlen = 0., zName = 'altitude')
str(sstTrack)

## ----dateline_track------------------------------------------------------
dataInfo <- rerddap::info('jplMURSST41mday')
parameter <- 'sst'
xcoord <- c(179.7, 179.8, 179.9, 180., 180.1, 180.2, 180.3, 180.4)
ycoord <- c(40, 40, 40, 40, 40, 40, 40, 40)
tcoord <- c('2018-03-16', '2018-03-16', '2018-03-16','2018-03-16','2018-03-16','2018-03-16','2018-03-16','2018-03-16')
xlen <- .05
ylen <- .05
extract <- rxtracto(dataInfo, parameter = parameter, xcoord = xcoord,
                    ycoord = ycoord, tcoord = tcoord,
                    xlen = xlen, ylen = ylen)
str(extract)

## ----VIIRSchla, warning = FALSE,  message = FALSE------------------------
require("rerddap")
require("rerddapXtracto")

xpos <- c(-125, -120) 
ypos <- c(39, 36)
tpos <- c("last", "last")
tpos <- c("2017-04-15", "2017-04-15")
VIIRSInfo <- rerddap::info('erdVH3chlamday')
VIIRS <- rxtracto_3D(VIIRSInfo, parameter = 'chla', xcoord = xpos, ycoord = ypos, tcoord = tpos)

## ----VIIRSLogPlot, fig.width = 5, fig.height = 5, fig.align = 'center', warning = FALSE----
require("ggplot2")
require("plotdap")
myFunc <- function(x) log(x)
chlalogPlot <- plotBBox(VIIRS, plotColor = 'chlorophyll', myFunc = myFunc)
chlalogPlot

## ----dateline_3D, echo = TRUE,  eval = FALSE-----------------------------
#  dataInfo <- rerddap::info('jplMURSST41mday')
#  parameter <- 'sst'
#  xcoord <- c(175, 185)
#  ycoord <- c(40, 50)
#  tcoord <- c('2019-03-16', '2019-03-16')
#  mur_dateline <- rxtracto_3D(dataInfo, parameter, xcoord = xcoord, ycoord = ycoord,
#                         tcoord = tcoord)

## ----world2hires, echo = TRUE, eval = FALSE------------------------------
#  xlim <- c(170, 190)
#  ylim <- c(40, 55)
#  remove <- c("UK:Great Britain", "France", "Spain", "Algeria", "Mali", "Burkina Faso", "Ghana", "Togo")
#  w <- map("world2Hires", xlim = xlim, ylim = ylim, fill = TRUE, plot = FALSE)
#  w <- map("mapdata::world2Hires", regions = w$names[!(w$names %in% remove)], plot = FALSE, fill = TRUE, ylim = ylim, xlim = xlim)
#  

## ----world2hires_map, echo = TRUE, eval = FALSE--------------------------
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
#  mycolor <- cmocean$thermal
#    myplot <- ggplot(data = mur_frame, aes(x = x, y = y, fill = sst)) +
#    geom_polygon(data = w, aes(x = long, y = lat, group = group), fill = "grey80") +     geom_raster(interpolate = FALSE) +
#      scale_fill_gradientn(colours = mycolor, na.value = NA) +
#      theme_bw() + ylab("latitude") + xlab("longitude") +
#      coord_fixed(1.3, xlim = xlim, ylim = ylim)
#  myplot

## ----mbnmsChla-----------------------------------------------------------
require("rerddapXtracto")
dataInfo <- rerddap::info('erdVH3chlamday')
parameter = 'chla'
tpos <- c("2014-09-01", "2014-10-01")
#tpos <-as.Date(tpos)
xpos <- mbnms$Longitude
ypos <- mbnms$Latitude
sanctchl <- rxtractogon(dataInfo, parameter = parameter, xcoord = xpos, ycoord = ypos,  tcoord = tpos)
str(sanctchl)

## ----mbnmsChlaPlot, fig.width = 6, fig.height = 3, fig.align = 'center', warning = FALSE----
require("ggplot2")
require("plotdap")
myFunc <- function(x) log(x)
sanctchl1 <- sanctchl
sanctchl1$chla <- sanctchl1$chla[, , 2]
sanctchl1$time <- sanctchl1$time[2]
sanctchlPlot <- plotBBox(sanctchl1, plotColor = 'chlorophyll', myFunc = myFunc)
sanctchlPlot

## ----animate, eval = FALSE-----------------------------------------------
#  require("gganimate")
#  #> Loading required package: gganimate
#  require("ggplot2")
#  require("plotdap")
#  myFunc <- function(x) log(x)
#  sanctchlPlot <- plotBBox(sanctchl, plotColor = 'chlorophyll', myFunc = myFunc, time = identity, animate = TRUE)

## ----mbnmsBathy, warning = FALSE-----------------------------------------
require("rerddap")
dataInfo <- rerddap::info('etopo180')
xpos <- mbnms$Longitude
ypos <- mbnms$Latitude
bathy <- rxtractogon(dataInfo, parameter = 'altitude', xcoord = xpos, ycoord = ypos)
str(bathy)

## ----mbnmsBathyPlot, fig.width = 5, fig.height = 5, fig.align = 'center', warning = FALSE, message = FALSE----
require("ggplot2")
require("mapdata")
myFunc = function(x) -x
bathyPlot <- plotBBox(bathy, plotColor = 'density', myFunc = myFunc, name = 'Depth')
bathyPlot

## ----soda70--------------------------------------------------------------
require("rerddap")
dataInfo <- rerddap::info('erdSoda331oceanmday')
xpos <- c(185.25, 240.25)
ypos <- c(20.25, 60.25)
zpos <- c(76.80285, 76.80285)
tpos <- c('2010-12-15', '2010-12-15')
soda70 <- rxtracto_3D(dataInfo, parameter = 'temp', xcoord = xpos, ycoord = ypos, tcoord = tpos, zcoord = zpos, zName = 'depth')
str(soda70)

## ----soda70Plot, fig.width = 6, fig.height = 3, fig.align = 'center', warning = FALSE----
require("ggplot2")
require("plotdap")
sodaPlot <- plotBBox(soda70, plotColor = 'temperature', name = 'temp_at_70m', maxpixels = 30000)
sodaPlot


## ----NAtlSSS, eval = FALSE, echo = TRUE----------------------------------
#  require("rerddap")
#  urlBase <- "https://erddap.marine.ie/erddap/"
#  parameter <- "sea_surface_salinity"
#  sssTimes <- c("last", "last")
#  sssLats <- c(48.00625, 57.50625)
#  sssLons <- c(-17.99375, -1.00625)
#  dataInfo <- rerddap::info("IMI_NEATL", url = urlBase)
#  NAtlSSS <- rxtracto_3D(dataInfo, parameter = parameter, xcoord = sssLons, ycoord = sssLats, tcoord = sssTimes)
#  

## ----NAtlSSSplot, eval = FALSE, echo = TRUE------------------------------
#  require("ggplot2")
#  require("plotdap")
#  NAtlSSSPlot <- plotBBox(NAtlSSS, plotColor = 'salinity', name = "salinity", maxpixels = 30000)
#  NAtlSSSPlot

## ----NAtlSSSplot1, eval = FALSE, echo = TRUE-----------------------------
#  require("ggplot2")
#  require("plotdap")
#  add_ggplot(NAtlSSSPlot, scale_colour_gradientn(colours = cmocean$haline, na.value = NA, limits = c(32, 36)), scale_fill_gradientn(colours = cmocean$haline, na.value = NA, limits = c(32, 36)))

## ----IFREMER-------------------------------------------------------------
require("rerddap")
urlBase <- "https://www.ifremer.fr/erddap/"
parameter <- "PSAL"
ifrTimes <- c("2013-09-15", "2013-09-15")
ifrLats <- c(30., 50.)
ifrLons <- c(-140., -110.)
ifrDepth <- c(75., 75.)
dataInfo <- rerddap::info("ifremer_tds0_6080_109e_ed80", url = urlBase)
ifrPSAL <- rxtracto_3D(dataInfo, parameter = parameter, xcoord = ifrLons, ycoord = ifrLats, tcoord = ifrTimes, zcoord = ifrDepth, zName = 'depth')
str(ifrPSAL)

## ----ifrPSALplot, fig.width = 6, fig.height = 3, fig.align='center', warning = FALSE----
require("ggplot2")
require("plotdap")
ifrPSALPlot <- plotBBox(ifrPSAL, plotColor = 'salinity', name = "salinity", maxpixels = 30000)
ifrPSALPlot

## ----nearGrid, eval = FALSE----------------------------------------------
#  latitude[which.min(abs(latitude - ypos[1]))]  # minimum latitude
#  latitude[which.min(abs(latitude - ypos[2]))]  # maximum latitude
#  longitude[which.min(abs(longitude- xpos[1]))] # minimum longitude
#  longitude[which.min(abs(longitude - xpos[2]))] # maximum longitude
#  isotime[which.min(abs(time - tpos[1]))] # minimum time
#  isotime[which.min(abs(time - tpos[2]))] # maximum time


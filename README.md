# rerddapXtracto
R package for accessing environmental data using 'rerddap' and 'ERDDAP' servers

`rerddapXtracto` is an 'R' developed to subset and extract satellite and other oceanographic related data from a remote 'ERDDAP' server. The program can extract data for a moving point in time along a user-supplied set of longitude, latitude and time points; in a 3D bounding box; or within a polygon (through time).  

These functions differ from those in the 'xtractomatic' package in that they use the 'rerddap' package to access gridded data on any 'ERDDAP' server, but they require the user to provide initial information about the data to be extracted.


There are three main data extraction functions in the package: 

- rxtracto <- function(dataInfo, parameter = NULL, 
                  xcoord=NULL, ycoord=NULL, zcoord = NULL, tcoord=NULL,
                  xlen=0., ylen=0., 
                  xName='longitude', yName='latitude', zName='altitude', tName='time',
                  verbose = FALSE)

- rxtracto_3D <- function(dataInfo, parameter = NULL, 
                     xcoord=NULL, ycoord=NULL, zcoord = NULL, tcoord=NULL,
                     xName='longitude', yName='latitude', zName='altitude', tName='time',
                     verbose = FALSE)

- rxtractogon <- function(dataInfo, parameter, 
                     xcoord=NULL, ycoord=NULL, zcoord = NULL, tcoord=NULL, 
                     xName='longitude', yName='latitude', zName='altitude', tName='time', 
                     verbose = FALSE)


and two functions for producing maps using the 'plotdap' package:

- plotTrack <- function(resp, xcoord, ycoord, plotColor = 'viridis', 
                   animate = FALSE, cumulative = FALSE, 
                   name = NA, myFunc = NA, shape = 20, size = .5)

- plotBBox <- function(resp, plotColor = 'viridis', time = NA,
                  animate = FALSE, cumulative = FALSE,
                  name = NA, myFunc = NA, maxpixels = 10000)


`rerddapXtracto` uses the `ncdf4`, `plotdap`, `rerddap`,  and `sp` packages , and these packages (and the packages imported by these packages) must be installed first or `rerddapXtracto` will fail to install.   

```r
install.packages("ncdf4") 
install.packages("plotdap") 
install.packages("rerddap", dependencies = TRUE)
install.packages("sp")
```

The `rerddapXtracto` package can be installed from CRAN:

```r 
install.packages("rerddapXtracto")
```

or the development version `rerddapXtracto` package at the moment can be installed from Github using the devtools package:

```r
install.packages("devtools")
devtools::install_github("rmendels/rerddapXtracto")
```



# Required legalese

“The United States Department of Commerce (DOC) code is provided
on an "as is" basis and the user assumes responsibility for its use.
DOC has relinquished control of the information and no longer has responsibility
to protect the integrity, confidentiality, or availability of the information.
Any claims against the Department of Commerce stemming from the use of its
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their endorsement,
recommendation or favoring by the Department of Commerce. The Department of
Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used
in any manner to imply endorsement of any commercial product or activity by DOC
or the United States Government.”



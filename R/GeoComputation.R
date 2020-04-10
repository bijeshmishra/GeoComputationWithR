# Set working directory:
setwd("~/Documents/GitHub/GeoComputationWithR")
# ##### Contribute  (Not Necesssary for Class) #####
# devtools::install_github("geocompr/geocompkg") #Reporduce code in the book.
# library(bookdown)
# bookdown::render_book("index.Rmd") # to build the book locally.
# browseURL("_book/index.html") # to view it.


# ##### Chapter 1 #####
# Introduction
# 1.1 What is geocomputation?: Closely related to Geographic Information science, geomatics, geoinformatics, spatial information science, geoinfomatiion engineering and geographical data science (GDS).

# 1.2 Why use R for geocomputation?
# install.packages("leaflet")
library(leaflet)
popup = c("Bijesh", "Rajesh", "Shyam")
leaflet() %>%
  addProviderTiles("NASAGIBS.ViirsEarthAtNight2012") %>%
  addMarkers(lng = c(-3, 23, 11),
             lat = c(52, 53, 49),
             popup = popup)

# 1.3 Software for geocomputation: QGIS, GRASS and SAGA built on C++, C++ is fasst but hard to learn compared to python. Rcpp made C++ more accessible. Java is another important program.Python is important program for geocomputing. Packages: osgeo, Shapely, NumPy, PyGeoProcessing.

# 1.4 R’s spatial ecosystem
# 1.5 The history of R-spatial

# 1.6 Exercises
# Q1. Think about the terms "GIS" and "geocomputation" described above. Which (if any) best describe the work you would like to do  using geo* methods and software and why?
# Q2. Provide three reasons for using a scriptable language such as R for geocomputation instead of using an established GIS program such as QGIS?
# Q3. Name two advantages and two disadvantages of using mature vs recent packages for geographic data analysis (for example sf vs sp).
# Note: Python modules providing access to geoalgorithms include grass.script for GRASS, saga-python for SAGA-GIS, processing for QGIS and arcpy for ArcGIS.

# An overview of R’s spatial ecosystem can be found in the CRAN Task View on the Analysis of Spatial Data (see https://cran.r-project.org/web/views/Spatial.html).↩

# ##### Chapter 2: Geographic Data in R #####
# # Prerequisites:
# # Install Packages:
# # this step is to install library from github:
# install.packages("devtools") #install devtools package. Do not run of devtools is alread present in your machine.
# library(devtools) #load devtools before downloading packages from github.
# install_github("r-spatial/sf") #isntall sf package from github
# devtools::install_github("r-spatial/sf", type = "binary", force = TRUE) #This works without attaching devtools library.

# install.packages("raster")
# install.packages("spData")
# install.packages("devtools") #Install this package if you are installing packages from GitHub.
# devtools::install_github("Nowosad/spDataLarge", force = TRUE) #Need to add "Force = TRUE" to install this package in my MAC.
# install.packages("sf", type = "binary") #OSU IT People response: There is no such non-binary version. Homebrew and gdal packages are only necessary if you want to modify package. For, user, binary version is enough.

# Load Libraries.
library(raster) #class ad function for raster data.
library(spData) #load geographic data.
library(devtools) #used to install packages from GitHub.
library(spDataLarge) #Load larger geograpic Data.
library(sf) #class and function for vector data

# # Load data from SpDataLarge if you cannot install the package. Code from Github.
# if(!require(spDataLarge)) {
#   download.file("https://github.com/Nowosad/spDataLarge/archive/master.zip", "spDataLarge.zip")
#   unzip("spDataLarge.zip")
#   files_rda = list.files("spDataLarge-master/data/", full.names = TRUE)
#   sapply(files_rda, load, envir = .GlobalEnv)
# }

# # 2.1 Introduction
# # 2.2 Vector data

# 2.2.1 An introduction to simple features
vignette(package = "sf") # see which vignettes/descriptions are available
vignette("sf1")          # an introduction to the package
names(world)
plot(world)
summary(world["lifeExp"])
world_mini = world[1:2, 1:3]
world_mini

# 2.2.2 Why simple features?
library(sp)
world_sp = as(world, Class = "Spatial")
world_sf = st_as_sf(world_sp)

# 2.2.3 Basic map making
plot(world[3:6])
plot(world["pop"])
world_asia = world[world$continent == "Asia", ]
asia = st_union(world_asia)
plot(world["pop"], reset = FALSE)
plot(asia, add = TRUE, col = "red")

# 2.2.4 Base plot arguments
plot(world["continent"], reset = FALSE)
cex = sqrt(world$pop) / 10000

world_cents = st_centroid(world, of_largest = TRUE)
# # I got this warning message after running above line of code. This might the resaon why map of India shrunk down and moved to mid-bottom of plot.
# Warning messages:
#   1: In st_centroid.sf(world, of_largest = TRUE) :
#   st_centroid assumes attributes are constant over geometries of x
# 2: In st_centroid.sfc(st_geometry(x), of_largest_polygon = of_largest_polygon) :
#   st_centroid does not give correct centroids for longitude/latitude data

plot(st_geometry(world_cents), add = TRUE, cex = cex)
india = world[world$name_long == "India", ]
plot(st_geometry(india), expandBB = c(0, 0.2, 0.1, 1), col = "gray", lwd = 3) #When I run this code, the map shrink down to really tiny size and moves down to mid-bottom of plotting region.
plot(world_asia[0], add = TRUE)

# # 2.2.5 Geometry types
# # seven most commonly used types: POINT, LINESTRING, POLYGON, MULTIPOINT, MULTILINESTRING, MULTIPOLYGON and GEOMETRYCOLLECTION.

# # 2.2.6 Simple feature geometries (sfg)
# # A point: st_point()
# # A linestring: st_linestring()
# # A polygon: st_polygon()
# # A multipoint: st_multipoint()
# # A multilinestring: st_multilinestring()
# # A multipolygon: st_multipolygon()
# # A geometry collection: st_geometrycollection()

st_point(c(5, 2))                 # XY point
st_point(c(5, 2, 3))              # XYZ point
st_point(c(5, 2, 1), dim = "XYM") # XYM point
st_point(c(5, 2, 3, 1))           # XYZM point

# the rbind function simplifies the creation of matrices
# MULTIPOINT
multipoint_matrix = rbind(c(5, 2), c(1, 3), c(3, 4), c(3, 2))
st_multipoint(multipoint_matrix)
# LINESTRING
linestring_matrix = rbind(c(1, 5), c(4, 4), c(4, 1), c(2, 2), c(3, 2))
st_linestring(linestring_matrix)
# POLYGON
polygon_list = list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5)))
st_polygon(polygon_list)
# POLYGON with a hole
polygon_border = rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5))
polygon_hole = rbind(c(2, 4), c(3, 4), c(3, 3), c(2, 3), c(2, 4))
polygon_with_hole_list = list(polygon_border, polygon_hole)
st_polygon(polygon_with_hole_list)
# MULTILINESTRING
multilinestring_list = list(rbind(c(1, 5), c(4, 4), c(4, 1), c(2, 2), c(3, 2)),
                            rbind(c(1, 2), c(2, 4)))
st_multilinestring((multilinestring_list))
# MULTIPOLYGON
multipolygon_list = list(list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5))),
                         list(rbind(c(0, 2), c(1, 2), c(1, 3), c(0, 3), c(0, 2))))
st_multipolygon(multipolygon_list)
# GEOMETRYCOLLECTION
gemetrycollection_list = list(st_multipoint(multipoint_matrix),
                              st_linestring(linestring_matrix))
st_geometrycollection(gemetrycollection_list)

# 2.2.7 Simple feature columns (sfc)
# sfc POINT
point1 = st_point(c(5, 2))
point2 = st_point(c(1, 3))
points_sfc = st_sfc(point1, point2)
points_sfc

# sfc POLYGON
polygon_list1 = list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5)))
polygon1 = st_polygon(polygon_list1)
polygon_list2 = list(rbind(c(0, 2), c(1, 2), c(1, 3), c(0, 3), c(0, 2)))
polygon2 = st_polygon(polygon_list2)
polygon_sfc = st_sfc(polygon1, polygon2)
st_geometry_type(polygon_sfc)

# sfc MULTILINESTRING
multilinestring_list1 = list(rbind(c(1, 5), c(4, 4), c(4, 1), c(2, 2), c(3, 2)),
                             rbind(c(1, 2), c(2, 4)))
multilinestring1 = st_multilinestring((multilinestring_list1))
multilinestring_list2 = list(rbind(c(2, 9), c(7, 9), c(5, 6), c(4, 7), c(2, 7)),
                             rbind(c(1, 7), c(3, 8)))
multilinestring2 = st_multilinestring((multilinestring_list2))
multilinestring_sfc = st_sfc(multilinestring1, multilinestring2)
st_geometry_type(multilinestring_sfc)

# sfc GEOMETRY
point_multilinestring_sfc = st_sfc(point1, multilinestring1)
st_geometry_type(point_multilinestring_sfc)

st_crs(points_sfc) #Coordinate Reference System: NA.

# EPSG definition
points_sfc_wgs = st_sfc(point1, point2, crs = 4326)
st_crs(points_sfc_wgs)

# PROJ4STRING definition
st_sfc(point1, point2, crs = "+proj=longlat +datum=WGS84 +no_defs") # Note: Sometimes st_crs() will return a proj4string but not an epsg code.

# 2.2.8 The sf class
lnd_point = st_point(c(0.1, 51.5))                 # sfg object
lnd_geom = st_sfc(lnd_point, crs = 4326)           # sfc object
lnd_attrib = data.frame(                           # data.frame object
  name = "London",
  temperature = 25,
  date = as.Date("2017-06-21")
)
lnd_sf = st_sf(lnd_attrib, geometry = lnd_geom)    # sf object
lnd_sf
class(lnd_sf)

# 2.3 Raster data
# 2.3.1 An introduction to raster
raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
new_raster = raster(raster_filepath)
new_raster
dim(new_raster) #Number of rw=ows, columns and layers.
ncell(new_raster) #number of cells (pixels)
extent(new_raster) #raster's spatial resolution
crs(new_raster) #coordinate reference system of raster
help("raster-package") # returns a full list of all available raster functions.

# 2.3.2 Basic map making
plot(new_raster)
# 2.3.3 Raster classes
# # Raster package supports numerous drivers with the help of rdgal. To check which drivers are availale, run two commands below:
# raster::writeFormats()
# rgdal::gdalDrivers()
raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
new_raster = raster(raster_filepath)
new_raster2 = raster(nrows = 6, ncols = 6, res = 0.5,
                     xmn = -1.5, xmx = 1.5, ymn = -1.5, ymx = 1.5,
                     vals = 1:36)
multi_raster_file = system.file("raster/landsat.tif", package = "spDataLarge")
r_brick = brick(multi_raster_file)
r_brick
nlayers(r_brick)
raster_on_disk = raster(r_brick, layer = 1)
raster_in_memory = raster(xmn = 301905, xmx = 335745,
                          ymn = 4111245, ymx = 4154085,
                          res = 30)
values(raster_in_memory) = sample(seq_len(ncell(raster_in_memory)))
crs(raster_in_memory) = crs(raster_on_disk)
r_stack = stack(raster_in_memory, raster_on_disk)
r_stack

# 2.4 Coordinate Reference Systems
# 2.4.1 Geographic coordinate systems
# 2.4.2 Projected coordinate reference systems
# 2.4.3 CRSs in R
crs_data = rgdal::make_EPSG()
View(crs_data)
vector_filepath = system.file("vector/zion.gpkg", package = "spDataLarge")
new_vector = st_read(vector_filepath)
st_crs(new_vector) # get CRS
new_vector = st_set_crs(new_vector, 4326) # set CRS
projection(new_raster) # get CRS
# 2.4
# 2.5
# Excercise 4 and 5
#Answer question in R or in word file.

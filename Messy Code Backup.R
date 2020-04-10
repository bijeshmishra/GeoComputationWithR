# Set working directory:
setwd("~/Documents/GitHub/GeoComputationWithR")
# ###### Contribute  (Not Necesssary for Class) #####
# devtools::install_github("geocompr/geocompkg") #Reporduce code in the book.
# library(bookdown)
# bookdown::render_book("index.Rmd") # to build the book locally.
# browseURL("_book/index.html") # to view it.


# ##### Install and Load Packages (Part of Ex.12; needed later as well) #####
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
library(sp) #spatial feature.

# # Load data from SpDataLarge if you cannot install the package. Code from Github.
# if(!require(spDataLarge)) {
#   download.file("https://github.com/Nowosad/spDataLarge/archive/master.zip", "spDataLarge.zip")
#   unzip("spDataLarge.zip")
#   files_rda = list.files("spDataLarge-master/data/", full.names = TRUE)
#   sapply(files_rda, load, envir = .GlobalEnv)
# }


# ##### Excercise 12: Spatial Data in R: Chapter 1 #####
# Introduction
# 1.1 What is geocomputation?: Closely related to Geographic Information science, geomatics, geoinformatics, spatial information science, geoinfomatiion engineering and geographical data science (GDS).

# 1.2 Why use R for geocomputation?
# install.packages("leaflet")
library(leaflet)
popup = c("Robin", "Jakub", "Jannes")
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


# ##### Excercise 12: Chapter 2: Geographic Data in R #####
# # 2.1 Introduction
# # 2.2 Vector data
# London coordinates: c(-0.1, 51.5) is (longitude/Prime Meridiaan = degree east, Latitude/Equator = degree north) of the origin.
# In Projected CRS, above point = c(530000, 180000) is (Easting/Northing) means London is 530 km East and 180 km North of the norigin of CRS. This is British National Grid System.

# 2.2.1 An introduction to simple features:
vignette(package = "sf") # see which vignettes/descriptions are available
vignette("sf1")          # an introduction to the package
names(world) #world is spatial object containing spatial and attribute columns returned by this function.
plot(world) #plot world.
world$geom # list of columns that contains all coordinates of the country polygons.
summary(world["lifeExp"]) # geom is stick with the desired variable unless it is removed. MULTIPOLYGON is seen for countries with Island.
world_mini = world[1:2, 1:3] #A subset of first two row (from row 1 to rown 2) and three columns (from column 1 to column 3). This adds additional geographical data (geometry type, dimension, bbox) and CRS information ( epsg(SRID), proj4string ) and the presence of a geometry column (geom).
world_mini
View(world)

# 2.2.2 Why simple features (sf)? #Advantage: Cross transferrable.Also supports %>%, tmap, mapview, tidycensus.
library(sp)
world_sp = as(world, Class = "Spatial") # convert to spatial class.
world_sf = st_as_sf(world_sp, "sf") # spatial object convert back to sf.
world_sf

# 2.2.3 Basic map making
plot(world[3:6]) # similar to spplot() in sp package. Set color with "col = <color name>"
plot(world["pop"])
world_asia = world[world$continent == "Asia", ]
asia = st_union(world_asia)
plot(world["pop"], reset = FALSE)
plot(asia, add = TRUE, col = "red")

# 2.2.4 Base plot arguments
?graphics ::plot #help
?par #help

plot(world["continent"], reset = FALSE) #main = "abc" gives as title name as abc.
cex = sqrt(world$pop) / 10000 #set diameter.
world_cents = st_centroid(world, of_largest = TRUE)
plot(st_geometry(world_cents), add = TRUE, cex = cex)

india = world[world$name_long == "India", ]
plot(st_geometry(india), expandBB = c(0, 0.2, 0.1, 1), col = "gray", lwd = 3) # When I run this code, the map shrink down to really tiny size and moves down to mid-bottom of plotting region.
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
points_sfc = st_sfc(point1, point2) #combine two simple features into one object with two features.
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

#Create an sfc object from sfg object with differnet geometry types:
# sfc GEOMETRY
point_multilinestring_sfc = st_sfc(point1, multilinestring1)
st_geometry_type(point_multilinestring_sfc)
# epsg(SRID)
# proj4string()
st_crs(points_sfc) #Coordinate Reference System: NA.

# EPSG definition: add coordinate reference system.
points_sfc_wgs = st_sfc(point1, point2, crs = 4326)
st_crs(points_sfc_wgs)

# PROJ4STRING definition
st_sfc(point1, point2, crs = "+proj=longlat +datum=WGS84 + no_defs") # Note: Sometimes st_crs() will return a proj4string but not an epsg code.

# 2.2.8 The sf class: the coordinate were used to create simple feature geometry, geometry was converted into a simple feature geometry column with a CRS and attributes were stored in a data.frame, which was combined with sfc object with st_fc().
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
vignette(package = "raster") # ignette("functions", package = "raster") did not work.
raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
new_raster = raster(raster_filepath)
new_raster #srtm.tif is digital elevation model of this area.
dim(new_raster) #Number of rows, columns and layers.
ncell(new_raster) #number of cells (pixels)
extent(new_raster) #raster's spatial resolution
crs(new_raster) #coordinate reference system of raster
inMemory(new_raster) #gives whether data is stored in memory (default) or on disk.
help("raster-package") # returns a full list of all available raster functions.

# 2.3.2 Basic map making
plot(new_raster) #other ways to plot raster data are ssplot(), levelplot().

# 2.3.3 Raster classes
# # Raster package supports numerous drivers with the help of rdgal. To check which drivers are availale, run two commands below:
# raster::writeFormats()
# rgdal::gdalDrivers()
raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
new_raster = raster(raster_filepath)
# this creates new raster file with 36 cells, centered around the prime neridian and equator (xmn, xmx, ymn, ymx).
new_raster2 = raster(nrows = 6, ncols = 6, res = 0.5,
                     xmn = -1.5, xmx = 1.5, ymn = -1.5, ymx = 1.5,
                     vals = 1:36)
?raster
# RasterBrick, given by brick() is a multiple layers which is similar to single multispectral stellite file.
multi_raster_file = system.file("raster/landsat.tif", package = "spDataLarge")
r_brick = brick(multi_raster_file)
r_brick
nlayers(r_brick) #numbers of layers tored in a Raster* object.

#RasterStack: similar to RasterLayer; stack of RasterLayer objects with same extent and resoultion. It allows to connect several raster objects stored in different files or multiply objects in memory.
raster_on_disk = raster(r_brick, layer = 1)
raster_in_memory = raster(xmn = 301905, xmx = 335745,
                          ymn = 4111245, ymx = 4154085,
                          res = 30)
values(raster_in_memory) = sample(seq_len(ncell(raster_in_memory)))
crs(raster_in_memory) = crs(raster_on_disk)
r_stack = stack(raster_in_memory, raster_on_disk)
r_stack


# ##### Excercise 13 (Coordinate Reference System and Onward) #####
# # 2.4 Coordinate Reference Systems (CRS): How spatial elements of data relate to the surface of earth. Either Gegographic or Projected.

# # 2.4.1 Geographic coordinate systems: Distance measured in degree of angle from the Prime Meridian Plane. Idenfity location using Longitude (East-West Direction or Prime Meridian Plane) and Latitude (North South Direction or Equatorial Plane).
# # Represented by spherical model (a perfect sphere of given radius; simple but rarely used) or ellipsodial model (defined by two parameters: equatorial radius (~ 11.5 km longer than polar radius) & polar radius).
# # Datum: Two types 1) Local (Eg. NAD83) and 2) Geocentric (Eg. WSG84).
st_proj_info (type = "datum") #available datum definitions.

# # 2.4.2 Projected coordinate reference systems
# # Based on caretesian coordinates on flat surface. Have origin x and y axis and limear unit of measurement (meter). Based on Geographical coordinate system.
# # Projected coordinate system are named based on what they preserve: equal-area preserves area, azimuthal preserves direction, equidistnce preserves distance and conformal preserves local shape.
# # Three main group of Projected coordinate system: conical, cylindrical and planar. Conical: earth surface is projected onto a cone (best for mid latitude), cylindrical projection maps surface into cylinder (good for mapping entire earth), planar projection projects surface onto flat surface (used for mapping polar region.)
st_proj_info(type = "proj") #list of available projections supported by PROJ library.

# # 2.4.3: CRS in R:
# # two ways: epsg (shorter, easier, one CRS) or proj4string (more flexible to speficy projection type, the datum and the ellipsoid, more complicated).
crs_data = rgdal::make_EPSG()  # Shows available CRSs.
View(crs_data)
vector_filepath = system.file("vector/zion.gpkg", package = "spDataLarge")
new_vector = st_read(vector_filepath) #polygon representing borers of Zion National Park (?zion)
st_crs(new_vector) # get CRS of an object.
new_vector = st_set_crs(new_vector, 4326) # set CRS if CRS is wrong or missing). This does not reproject data or transoform data into another projection. (Only fly-on-projection).
projection(new_raster) # get CRS from  Raster* object.Raster object only accept proj4 definitions.
projection(new_raster) = "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0
                            +units=m +no_defs" # set CRS.

# # 2.5: Units:
# The distance, area and otehr geometric calculations in sf return values that comes with units defined by the units package.
# install.packages("lwgeom") # Note: I had to install new package "lwgeom".
library(lwgeom)
luxembourg = world[world$name_long == "Luxembourg",]
st_area(luxembourg)
attributes(st_area(luxembourg)) #
st_area(luxembourg)/1000000 # change square meters into square kilometer.
units :: set_units(st_area(luxembourg), km^2) #automatically set unit to sq. kilometer.
res(new_raster)
repr = projectRaster(new_raster, crs = "+init=epsg:26912")
res(repr) #gives vector without unit. Need to know that the unit of the UTM projection is Meters.

# # 2.6 Excercises:
# ##### Excluded in Homework #####
# # 1. Use summary() on the geometry column of the world data object. What does the output tell us about: Its geometry type? The number of countries? Its coordinate reference system (CRS)?
# # 2. Run the code that ‘generated’ the map of the world in Figure 2.5 at the end of Section 2.2.4. Find two similarities and two differences between the image on your computer and that in the book. What does the cex argument do (see ?plot)? Why was cex set to the sqrt(world$pop) / 10000? Bonus: experiment with different ways to visualize the global population.
# # Use plot() to create maps of Nigeria in context (see Section 2.2.4). Adjust the lwd, col and expandBB arguments of plot(). Challenge: read the documentation of text() and annotate the map.
# ###### Excercise 4 and 5 #####
# # 4. Create an empty RasterLayer object called my_raster with 10 columns and 10 rows. Assign random values between 0 and 10 to the new raster and plot it.
my_raster = raster(nrows = 10, ncols = 10, res = 0.5,
                   xmn = -1.5, xmx = 1.5, ymn = -1.5, ymx = 1.5,
                   vals = 1:36)
plot(my_raster)

# # 5. Read-in the raster/nlcd2011.tif file from the spDataLarge package. What kind of information can you get about the properties of this file?

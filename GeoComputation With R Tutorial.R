# # Lab 12
# ##### Chapter 1: Introduction #####
# 1.2 Why use R for geocomputation?
# ##### Chapter 2: Geographic Data in R #####
# Set working directory:

# # ##### Install Packages #####
# install.packages("devtools", dependencies = TRUE)
# install.packages("tidyverse", dependencies = TRUE)
# install.packages("sf", dependencies = TRUE)
# install.packages("raster", dependencies = TRUE)
# install.packages("spData", dependencies = TRUE)
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("Nowosad/spDataLarge")
# install.packages("sp", dependencies = TRUE)
# devtools::install_github("r-spatial/sf", type = "binary", force = TRUE)
# devtools::install_github("Nowosad/spDataLarge", force = TRUE)
# install.packages("sf", type = "binary")
# install.packages("leaflet")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages(tidyr)
# install.packages("RQGIS", dependencies = TRUE, force = TRUE, type = "binary")

# Load Libraries.
library(leaflet)
library(raster) #class ad function for raster data.
library(spData) #load geographic data.
library(devtools) #used to install packages from GitHub.
library(spDataLarge) #Load larger geograpic Data.
library(sf) #class and function for vector data
library(sp) #spatial feature
library(dplyr)
library(stringr) # for working with strings (pattern matching)
library(tidyr) # for unite() and separate()
library(RQGIS)

# Chapter 1: Introduction.
popup = c("Robin", "Jakub", "Jannes")
leaflet() %>%
  addProviderTiles("NASAGIBS.ViirsEarthAtNight2012") %>%
  addMarkers(lng = c(-3, 23, 11),
             lat = c(52, 53, 49),
             popup = popup)

# 2.1 Introduction
# 2.2 Vector data
# 2.2.1 An introduction to simple features:
# vignette(package = "sf") # see which vignettes/descriptions are available
# vignette("sf1")          # an introduction to the package
names(world)
plot(world)
summary(world["lifeExp"])
world_mini = world[1:2, 1:3]
world_mini

# 2.2.2 Why simple features?
world_sp = as(world, Class = "Spatial")
world_sf = st_as_sf(world_sp)

# 2.2.3 Basic Map Making:
plot(world[3:6])
plot(world["pop"])
world_asia = world[world$continent == "Asia",]
asia = st_union(world_asia)
plot(world["pop"], reset = FALSE)
plot(asia, add = TRUE, col = "red")

# 2.2.4: Base Plot Arguments:
plot(world["continent"], reset = FALSE)
cex = sqrt(world$pop) / 10000
world_cents = st_centroid(world, of_largest = TRUE)
plot(st_geometry(world_cents), add = TRUE, cex = cex)
india = world[world$name_long == "India",]
plot(
  st_geometry(india),
  expandBB = c(0, 0.2, 0.1, 1),
  col = "gray",
  lwd = 3
)
plot(world_asia[0], add = TRUE)
# # 2.2.5: Geometry Types: Total 17 types. Common are POINT, LINESTRING, POLYGON, MULTIPOINT, MULTILINESTRING, MULTIPOLYGON and GEOMETRYCOLLECTION.

# # 2.2.6: Simple feature geometrics (sfg):
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

multipoint_matrix = rbind(c(5, 2), c(1, 3), c(3, 4), c(3, 2)) ## MULTIPOINT
st_multipoint(multipoint_matrix)

linestring_matrix = rbind(c(1, 5), c(4, 4), c(4, 1), c(2, 2), c(3, 2)) ## LINESTRING
st_linestring(linestring_matrix)

polygon_list = list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5))) ## POLYGON
st_polygon(polygon_list)

polygon_border = rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5)) ## POLYGON with a hole
polygon_hole = rbind(c(2, 4), c(3, 4), c(3, 3), c(2, 3), c(2, 4))
polygon_with_hole_list = list(polygon_border, polygon_hole)
st_polygon(polygon_with_hole_list)

multilinestring_list = list(rbind(c(1, 5), c(4, 4), c(4, 1), c(2, 2), c(3, 2)),
                            rbind(c(1, 2), c(2, 4))) ## MULTILINESTRING
st_multilinestring((multilinestring_list))

multipolygon_list = list(list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5))),
                         list(rbind(c(0, 2), c(1, 2), c(1, 3), c(0, 3), c(0, 2)))) ## MULTIPOLYGON
st_multipolygon(multipolygon_list)

gemetrycollection_list = list(st_multipoint(multipoint_matrix),
                              st_linestring(linestring_matrix)) ## GEOMETRYCOLLECTION
st_geometrycollection(gemetrycollection_list)

# 2.2.7: Simple Feature columns (sfc):
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

st_crs(points_sfc) #coordinate reference system.

# EPSG definition
points_sfc_wgs = st_sfc(point1, point2, crs = 4326)
st_crs(points_sfc_wgs)

# PROJ4STRING definition
st_sfc(point1, point2, crs = "+proj=longlat +datum=WGS84 +no_defs")

# 2.2.8: The sf Class:
lnd_point = st_point(c(0.1, 51.5))                 # sfg object
lnd_geom = st_sfc(lnd_point, crs = 4326)           # sfc object
lnd_attrib = data.frame(# data.frame object
  name = "London",
  temperature = 25,
  date = as.Date("2017-06-21"))
lnd_sf = st_sf(lnd_attrib, geometry = lnd_geom)    # sf object
lnd_sf
class(lnd_sf)

# Lab 13:
# 2.3: Raster Data:
# 2.3.1: Introduction to Raster:
# vignette("functions", package = "raster")
raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
new_raster = raster(raster_filepath)
new_raster
help("raster-package")

# 2.3.2: Basic Map Making:
plot(new_raster)
spplot(new_raster)
levelplot(new_raster)

# 2.3.3: Raster Classes:
raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
new_raster = raster(raster_filepath)

new_raster2 = raster(
  nrows = 6,
  ncols = 6,
  res = 0.5,
  xmn = -1.5,
  xmx = 1.5,
  ymn = -1.5,
  ymx = 1.5,
  vals = 1:36
)
multi_raster_file = system.file("raster/landsat.tif", package = "spDataLarge")
r_brick = brick(multi_raster_file)
r_brick
nlayers(r_brick)

raster_on_disk = raster(r_brick, layer = 1)
raster_in_memory = raster(
  xmn = 301905,
  xmx = 335745,
  ymn = 4111245,
  ymx = 4154085,
  res = 30
)
values(raster_in_memory) = sample(seq_len(ncell(raster_in_memory)))
crs(raster_in_memory) = crs(raster_on_disk)

r_stack = stack(raster_in_memory, raster_on_disk)
r_stack

# 2.4: Coordinate Reference Systems
# 2.4.1: Geographic Coordinate Systems:
st_proj_info(type = "datum")

# 2.4.2: Projected Coordinate Reference Systems:
st_proj_info(type = "proj")

# 2.4.3: CRS in R:
crs_data = rgdal::make_EPSG()
View(crs_data)

vector_filepath = system.file("vector/zion.gpkg", package = "spDataLarge")
new_vector = st_read(vector_filepath)

st_crs(new_vector) # get CRS
new_vector = st_set_crs(new_vector, 4326) # set CRS #Fly-on-projection only.
projection(new_raster) # get CRS

projection(new_raster) = "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0
                            +units=m +no_defs" # set CRS

# 2.5: Units:
luxembourgst_area(luxembourg)
st_area(luxembourg)
st_area(luxembourg) / 1000000
units::set_units(st_area(luxembourg), km ^ 2)
res(new_raster)
repr = projectRaster(new_raster, crs = "+init=epsg:26912")
res(repr)

# 2.6 Excercises:
# Q1.
summary(world)
# Geometry: Mulitpolygon.
#Number of Countries: 177
# Coordinate Reference System (CRS): Null.
crs(world) #Also gives NULL.
# Q2:
# Two differences:
# No lattitudinal and longitudial line referencing earth surface in map I made.
# Different shape and size of countries.

# Two Similarities:
# Population circles
# Same color
? plot # Cex: Set symbol size for scatter plots
# Cex Setting: Otherwise the circle disappers from the map. If only world$pop is used, it appears as small dots.
# Bonus: experiment with alternative method of visualization.
plot(world["continent"], reset = FALSE)
cex = sqrt(world$pop) / world$area_km2
world_cents = st_centroid(world, of_largest = TRUE)
plot(st_geometry(world_cents), add = TRUE, cex = cex)

# Q3: Can be improved.
plot(world["continent"], reset = FALSE)
cex = sqrt(world$pop) / 10000
world_cents = st_centroid(world, of_largest = TRUE)
plot(st_geometry(world_cents), add = TRUE, cex = cex)

Nigeria = world[world$name_long == "Nigeria",]
plot(
  st_geometry(Nigeria),
  expandBB = c(0, 2, 0.1, 2),
  col = "green",
  lwd = 2
)
# Challenge: Use text() and annotate map. (Can be improved)
text(
  x = 1,
  y = NULL,
  labels = "Map of Nigeria",
  adj = 1,
  pos = NULL,
  offset = 0.5,
  vfont = NULL,
  cex = 1,
  col = NULL,
  font = NULL
)

# Q4: Create a new raster.
my_raster = raster(
  nrow = 10,
  ncol = 10,
  xmn = -3,
  xmx = 3,
  ymn = -3,
  ymx = 3,
  vals = runif(n = 100, min = 0, max = 10)
)
plot(my_raster)

# Q5: Read raster/nlcd2011.tif
multi_raster_nlcd2011 = system.file("raster/nlcd2011.tif", package = "spDataLarge")
r_brick = brick(multi_raster_nlcd2011)
r_brick
nlayers(r_brick)

# Lab 14:
#Read Chapter 3 in the Geocomputation with R book. Complete the exercises at the end of Chapter 3 (Section 3.4). Submit your work here as an R script. You may type your written answers into the R script (use #) or turn in a separate Word doc if you prefer.

# Chapter 3: Attribute Data Operation:
# 3.1: Introduction:
# 3.2: Vector Attribute Manipulation:
methods(class = "sf") # methods for sf objects, first 12 shown
st_sf(data.frame(n = world$name_long), g = world$geom) #creates a geometry column named g.
dim(world) # it is a 2 dimensional object, with rows and columns
nrow(world) # how many rows?
ncol(world) # how many columns?
world_df = st_drop_geometry(world)
class(world_df)

# 3.2.1: Vector Attribute Subsetting:
world[1:6,] # subset rows by position
world[, 1:3] # subset columns by position
world[, c("name_long", "lifeExp")] # subset columns by name
# select(), filter(), pull(), slice() functions are used with dplyr package.
sel_area = world$area_km2 < 10000
summary(sel_area) # a logical vector
small_countries = world[sel_area,]
small_countries = world[world$area_km2 < 10000,]
small_countries = subset(world, area_km2 < 10000) #base package.

world1 = dplyr::select(world, name_long, pop)
names(world1)
world2 = dplyr::select(world, name_long:pop) # all columns between name_long and pop (inclusive)
world3 = dplyr::select(world,-subregion,-area_km2) # all columns except subregion and area_km2 (inclusive)
world4 = dplyr::select(world, name_long, population = pop)
names(world4)
world5 = world[, c("name_long", "pop")] # subset columns by name
names(world5)[names(world5) == "pop"] = "population" # rename column manually

d = data.frame(pop = 1:10, area = 1:10) # create throw-away data frame
d[, "pop", drop = FALSE] # equivalent to d["pop"] # return data frame object when selecting a single column
select(d, pop)
d[, "pop"] # return a vector when selecting a single column
pull(d, pop)

world[, "pop"] # data frame object
world$pop # vector objects
pull(world, pop)
world[, "pop"] # data frame object
world$pop # vector objects
pull(world, pop)
slice(world, 3:5)
world6 = filter(world, lifeExp > 82) # Countries with a life expectancy longer than 82 years

# # Symbol	Name
# # ==	Equal to
# # !=	Not equal to
# # >, <	Greater/Less than
# # >=, <=	Greater/Less than or equal
# # &, |, !	Logical operators: And, Or, Not

world7 = world %>%
  filter(continent == "Asia") %>%
  dplyr::select(name_long, continent) %>%
  slice(1:5) #dplyr pipe (%>%) operator.

world8 = slice(dplyr::select(filter(world, continent == "Asia"),
                             name_long, continent),
               1:5)

# 3.2.2: Vector Attribute Aggregation:
world_agg1 = aggregate(pop ~ continent,
                       FUN = sum,
                       data = world,
                       na.rm = TRUE)
class(world_agg1)

world_agg2 = aggregate(
  world["pop"],
  by = list(world$continent),
  FUN = sum,
  na.rm = TRUE
)
class(world_agg2)

world_agg3 = world %>%
  group_by(continent) %>%
  summarize(pop = sum(pop, na.rm = TRUE))

world %>%
  summarize(pop = sum(pop, na.rm = TRUE), n = n())

world %>%
  dplyr::select(pop, continent) %>%
  group_by(continent) %>%
  summarize(pop = sum(pop, na.rm = TRUE), n_countries = n()) %>%
  top_n(n = 3, wt = pop) %>%
  arrange(desc(pop)) %>%
  st_drop_geometry()
# vignette((package = "dplyr"))

# 3.2.3: Vector Attribute Joining
world_coffee = left_join(world, coffee_data)
class(world_coffee)
names(world_coffee)
plot(world_coffee["coffee_production_2017"])

coffee_renamed = rename(coffee_data, nm = name_long)
world_coffee2 = left_join(world, coffee_renamed, by = c(name_long = "nm"))
world_coffee_inner = inner_join(world, coffee_data)
nrow(world_coffee_inner)

setdiff(coffee_data$name_long, world$name_long)
str_subset(world$name_long, "Dem*.+Congo")

coffee_data$name_long[grepl("Congo,", coffee_data$name_long)] =
  str_subset(world$name_long, "Dem*.+Congo")
world_coffee_match = inner_join(world, coffee_data)
nrow(world_coffee_match)

coffee_world = left_join(coffee_data, world)
class(coffee_world)

# 3.2.4: Creating attributes and removing spatial information:
world_new = world # do not overwrite our original data
world_new$pop_dens = world_new$pop / world_new$area_km2

world %>%
  mutate(pop_dens = pop / area_km2)

world %>%
  transmute(pop_dens = pop / area_km2)

world_unite = world %>%
  unite("con_reg",
        continent:region_un,
        sep = ":",
        remove = TRUE)

world_separate = world_unite %>%
  separate(con_reg, c("continent", "region_un"), sep = ":")

world %>%
  rename(name = name_long)
new_names = c("i", "n", "c", "r", "s", "t", "a", "p", "l", "gP", "geom")

world %>%
  set_names(new_names)
world_data = world %>% st_drop_geometry()
class(world_data)

# 3.3 Manipulating raster objects:
elev = raster(
  nrows = 6,
  ncols = 6,
  res = 0.5,
  xmn = -1.5,
  xmx = 1.5,
  ymn = -1.5,
  ymx = 1.5,
  vals = 1:36
)
grain_order = c("clay", "silt", "sand")
grain_char = sample(grain_order, 36, replace = TRUE)
grain_fact = factor(grain_char, levels = grain_order)
grain = raster(
  nrows = 6,
  ncols = 6,
  res = 0.5,
  xmn = -1.5,
  xmx = 1.5,
  ymn = -1.5,
  ymx = 1.5,
  vals = grain_fact
)
levels(grain)[[1]] = cbind(levels(grain)[[1]], wetness = c("wet", "moist", "dry"))
levels(grain)
factorValues(grain, grain[c(1, 11, 35)])

# 3.3.1: Raser Subsetting:
elev[1, 1] # row 1, column 1
elev[1] # cell ID 1

r_stack = stack(elev, grain)
names(r_stack) = c("elev", "grain")
raster::subset(r_stack, "elev") # three ways to extract a layer of a stack
r_stack[["elev"]]
r_stack$elev

elev[1, 1] = 0
elev[]
elev[1, 1:2] = 0
elev[1, 1:2] = 0

# 3.3.2 : Summarizing Raster Objects:
cellStats(elev, sd)
hist(elev)

# 3.4. Ecercise:
library(spData)
data(us_states)
data(us_states_df)
# # Q1: Create a new object called us_states_name that contains only the NAME column from the us_states object. What is the class of the new object and what makes it geographic?
us_states_name = us_states[, c("NAME")] #subset by NAME
View(us_states_name) #View after subsetting data.
class(us_states_name) #Class = sf, data.frame.

us_states_name = us_states %>% dplyr::select(NAME)
class(us_states_name)

# # Q2: Select columns from the us_states object which contain population data. Obtain the same result using a different command (bonus: try to find three ways of obtaining the same result). Hint: try to use helper functions, such as contains or starts_with from dplyr (see ?contains).

# Five Different ways to extract same information:
us_states_pop1 = us_states[, c("total_pop_10", "total_pop_15")] #subset by total population 10 and 15.
us_states_pop2 = us_states[, 5:6]
us_states_pop3 = us_states %>% dplyr::select(total_pop_10, total_pop_15)
us_states_pop4 = us_states %>% dplyr::select(starts_with("total_pop"))
us_states_pop5 = us_states %>% dplyr::select(contains("total_pop"))

# # Q3: Find all states with the following characteristics (bonus find and plot them): 1) Belong to the Midwest region. 2) Belong to the West region, have an area below 250,000 km2 and in 2015 a population greater than 5,000,000 residents (hint: you may need to use the function units::set_units() or as.numeric()). 3) Belong to the South region, had an area larger than 150,000 km2 or a total population in 2015 larger than 7,000,000 residents.
#Midwest US:
us_states_widwest1 = us_states %>%
  filter(REGION == "Midwest")
plot(us_states_widwest1)
us_states_midwest2 = subset(us_states, REGION == "Midwest") #Use == for equal.
plot(us_states_midwest2)

identical(us_states_midwest$REGION, us_states_midwest1$REGION). = TRUE
identical(us_states_midwest, us_states_midwest1). # = FALSE => retains ID Number.

#West, Area and Population:
us_states_west1 = us_states %>% filter(REGION == "West",
                                       AREA < units::set_units(250000, km ^ 2),
                                       total_pop_15 > 5000000)
View(us_states_west1)

us_states_west2 = us_states %>% filter(REGION == "West",
                                       as.numeric(AREA) < 250000,
                                       total_pop_15 > 5000000)
View(us_states_west2)

#South, Area and Population:
us_states_south1 = us_states %>% filter(REGION == "South",
                                        AREA > units::set_units(150000, km ^ 2),
                                        total_pop_15 > 7000000)
View(us_states_south1)
us_states_south2 = us_states %>% filter(REGION == "South",
                                        as.numeric(AREA) > 150000,
                                        total_pop_15 > 7000000)
View(us_states_south2)

# # Q4: What was the total population in 2015 in the us_states dataset? What was the minimum and maximum total population in 2015?
us_states %>% summarize(
  total_pop = sum(total_pop_15),
  min_pop = min(total_pop_15),
  max_pop = max(total_pop_15)
)

# # Q5: How many states are there in each region?
us_states %>%
  group_by(REGION) %>%
  summarize(nr_of_states = n())

# # Q6: What was the minimum and maximum total population in 2015 in each region? What was the total population in 2015 in each region?
us_states %>%
  group_by(REGION) %>%
  summarize(
    min_pop = min(total_pop_15),
    max_pop = max(total_pop_15),
    tot_pop = sum(total_pop_15)
  )

# # Q7: Add variables from us_states_df to us_states, and create a new object called us_states_stats. What function did you use and why? Which variable is the key in both datasets? What is the class of the new object?
us_states_stats = us_states %>%
  left_join(us_states_df, by = c("NAME" = "state"))
class(us_states_stats)

# # Q8: us_states_df has two more rows than us_states. How can you find them? (hint: try to use the dplyr::anti_join() function)
us_states_df %>%
  anti_join(st_drop_geometry(us_states), by = c("state" = "NAME"))

# # Q9: What was the population density in 2015 in each state? What was the population density in 2010 in each state?
us_states2 = us_states %>%
  mutate(pop_dens_15 = total_pop_15 / AREA,
         pop_dens_10 = total_pop_10 / AREA)
View(us_states2)

# # Q10: How much has population density changed between 2010 and 2015 in each state? Calculate the change in percentages and map them
us_popdens_change = us_states2 %>%
  mutate(
    pop_dens_diff_10_15 = pop_dens_15 - pop_dens_10,
    pop_dens_diff_10_15p = (pop_dens_diff_10_15 / pop_dens_15) * 100
  )
plot(us_popdens_change["pop_dens_diff_10_15p"])

# # Q11: Change the columns' names in us_states to lowercase. (Hint: helper functions - tolower() and colnames() may help.)
us_states %>%
  setNames(tolower(colnames(.)))

# # Q12: Using us_states and us_states_df create a new object called us_states_sel. The new object should have only two variables - median_income_15 and geometry. Change the name of the median_income_15 column to Income
us_states_sel = us_states %>%
  left_join(us_states_df, by = c("NAME" = "state")) %>%
  dplyr::select(Income = median_income_15)
View(us_states_sel)

# # Q13: Calculate the change in median income between 2010 and 2015 for each state. Bonus: What was the minimum, average and maximum median income in 2015 for each region? What is the region with the largest increase of the median income?

us_income_change = us_states %>%
  left_join(us_states_df, by = c("NAME" = "state")) %>%
  mutate(income_change = median_income_15 - median_income_10)

us_income_change_reg = us_income_change %>%
  group_by(REGION) %>%
  summarize(
    min_income_change = min(income_change),
    mean_income_change = mean(income_change),
    max_income_change = max(income_change)
  )

us_income_change_reg %>%
  filter(mean_income_change == max(mean_income_change)) %>%
  pull(REGION) %>%
  as.character()

# # Q14: Create a raster from scratch with nine rows and columns and a resolution of 0.5 decimal degrees (WGS84). Fill it with random numbers. Extract the values of the four corner cells.
r = raster(nrow = 9, ncol = 9, res = 0.5, xmn = 0, xmx = 4.5,
           ymn = 0, ymx = 4.5, vals = rnorm(81))

# using cell IDs
r[c(1, 9, 81 - 9, 81)]

# using indexing
r[c(1, nrow(r)), c(1, ncol(r))]

# # Q15:What is the most common class of our example raster grain (hint: modal())?
grain_size = c("clay", "silt", "sand")
grain = raster(nrow = 6, ncol = 6, res = 0.5,
               xmn = -1.5, xmx = 1.5, ymn = -1.5, ymx = 1.5,
               vals = factor(sample(grain_size, 36, replace = TRUE),
                             levels = grain_size))
cellStats(grain, modal) %>%
  factorValues(grain, .)
factorValues(grain, modal(values(grain)))

# # Q16: Plot the histogram and the boxplot of the data(dem, package = "RQGIS") raster.
data(dem, package = "RQGIS")
par(mfrow = c(1, 2))
hist(dem)
boxplot(dem)
# 17: Now attach also data(ndvi, package = "RQGIS"). Create a raster stack using dem and ndvi, and make a pairs() plot.
data(ndvi, package = "RQGIS")
s = stack(dem, ndvi)
pairs(s)

# END OF CHAPTER 3 #

# Chapter 4: Spatial Data Operation:
# 4.1 Introduction
# 4.2 Spatial operations on vector data
# 4.2.1 Spatial subsetting
canterbury = nz %>% filter(Name == "Canterbury")
canterbury_height = nz_height[canterbury,]

nz_height[canterbury, , op = st_disjoint]

sel_sgbp = st_intersects(x = nz_height, y = canterbury)
class(sel_sgbp)
sel_logical = lengths(sel_sgbp) > 0
canterbury_height2 = nz_height[sel_logical,]

canterbury_height3 = nz_height %>%
  filter(st_intersects(x = ., y = canterbury, sparse = FALSE))

# 4.2.2 Topological relations
# create a polygon
a_poly = st_polygon(list(rbind(c(-1,-1), c(1,-1), c(1, 1), c(-1,-1))))
a = st_sfc(a_poly)

# create a line
l_line = st_linestring(x = matrix(c(-1,-1,-0.5, 1), ncol = 2))
l = st_sfc(l_line)

# create points
p_matrix = matrix(c(0.5, 1,-1, 0, 0, 1, 0.5, 1), ncol = 2)
p_multi = st_multipoint(x = p_matrix)
p = st_cast(st_sfc(p_multi), "POINT")

st_intersects(p, a)
st_intersects(p, a, sparse = FALSE)
st_disjoint(p, a, sparse = FALSE)[, 1]
st_touches(p, a, sparse = FALSE)[, 1]
st_within(p, a, sparse = FALSE)[, 1]
sel = st_is_within_distance(p, a, dist = 0.9) # can only return a sparse matrix
lengths(sel) > 0

# 4.2.3 Spatial joining
set.seed(2018) # set seed for reproducibility
(bb_world = st_bbox(world)) # the world's bounds
random_df = tibble(
  x = runif(n = 10, min = bb_world[1], max = bb_world[3]),
  y = runif(n = 10, min = bb_world[2], max = bb_world[4])
)
random_points = random_df %>%
  st_as_sf(coords = c("x", "y")) %>% # set coordinates
  st_set_crs(4326) # set geographic CRS

world_random = world[random_points,]
nrow(world_random)
random_joined = st_join(random_points, world["name_long"])

# 4.2.4 Non-overlapping joins
plot(st_geometry(cycle_hire), col = "blue")
plot(
  st_geometry(cycle_hire_osm),
  add = TRUE,
  pch = 3,
  col = "red"
)

any(st_touches(cycle_hire, cycle_hire_osm, sparse = FALSE))

cycle_hire_P = st_transform(cycle_hire, 27700)
cycle_hire_osm_P = st_transform(cycle_hire_osm, 27700)
sel = st_is_within_distance(cycle_hire_P, cycle_hire_osm_P, dist = 20)
summary(lengths(sel) > 0)

z = st_join(cycle_hire_P,
            cycle_hire_osm_P,
            join = st_is_within_distance,
            dist = 20)
nrow(cycle_hire)
nrow(z)

z = z %>%
  group_by(id) %>%
  summarize(capacity = mean(capacity))
nrow(z) == nrow(cycle_hire)

plot(cycle_hire_osm["capacity"])
plot(z["capacity"])

# 4.2.5 Spatial data aggregation
nz_avheight = aggregate(x = nz_height, by = nz, FUN = mean)

nz_avheight2 = nz %>%
  st_join(nz_height) %>%
  group_by(Name) %>%
  summarize(elevation = mean(elevation, na.rm = TRUE))

agg_aw = st_interpolate_aw(incongruent[, "value"], aggregating_zones,
                           extensive = TRUE)
agg_aw$value # show the aggregated result

# 4.2.6 Distance relations
nz_heighest = nz_height %>% top_n(n = 1, wt = elevation)
canterbury_centroid = st_centroid(canterbury)
st_distance(nz_heighest, canterbury_centroid)

co = filter(nz, grepl("Canter|Otag", Name))
st_distance(nz_height[1:3,], co)

plot(st_geometry(co)[2])
plot(st_geometry(nz_height)[2:3], add = TRUE)

# 4.3 Spatial operations on raster data
# 4.3.1 Spatial subsetting
id = cellFromXY(elev, xy = c(0.1, 0.1))
elev[id] # is the same as
raster::extract(elev, data.frame(x = 0.1, y = 0.1))

clip = raster(
  xmn = 0.9,
  xmx = 1.8,
  ymn = -0.45,
  ymx = 0.45,
  res = 0.3,
  vals = rep(1, 9)
)
elev[clip]
extract(elev, extent(clip)) # we can also use extract

elev[1:2, drop = FALSE]    # spatial subsetting with cell IDs
elev[1, 1:2, drop = FALSE] # spatial subsetting by row,column indices

# create raster mask
rmask = elev
values(rmask) = sample(c(NA, TRUE), 36, replace = TRUE)

# spatial subsetting
elev[rmask, drop = FALSE]           # with [ operator
mask(elev, rmask)                   # with mask()
overlay(elev, rmask, fun = "max")   # with overlay

# 4.3.2 Map algebra:
# vignette("Raster")

# 4.3.3 Local operations
rcl = matrix(c(0, 12, 1, 12, 24, 2, 24, 36, 3), ncol = 3, byrow = TRUE)
recl = reclassify(elev, rcl = rcl)

elev + elev
elev ^ 2
log(elev)
elev > 5

# 4.3.4 Focal operations
r_focal = focal(elev, w = matrix(1, nrow = 3, ncol = 3), fun = min)

# 4.3.5 Zonal operations
z = zonal(elev, grain, fun = "mean") %>%
  as.data.frame()
z

# 4.3.6 Global operations and distances
# 4.3.7 Map algebra counterparts in vector processing
# 4.3.8 Merging rasters
aut = getData("alt", country = "AUT", mask = TRUE)
ch = getData("alt", country = "CHE", mask = TRUE)
aut_ch = merge(aut, ch)

# 4.3.7 Map algebra counterparts in vector processing
# 4.3.8 Merging rasters
aut = getData("alt", country = "AUT", mask = TRUE)
ch = getData("alt", country = "CHE", mask = TRUE)
aut_ch = merge(aut, ch)

# 4.4 Exercises:
# Q1: It was established in Section 4.2 that Canterbury was the region of New Zealand containing most of the 100 highest points in the country. How many of these high points does the Canterbury region contain?

library(tmap)
tmap_mode("view")
# tmap mode set to interactive viewing
qtm(nz) + qtm(nz_height)

canterbury = nz %>% filter(Name == "Canterbury")
canterbury_height = nz_height[canterbury, ]
nrow(canterbury_height) # answer: 70

# Q2:Which region has the second highest number of nz_height points in, and how many does it have?
nz_height_count = aggregate(nz_height, nz, length)
nz_height_combined = cbind(nz, count = nz_height_count$elevation)
nz_height_combined %>%
  st_set_geometry(NULL) %>%
  dplyr::select(Name, count) %>%
  arrange(desc(count)) %>%
  slice(2)

# Q3:Generalizing the question to all regions: how many of New Zealand's 16 regions contain points which belong to the top 100 highest points in the country? Which regions? Bonus: create a table listing these regions in order of the number of points and their name.
nz_height_count = aggregate(nz_height, nz, length)
nz_height_combined = cbind(nz, count = nz_height_count$elevation)
nz_height_combined %>%
  st_set_geometry(NULL) %>%
  dplyr::select(Name, count) %>%
  arrange(desc(count)) %>%
  na.omit()

# Q4:Use data(dem, package = "RQGIS"), and reclassify the elevation in three classes: low, medium and high. Secondly, attach the NDVI raster (data(ndvi, package = "RQGIS")) and compute the mean NDVI and the mean elevation for each altitudinal class.
library(classInt)
data(dem, package = "RQGIS")
data(ndvi, package = "RQGIS")
summary(dem)

brk = classIntervals(values(dem), n = 3)$brk
# also try
# breask = classIntervals(values(dem), n = 3, style = "fisher")
# construct reclassification matrix

rcl = matrix(c(brk[1] - 1, brk[2], 1, brk[2], brk[3], 2, brk[3], brk[4], 3),
             ncol = 3, byrow = TRUE)

# reclassify altitudinal raster
recl = reclassify(dem, rcl = rcl)

# compute the mean dem and ndvi values for each class
zonal(stack(dem, ndvi), recl, fun = "mean")

# Q5: Apply a line detection filter to raster(system.file("external/rlogo.grd", package = "raster")). Plot the result. Hint: Read ?raster::focal().
# from the focal help page (?raster::focal()):
# Laplacian filter: filter=matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)
# Sobel filter: filter=matrix(c(1,2,1,0,0,0,-1,-2,-1) / 4, nrow=3)

# just retrieve the first channel of the R logo
r = raster(system.file("external/rlogo.grd", package = "raster"))
# compute the Sobel filter
filter = matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1) / 4, nrow = 3)
sobel = focal(r, w = filter)
plot(sobel, col = c("black", "white"))

# Q6: Calculate the NDVI of a Landsat image. Use the Landsat image provided by the spDataLarge package (system.file("raster/landsat.tif", package="spDataLarge")).
file = system.file("raster/landsat.tif", package = "spDataLarge")
r = stack(file)
# compute NDVI manually
ndvi = (r[["landsat.4"]] - r[["landsat.3"]]) / (r[["landsat.4"]] + r[["landsat.3"]])

# compute NDVI with the help of RStoolbox
ndvi_rstoolbox = RStoolbox::spectralIndices(r, red = 3, nir = 4, indices = "NDVI")
all.equal(ndvi, ndvi_rstoolbox)

# compute NDVI with the help of RStoolbox
ndvi_rstoolbox = RStoolbox::spectralIndices(r, red = 3, nir = 4, indices = "NDVI")
all.equal(ndvi, ndvi_rstoolbox)

# compute NDVI with the help of RStoolbox
ndvi_rstoolbox = RStoolbox::spectralIndices(r, red = 3, nir = 4, indices = "NDVI")
all.equal(ndvi, ndvi_rstoolbox)

# Q7: A StackOverflow post shows how to compute distances to the nearest coastline using raster::distance(). Retrieve a digital elevation model of Spain, and compute a raster which represents distances to the coast across the country (hint: use getData()). Second, use a simple approach to weight the distance raster with elevation (other weighting approaches are possible, include flow direction and steepness); every 100 altitudinal meters should increase the distance to the coast by 10 km. Finally, compute the difference between the raster using the Euclidean distance and the raster weighted by elevation. Note: it may be wise to increase the cell size of the input raster to reduce compute time during this operation.
# find out the ISO_3 code of Spain
dplyr::filter(ccodes(), NAME %in% "Spain")

# retrieve a dem of Spain
dem = getData("alt", country = "ESP", mask = FALSE)

# change the resolution to decrease computing time
agg = aggregate(dem, fact = 5)

# poly = getData("GADM", country = "ESP", level = 1) # fails
plot(dem)
plot(poly, add = TRUE)

# visualize NAs
plot(is.na(agg))

# construct a distance input raster
# we have to set the land cells to NA and the sea cells to an arbitrary value since 
# raster::distance computes the distance to the nearest non-NA cell
dist = is.na(agg)
cellStats(dist, summary)

# convert land cells into NAs and sea cells into 1s
dist[dist == FALSE] = NA
dist[dist == TRUE] = 1
plot(dist)

# compute distance to nearest non-NA cell
dist = raster::distance(dist)

# just keep Spain
dist = mask(dist, poly)

# convert distance into km
dist = dist / 1000

# now let's weight each 100 altitudinal meters by an additionaly distance of 10 km
agg = mask(agg, poly)
agg[agg < 0] = 0
weight = dist + agg / 100 * 10
plot(weight - dist)

#END OF CHAPTER 4 #
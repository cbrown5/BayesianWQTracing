#Calculate run-off and sediment yield
#
#CJ Brown 28 April 2016
#
# v2: uses catchment averages to calculate runoff
# Catchment averages are almost identical to pixel by pixel (r2 = 0.9925)
#

library(dplyr)
library(raster)
library(rgdal)
library(sp)

setwd('/Users/s2973410/Databases/Fiji reefs analysis/Fiji catchments/Catchments_processed_v2')

#
# PARAMETERS
#
# Using Estimates of sediment yield from Table 1 of Neil et al. 2002, mg L-1 for Northern Monsoon catchments
yield_natural <- 32/1000 #convert 32 mg to grams per litre
yield_degrade <- 99/1000

#
# Load catchments and data layers
#
fijiproj <- "+proj=utm +zone=60 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

rrain <- raster('RainfallVanuaLevu_wetseason.grd')
# dat <- read.csv('River_mouths_Basin_cover.csv', header = T)
dat <- read.csv('River_mouths_Basin_cover_densevegonly_12Apr2016.csv', header = T)
rcatch <- raster('VanuaLevuBasins.grd')
projection(rcatch) <- fijiproj


#
# Land cover grid
#
fijiproj <- "+proj=utm +zone=60 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
rland <- raster('VanuaLevu_Land0_Mangroves1_Water2.grd')


#
# Functions
#

#Function to calculate proportional run-off from Zhang et al. 2001
prop_runoff <- function(f, P){
	#Where
	# f is proportion forested
	# P is precipitation (mm)
	#
	#Returns proportion of rainfall that runs off
	#
	
	#Parameters from Zhang
	wforest <- 2
	E0forest <- 1410

	wgrass <- 0.5
	E0grass <- 1100
	
	#calculate Evapotransipriation
	forest_ET <- P*(1 + wforest*(E0forest/P)) / (1 + (wforest*(E0forest/P)) + (P/E0forest))
	grass_ET <- P*(1 + wgrass*(E0grass/P)) / (1 + (wgrass*(E0grass/P)) + (P/E0grass))
	
	total_ET <- (f * forest_ET) + ((1-f) * grass_ET)
	
	#Calculate proportion run-off
	totalro <- (P - total_ET)/P
	totalro	
	}

# p <- seq(0, 9600, length.out = 100)
# plot(p, prop_runoff(0, p))
# lines(p, prop_runoff(1, p))


#
# Cell areas
#

#Cell areas
cellarea <- prod(res(rcatch))/10000 #hectares
sqrm_perpix <- prod(res(rcatch)) #square metres per pixel,

#
# Calculate catchment average precipitation and forest cover
#
rain.tab <- zonal(rrain, rcatch, fun = 'mean', na.rm = T)
raindf <- data.frame(basin = rain.tab[,1], av_rain_mm = rain.tab[,2])

datout <- dat %>% inner_join(raindf)

#Rainfall in litres - rainfall_mm * [basinarea_ha * (m2/ha)] * proportion forest cover
#1 mm rainfall = 1L/m2

datout$rain_L_forest <- datout$av_rain_mm * datout$basin_area_ha * 10000 *datout$forest_cover 
datout$rain_L_uv <- datout$av_rain_mm * datout$basin_area_ha * 10000 * (1 - datout$forest_cover)

datout$annual_rainGL <- (datout$rain_L_forest + datout$rain_L_uv)/1000000000

#
# Calculate runoff and sediment yield
#

#proportional runoff for forest and unvegetated areas in each catchment
datout$propro_forest <- prop_runoff(1, datout$av_rain_mm)
datout$propro_uf <- prop_runoff(0, datout$av_rain_mm)

#Runoff in litres
datout$runoff_forest <- datout$propro_forest * datout$rain_L_forest
datout$runoff_uf <- datout$propro_uf * datout$rain_L_uv

datout$annual_runoffGL <- (datout$runoff_forest + datout$runoff_uf)/1000000000

#Sediment yield in tons. Divide by 1000^2 to convert grams into tons
datout$sedyld_tons <- ((datout$runoff_forest * yield_natural) + 
(datout$runoff_uf * yield_degrade)) / (1000^2)

with(datout, plot(basin_area_ha, sedyld_tons))
with(datout, plot(annual_rainGL, annual_runoffGL))

#
# Calculate sediment yield per km2 for comparison with Olley et al. 2015
#
plot(log10(datout$av_rain_mm*datout$propro_forest), log10(datout$sedyld_tons/(datout$basin_area_ha*0.01)))
#similar scales. 

#
# Select main columns and save
#
datsave <- datout %>% dplyr::select(-rain_L_forest, -rain_L_uv, -runoff_forest, -runoff_uf)

write.csv(datsave, 'Rivers_Basins_datv4.csv',row.names = F)

# #Compare to pixel by pixel calculations
# dat2 <- read.csv('Rivers_Basins_datv3.csv')

# plot(dat2$sedyld_tons, datsave$sedyld_tons)
# abline(lm(dat2$sedyld_tons~ datsave$sedyld_tons + 0))
# summary(lm(dat2$sedyld_tons~ datsave$sedyld_tons + 0))






################################################################################
# Project: Strategic surveillance of zoonotic pathogens
# Author: Marlon E. Cobos
# Date: 10/07/2024 (dd/mm/yyyy)
################################################################################



# Description ------------------------------------------------------------------
# The script contains code and instructions to load data, run all analyses, and 
# produce figures for the project.
#
# Part of the data can be obtained from online databases using the code or 
# following instructions in this script. Other data can be obtained from a 
# Figshare repository (see below).
# 
# Instructions to set working directory and get ready for analysis:
# 1. Download the data and scripts provided in the repository 
#    https://github.com/marlonecobos/Hanta_PAN_suitability/tree/main/Sampling_pathogens
# 2. Save the folders downloaded in your directory for work (unzip if needed).
# 3. Download manually environmental layers from the ChelsaClimate database:
#    https://chelsa-climate.org/; v2; climatologies/1981-2010/bio
#       Variables: bio1, bio5, bio6, bio7, bio12, bio15, bio16, bio17, pet_max, 
#                  pet_mean, pet_min, vpd_max, vpd_mean, vpd_min
# 4. Save downloaded variables in the folder "Data/chelsa" withing your directory
#    for work.
# 5. Define your working directory.
# 6. Run analyses according to comments and instructions.
#
# How to load points and layers to explore in Google Earth:
#   The MOP raster layer to load has been provided
#   Loading raster layer:
#     1. In the Google Earth menu: "File", click in "Open"
#     2. Select "GeoTIFF (*.tif)" as file format (close to the button "Open")
#     3. Select the file "rec_mop_distance1.tif"; click Open
#     4. Modify transparency to be able to check other geographic features
#     5. Click OK
#
#   The csv files with points are generated with this script.
#   Loading points:
#     1. In the Google Earth menu: "File", click in "Import"
#     2. Select "Generic Text (*.txt *.csv)" as file format (close to "Open")
#     3. Select the file needed (e.g., "summary_per_locality.csv"); click Open
#     4. In the section "Delimited", select "Comma"; Click Next
#     5. Define "x" and "y" as the Longitude and Latitude fields, respectively
#     6. Click Next
#     7. Define the type of variable for each field (e.g., "cell" is integer, 
#        but "prev" is floating point)
#     8. Click Finish
#     9. Select "Yes", when asked if you want to apply a style template
#     10. Select "Create new template"
#     11. Set "cell" as the name field.; Click OK
#     12. Save it with a name of your preference.
# ------------------------------------------------------------------------------


# Working directory ------------------------------------------------------------
setwd("YOUR/DIRECTORY")  # change and/or used as needed
# ------------------------------------------------------------------------------



# Packages ---------------------------------------------------------------------
# install
# install.packages("geodata")
# install.packages("mop")
# install.packages("remotes")
# 
# remotes::install_github("claununez/biosurvey")

# load
library(geodata)  # also loads terra
library(biosurvey) 
library(mop)
library(ks)
# ------------------------------------------------------------------------------



# Data -------------------------------------------------------------------------
# world layer
wrld <- world(path = "Data/spatial")

# Panama country border
pan <- gadm("PAN", level = 0, path = "Data/spatial")

# Panama provinces
panpro <- gadm("PAN", level = 1, path = "Data/spatial")

# road buffer 2.5 km
roads2 <- vect("Data/spatial/road_buffer_2.5km_un_dis_clip.gpkg")  # provided

# variables (downloaded previously, see "Description" above)
vars <- rast(list.files("Data/chelsa", pattern = ".tif$", full.names = TRUE))

# records from Arctos
rec <- read.csv("Data/records/Panama orthohanta positives_11 oct 2023.csv")
# ------------------------------------------------------------------------------



# Data preparation -------------------------------------------------------------
# project layers
roads2 <- project(roads2, pan)
  
# mask layers
## variables to PAN
varsp <- crop(vars, pan, mask = TRUE)
varsp <- match_na_raster(varsp)

## variables to PAN roads
varsr2 <- crop(vars, roads2, mask = TRUE)


# save masked layers
dir.create("Data/E_variables")

writeRaster(varsp, filename = "Data/E_variables/panama_var.tif")
writeRaster(varsr2, filename = "Data/E_variables/accessible_var.tif")


# filter localities to only the ones that are more accessible based on roads
## get values of variables in points
recev <- extract(varsr2, as.matrix(rec[, c("DEC_LONG", "DEC_LAT")]))

## filter records and extracted values
nona <- !is.na(recev$bio1)
recev <- recev[nona, ]
recr <- rec[nona, ]
# ------------------------------------------------------------------------------



# Explorations of sampling effort ----------------------------------------------
# records from Arctos
x11()

plot(pan, col = "gray85", main = "All records from Arctos")
points(recr[, c("DEC_LONG", "DEC_LAT")], pch = 3, 
       col = biosurvey:::make_alpha("blue", 0.3))


# analyses of sampling effort per pixel (1 km2)
# from now on, pixels are our candidate localities

## base raster layer for computation
denl <- varsp$bio1
denl <- (denl < 1000) * 1

## number of independent sampling days per pixel
### unique sampling days per locality
sdays <- recr[!duplicated(paste(recr$BEGAN_DATE, recr$DEC_LONG, recr$DEC_LAT)), ]

### calculations using a raster layer as a base for counting
sdv <- vect(sdays, geom = c("DEC_LONG", "DEC_LAT"), crs = crs(vars))
denlsd <- rasterize(sdv, denl, fun = function(i) {length(i)}, background = NA)
sdays_loc <- as.data.frame(denlsd, xy = TRUE, cells = TRUE)

## number of species per pixel
### unique species names per locality
sploc <- recr[!duplicated(paste(recr$SCIENTIFIC_NAME, recr$DEC_LONG, 
                                recr$DEC_LAT)), ]

### calculations using a raster layer as a base for counting
splv <- vect(sploc, geom = c("DEC_LONG", "DEC_LAT"), crs = crs(vars))
denlspl <- rasterize(splv, denl, fun = function(i) {length(i)}, background = NA)
denlspl <- rasterize(splv, denl, field = "SCIENTIFIC_NAME", 
                     fun = function(i) {length(unique(i))},
                      background = NA)
spl_loc <- as.data.frame(denlspl, xy = TRUE, cells = TRUE)

## number of samples per pixel
### counting number of records (samples) using a raster layer as a base
recv <- vect(recr, geom = c("DEC_LONG", "DEC_LAT"), crs = crs(vars))
denlr <- rasterize(recv, denl, fun = function(i) {length(i)}, background = NA)
rec_loc <- as.data.frame(denlr, xy = TRUE, cells = TRUE)

## hantavirus prevalence per pixel
### calculation of prevalences based on Arctos data
recrp <- recr
recrp$LOC <- paste(recr$DEC_LONG, recr$DEC_LAT)
recrp$TEST <- ifelse(recrp$DETECTED == "virus: Orthohantavirus", 1, 0)

recrp <- as.matrix(table(recrp$LOC, recrp$TEST))

### prepare data for computation in a raster layer
prevxy <- do.call(rbind, strsplit(rownames(recrp), " "))

recrpre <- data.frame(x = as.numeric(prevxy[, 1]), y = as.numeric(prevxy[, 2]), 
                      negative = recrp[, 1], positive = recrp[, 2], 
                      prevalence = recrp[, 2] / (recrp[, 1] + recrp[, 2]), 
                      row.names = NULL)

### calculating average prevalence per locality using a raster layer as a base
prevlv <- vect(recrpre, geom = c("x", "y"), crs = crs(vars))
denlprev <- rasterize(prevlv, denl, field = "prevalence", fun = mean,
                      background = NA)
prev_loc <- as.data.frame(denlprev, xy = TRUE, cells = TRUE)
prev_loc$meanr <- round(prev_loc$mean, 2)


# plot represent localities according to numbers
x11()
par(mfrow = c(4, 1), cex = 0.8)

plot(pan, col = "gray85", main = "Number of sampling days per locality")
points(sdays_loc[sdays_loc$V1 >= 5, 2:3])
points(sdays_loc[sdays_loc$V1 >= 45, 2:3], col = "red", pch = 19)
points(sdays_loc[sdays_loc$V1 >= 10, 2:3], pch = 3)
legend(-83, 8, legend = c(">5", ">10", ">45"), pch = c(1, 3, 19), 
       col = c(1, 1, 2), bty = "n", cex = 0.8)

plot(pan, col = "gray85", main = "Number of species per locality")
points(spl_loc[spl_loc$V1 >= 2, 2:3])
points(spl_loc[spl_loc$V1 >= 8, 2:3], col = "red", pch = 19)
points(spl_loc[spl_loc$V1 >= 5, 2:3], pch = 3)
legend(-83, 8, legend = c(">2", ">5", ">8"), pch = c(1, 3, 19), 
       col = c(1, 1, 2), bty = "n", cex = 0.8)

plot(pan, col = "gray85", main = "Number of samples per locality")
points(rec_loc[rec_loc$V1 >= 50, 2:3])
points(rec_loc[rec_loc$V1 >= 1000, 2:3], col = "red", pch = 19)
points(rec_loc[rec_loc$V1 >= 200, 2:3], pch = 3)
legend(-83, 8, legend = c(">50", ">200", ">1000"), pch = c(1, 3, 19), 
       col = c(1, 1, 2), bty = "n", cex = 0.8)

plot(pan, col = "gray85", main = "Hantavirus prevalence per locality")
points(prev_loc[prev_loc$meanr >= 0.1, 2:3])
points(prev_loc[prev_loc$meanr >= 0.5, 2:3], col = "red", pch = 19)
points(prev_loc[prev_loc$meanr >= 0.2, 2:3], pch = 3)
legend(-83, 8, legend = c(">0.10", ">0.20", ">0.50"), pch = c(1, 3, 19), 
       col = c(1, 1, 2), bty = "n", cex = 0.8)

## save results
dir.create("Results")

write.csv(sdays_loc, "Results/ndays_per_locality.csv", row.names = FALSE)
write.csv(spl_loc, "Results/nspp_per_locality.csv", row.names = FALSE)
write.csv(rec_loc, "Results/samples_per_locality.csv", row.names = FALSE)
write.csv(prev_loc, "Results/prev_per_locality.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# Prepare biosurvey objects for selection of sampling sites --------------------
# master matrix
pan_mm <- prepare_master_matrix(region = pan, mask = roads2, variables = varsp, 
                                do_pca = TRUE, center = TRUE, scale = TRUE)


## blocks in master matrix
pan_mm <- make_blocks(master_matrix = pan_mm, variable_1 = "PC1", 
                      variable_2 = "PC2", n_cols = 30, n_rows = 30,
                      block_type = "equal_area")


save_master(pan_mm, file_name = "Results/master_matrix_accessible.rds")
#pan_mm <- read_master("Results/master_matrix_accessible.rds") 

sums <- summary(pan_mm$PCA_results)

write.csv(sums$importance, file = "Results/summary_PCA.csv", row.names = TRUE)
write.csv(sums$rotation, file = "Results/loadings_PCA.csv", row.names = TRUE)


# check master matrix object
x11()
explore_data_EG(pan_mm, variable_1 = "PC1", variable_2 = "PC2")
# ------------------------------------------------------------------------------



# historically sampled points in geography and environment ---------------------
# transform raw variables to PCs (most of following analyses rely on PCs)
varspc <- predict(varsp, pan_mm$PCA_results, filename = "Results/pan_pcs.tif")
varsr2pc <- predict(varsr2, pan_mm$PCA_results, 
                    filename = "Results/accessible_pcs.tif")

recevpc <- predict(pan_mm$PCA_results, recev)
write.csv(recevpc, file = "Results/hist_localities_pcs.csv", row.names = FALSE)


# variables to be plotted
pcs <- c("PC1", "PC2")
xy <- c("x", "y")
lola <- c("Longitude", "Latitude")


# combine data to explore 
# explorations in environmental and geographic space where done graphically in R
# interactive explorations in geography were done in Google Earth (import csv)

## extract environments in sampled localities
rec_loce_pc <- extract(varsr2pc, rec_loc[, 2:3])[, 2:3]

## put information together
rec_loce_pc1 <- cbind(rec_loc[, 1:3], rec_loce_pc, days = sdays_loc$V1, 
                      species = spl_loc$SCIENTIFIC_NAME, samples = rec_loc$V1,
                      prev = prev_loc$mean)

## save in file
write.csv(rec_loce_pc1, "Results/summary_per_locality.csv", row.names = FALSE)


## plotting
x11()
par(mfrow = c(2, 1), cex = 0.6)

plot(pan, col = "gray85", main = "Localities sampled last 20 years")
points(rec_loce_pc1[, xy], pch = 3, col = "#0078FF")

par(mar = c(4, 4, 1, 1))
plot(as.data.frame(varspc)[, 1:2], pch = 19, col = "gray85")
points(rec_loce_pc1[, pcs], pch = 3, col = "#0078FF")
# ------------------------------------------------------------------------------



# MOP analyses -----------------------------------------------------------------
## records vs PAN
mop_rpan <- mop(m = recevpc[, 1:2], g = varspc[[1:2]], type = "simple", 
                calculate_distance = TRUE, where_distance = "all", 
                fix_NA = FALSE, rescale_distance = TRUE, parallel = TRUE, 
                n_cores = 16)

## roads vs PAN
mop_roadpan <- mop(m = varsr2pc[[1:2]], g = varspc[[1:2]], type = "simple", 
                   calculate_distance = TRUE, where_distance = "all", 
                   fix_NA = FALSE, rescale_distance = TRUE, parallel = TRUE, 
                   n_cores = 16)

## save results
writeRaster(mop_rpan$mop_basic, filename = "Results/rec_mop_basic.tif")
writeRaster(mop_rpan$mop_simple, filename = "Results/rec_mop_simple.tif")
writeRaster(mop_rpan$mop_distances, filename = "Results/rec_mop_distance.tif")
writeRaster(mop_roadpan$mop_basic, filename = "Results/road_mop_basic.tif")
writeRaster(mop_roadpan$mop_simple, filename = "Results/road_mop_simple.tif")
writeRaster(mop_roadpan$mop_distances, 
            filename = "Results/road_mop_distance.tif")


# plot MOP results
x11()
par(mfrow = c(2, 1), cex = 0.6)

## sampled localities vs PAN
plot(mop_rpan$mop_distances, col = darkros(length(levels(mopd))), 
     main = "Environmental comparison: sampled vs Panama")
points(rec_loce_pc1[, xy], pch = 3, col = "#0078FF")

## accessible areas vs PAN
plot(mop_roadpan$mop_distances, col = darkros(length(levels(mopdr))), 
     main = "Environmental comparison: accessible vs Panama")
plot(roads2, border = "gray", add = TRUE, 
     col = biosurvey:::make_alpha("gray", 0.3))
# ------------------------------------------------------------------------------



# selecting points based on historic sampling and other considerations ---------
# accessible areas in environmental conditions
x11()
plot(pan_mm$data_matrix[, pcs], pch = 19)


# historically sampled points
points(rec_loce_pc1[, pcs], pch = 3, col = "#0078FF")


# places selected based on sampling effort, and biogeographic and environmental 
# considerations
selected <- c(188029, 64043, 138570) 

rec_loce_pc1[rec_loce_pc1$cell %in% selected, ]

## check points 
points(rec_loce_pc1[rec_loce_pc1$cell %in% selected, pcs], pch = 16, 
       col = "red")


# places selected to complete environmental coverage
## check points in certain environmental ranges
rec_loce_pc1[rec_loce_pc1$PC1 > 2 & rec_loce_pc1$PC1 < 3.5 &
              rec_loce_pc1$PC2 > -3 & rec_loce_pc1$PC2 < -2, ]

rec_loce_pc1[rec_loce_pc1$PC1 > -1 & rec_loce_pc1$PC1 < 1 &
              rec_loce_pc1$PC2 > 3 & rec_loce_pc1$PC2 < 4, ]

rec_loce_pc1[rec_loce_pc1$PC1 > -1 & rec_loce_pc1$PC1 < 1 &
               rec_loce_pc1$PC2 > 3 & rec_loce_pc1$PC2 < 4, ]

## check these points in Google Earth 

## after checks and considering the number of sampling days per locality
ecomp <- c(163082, 29940, 84193)

rec_loce_pc1[rec_loce_pc1$cell %in% ecomp, ]

## check points 
x11()
### previous points
plot(pan_mm$data_matrix[, pcs], pch = 19)
points(rec_loce_pc1[, pcs], pch = 3, col = "#0078FF")
points(rec_loce_pc1[rec_loce_pc1$cell %in% selected, pcs], pch = 16, 
       col = "red")

## newly selected
points(rec_loce_pc1[rec_loce_pc1$cell %in% ecomp, pcs], pch = 16, 
       col = "yellow")


# these are the first set of points selected (6 sites in all Panama)
# all this points have been sampled before and are relatively easy to access
six_sites <- rec_loce_pc1[rec_loce_pc1$cell %in% c(selected, ecomp), ]

## save six sites
colnames(six_sites)[2:3] <- lola

write.csv(six_sites, "Results/six_selected_survey_sites.csv", row.names = FALSE)


# check 2nd step selected points in environment and in geography
x11()
par(mfrow = c(2, 1), cex = 0.6)

plot(roads2, col = "gray1")
points(six_sites[, lola], pch = 19, col = "red")

par(mar = c(4, 4, 1, 1))

plot(pan_mm$data_matrix[, pcs], pch = 19)
points(six_sites[, pcs], pch = 19, col = "red")
# ------------------------------------------------------------------------------



# select four more points to complete E coverage -------------------------------
# selecting points from poll of sampled localities

## exploring in environmental space using blocks, six sites, and localities
## to aid, we used Google Earth to check the points in geography
### colors for blocks
nblocks <- length(unique(pan_mm$data_matrix$Block))
col_blocks <- sample(gray.colors(nblocks), nblocks)

### plot to identify crucial points
x11()
plot_blocks_E(pan_mm, block_ID = FALSE, 
              col_all = col_blocks[as.factor(pan_mm$data_matrix$Block)])
points(rec_loce_pc1[, pcs], pch = 3, col = "blue")
points(six_sites[, pcs], pch = 19, col = "red")
identify(rec_loce_pc1[, pcs], labels = rec_loce_pc1$cell)  # see next comment

## click on the points to check which cell number they have, check those in 
## Google Earth to understand other geographic characteristics of the locality

## after explorations, these are the new sites to be included
more_ecomp <- c(59432, 28928, 109606)
more_ecomp <- rec_loce_pc1$cell %in% more_ecomp


# check 3rd step selected points in environment and in geography
x11()
par(mfrow = c(2, 1), cex = 0.6)

## geography
plot(roads2, col = "gray1")
points(six_sites[, lola], pch = 19, col = "red")

points(rec_loce_pc1[more_ecomp, xy], pch = 19, col = "green")


par(mar = c(4, 4, 1, 1))

## environment
plot(pan_mm$data_matrix[, pcs], pch = 19)
points(six_sites[, pcs], pch = 19, col = "red")

points(rec_loce_pc1[more_ecomp, pcs], pch = 19, col = "green")


# there is one big gap in sampling considering environmental conditions (top).
# none of the sampled localities goes close to that environmental region.

# to select select one more point based on environmental blocks from a 
# non-sampled E region
## check E blocks 
x11()
plot_blocks_E(pan_mm, block_ID = TRUE, col_ID = "blue", cex_ID = 0.5,
              col_all = col_blocks[as.factor(pan_mm$data_matrix$Block)])
points(six_sites[, pcs], pch = 19, col = "red")
points(rec_loce_pc1[more_ecomp, pcs], pch = 19, col = "green")

## export points from blocks of interest
int_blocks <- c(522, 523, 554, 553, 552, 584, 583, 614, 613, 645, 644)
write.csv(pan_mm$data_matrix[pan_mm$data_matrix$Block %in% int_blocks, 
                             c(1:2, 17:19)],
          "Results/blocks_interest.csv", row.names = TRUE)

## check the MOP_for_GEarth layer: sampled localities vs Panama in Google Earth
## also in Google Earth, check the points in the blocks that are of interest
## use the unnamed field as name

## as expected all these points are in yellow areas of the MOP layer, meaning
## they are dissimilar to previously sampled localities

## after exploring other geographic features: how close roads are, and distance 
## human settlements
## and the position of the block they belong to, we picked the point with the 
## cell ID "75354"


## select a locality based on novelty of environments represented and 
## accessibility from a road and a populated area 
sel_ce <- which(rownames(pan_mm$data_matrix) == "75354")


## check
x11()
par(mfrow = c(2, 1), cex = 0.7)

plot(roads2, col = "gray1")
points(six_sites[, lola], pch = 19, col = "red")
points(rec_loce_pc1[more_ecomp, xy], pch = 19, col = "green")
points(pan_mm$data_matrix[sel_ce, lola], pch = 19, col = "yellow")


par(mar = c(4, 4, 1, 1))

plot(pan_mm$data_matrix[, pcs], pch = 19)
points(six_sites[, pcs], pch = 19, col = "red")
points(rec_loce_pc1[more_ecomp, pcs], pch = 19, col = "green")
points(pan_mm$data_matrix[sel_ce, pcs], pch = 19, col = "yellow")


# save ten points selected
new_site <- cbind(cell = 75354, pan_mm$data_matrix[sel_ce, c(lola, pcs)],
                  days = NA, species = NA, samples = NA, prev = NA)

three_sites <- rec_loce_pc1[more_ecomp, ]
colnames(three_sites)[2:3] <- lola

ten_sites <- rbind(six_sites, three_sites, new_site)

## save six sites
write.csv(ten_sites, "Results/ten_selected_survey_sites.csv", row.names = FALSE)
# ------------------------------------------------------------------------------

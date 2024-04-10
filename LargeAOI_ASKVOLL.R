
library(raster)
library(rgdal)
library(sp)
library(lidR)
getwd()
las <- readLAS('AOI - Copy.laz')
plot(las)

LASfile <- system.file("extdata", "AOI - Copy.laz", package="lidR")
#las <- readLAS(las, select = "xyzr", filter = "-drop_z_below 0")
chm <- rasterize_canopy(las, 0.5, pitfree(subcircle = 0.2))
plot(las, bg = "white", size = 4)



dtm1 <- grid_terrain(las, res = 0.2, algorithm = tin())
plot(dtm1)


las <- normalize_height(las, dtm1)
chm <- grid_canopy(las, res = 0.25, p2r(0.3))
plot(chm)
chm


las <- segment_trees(las, li2012(), attribute = "treeID", uniqueness = "incremental")

ttops <- find_trees(chm, lmf(4, hmin = 2))
plot(ttops, add=TRUE)   ##plot the CHM and tree tops


convex_hulls <- delineate_crowns(las, type = c("convex"))
bbox_hulls <- delineate_crowns(las, type = c("bbox")) 
plot(convex_hulls)                


output_dir<-"path_to_save_clipped_treeswith21"
output_filename <- file.path(output_dir, paste0("AOI_ConvexHulls.tif"))
writeRaster(ttops, filename = output_filename, format = "GTiff")


#####################SAR 
tif_file_sar <- "AOI_SAR.tif" 
raster_data_sar1 <- raster(tif_file_sar)
raster_data_sar1
plot(raster_data_sar1)

#####################Optical
tif_file_optical <- "AOI_Optical.tif" 
raster_data_optical <- raster(tif_file_optical)
plot(raster_data_optical)
#####################TreeSpecies
tif_file <- "AOI_TreeSpecies.tif" 
raster_data_Species <- raster(tif_file)
plot(raster_data_Species)
las_extent <- extent(las)
las_extent
raster_data_Species <- crop(raster_data_Species, las_extent)
raster_data_Species
plot(raster_data_Species)





# Extract values from the clipped_raster for the tree coordinates
tree_coordinates <- coordinates(ttops)

tree_values1 <- extract(raster_data_Species, tree_coordinates)

tree_values1





####Replace 0 value of the tree 
# Assuming tree_values1 is your list
# Replace 0 values with their neighbors
for (i in 2:(length(tree_values1)-1)) {
  if (tree_values1[i] == 0) {
    tree_values1[i] <- ifelse(tree_values1[i-1] != 0, tree_values1[i-1], tree_values1[i+1])
  }
}

# Print the modified list
print(tree_values1)
##################################################DBH for every species

b1_spruce <- 2.2131
b2_spruce <- 0.3046
lna_spruce <- -0.5244
k_spruce <- 1.013
b_spruce <- 8.8563
m_spruce <- 19
c_spruce <- 0
d_spruce <- 0.3879



  
b1_pine <- 2.2845
b2_pine <- 0.3318
lna_pine <- -1.448
k_pine <- 1.009
b_pine <- 8.7399
m_pine <- 16
c_pine <- 0
d_pine <- 0.5624

b1_deciduous <- 1.649
b2_deciduous <- 0.373
lna_deciduous <- -2.1284
k_deciduous <- 1.004
b_deciduous <- 9.3375
m_deciduous <- 11
c_deciduous <- 0.0221
d_deciduous <- 0.2838
DBH=0
Biomass=0

calculate_DBH <- function(tree_value, Z) {
  if (tree_value == 0) {
    DBH <- 0
  } else if (tree_value == 1) {
    if (Z > asymptote_spruce) {
      Z <- 36
    }
    DBH <- ((Z - 1.3)^(1/3) * b1_spruce) / (1 - b2_spruce * (Z - 1.3)^(1/3))
  } else if (tree_value == 2) {
    if (Z > asymptote_pine) {
      Z <- 28
    }
    DBH <- ((Z - 1.3)^(1/3) * b1_pine )/ (1 - b2_pine * (Z - 1.3)^(1/3))
  } else if (tree_value == 3) {
    if (Z > asymptote_decidious) {
      Z <- 20
    }
    DBH <- ((Z - 1.3)^(1/3) * b1_deciduous) / (1 - b2_deciduous * (Z - 1.3)^(1/3))
  } else {
    stop("Invalid tree_value. Please use 0, 1, 2, or 3.")
  }
  return(DBH)
}

asymptote_spruce_Z <- 1.3 + (1 / b2_spruce)^3
asymptote_spruce_Z
asymptote_pine_Z <- 1.3 + (1 / b2_pine)^3
asymptote_pine_Z
asymptote_decidious_Z <- 1.3 + (1 / b2_deciduous)^3
asymptote_decidious_Z

#applied Asymptote 
asymptote_decidious = 20.56964
asymptote_pine = 28.67605
asymptote_spruce = 36.68428


calculate_Biomass <- function(tree_value, Z, DBH) {
  if (tree_value == 0) {
    Biomass<-0
  } else if (tree_value == 1) {
    if (Z > asymptote_spruce) {
      Z <- 36
      }
    Biomass <- k_spruce * exp(lna_spruce + b_spruce * (DBH/(DBH+m_spruce)) + c_spruce * Z + d_spruce * log (Z))
  } else if (tree_value == 2) {
    if (Z > asymptote_pine) {
      Z <- 28
      }
    Biomass <- k_pine * exp(lna_pine + b_pine * (DBH/(DBH+m_pine)) + c_pine * Z + d_pine * log (Z))
  } else if (tree_value == 3) {
    if (Z > asymptote_decidious) {
      Z <- 20
    }
    Biomass <- k_deciduous * exp(lna_deciduous + b_deciduous * (DBH/(DBH+m_deciduous)) + c_deciduous * Z + d_deciduous * log (Z))
  } else {
    stop("Invalid tree_value. Please use 0, 1, 2, or 3.")
  }
  return(Biomass)
}

for (i in 1:length(ttops$treeID)) {
  DBH[i] <- calculate_DBH(tree_values1[i], ttops$Z[i])
  Biomass[i]<- calculate_Biomass(tree_values1[i], ttops$Z[i], DBH[i])
}
DBH
Biomass



# Combine the tree ID, Z, and corresponding raster values into a data frame
tree_data_with_raster <- data.frame(treeID = ttops$treeID, Z = ttops$Z, TreeSpecies = tree_values1, tree_coordinates, DBH, Biomass)
tree_data_with_raster
# Print the resulting data frame
print(tree_data_with_raster)

# Define the path where you want to save the CSV file
csv_file_path <- "newAOI_tree_data11.csv"

# Save the data frame to a CSV file
write.csv(tree_data_with_raster, file = csv_file_path, row.names = FALSE)

# Confirm that the file has been saved
cat("Data saved to", csv_file_path, "\n")




#####################################################Patches Code

output_dir<-"AOI_clipped_Patches_optical16x16"
# Iterate through each tree and clip the optical image
for (i in 1:nrow(tree_coordinates)) {
  x <- tree_coordinates[i, 1]
  y <- tree_coordinates[i, 2]
  
  # Set the extent to a small window around the tree
  tree_extent <- extent(x - 8, x + 8, y - 8, y + 8)
  
  # Clip the optical image using the tree extent
  tree_clip <- crop(raster_data_optical, tree_extent)
  
  # Create a unique filename for each clipped image
  output_filename <- file.path(output_dir, paste0("tree_", i, ".tif"))
  
  # Save the clipped image as a GeoTIFF file
  writeRaster(tree_clip, filename = output_filename, format = "GTiff")
}



##SAR

output_dir<-"AOI_clipped_Patches_SAR16x16"
# Iterate through each tree and clip the optical image
for (i in 1:nrow(tree_coordinates)) {
  x <- tree_coordinates[i, 1]
  y <- tree_coordinates[i, 2]
  
  # Set the extent to a small window around the tree
  tree_extent <- extent(x - 8, x + 8, y - 8, y + 8)
  
  # Clip the optical image using the tree extent
  tree_clip <- crop(raster_data_sar1, tree_extent)
  
  # Create a unique filename for each clipped image
  output_filename <- file.path(output_dir, paste0("tree_", i, ".tif"))
  
  # Save the clipped image as a GeoTIFF file
  writeRaster(tree_clip, filename = output_filename, format = "GTiff")
}





##SARLBand

#####################SAR 
tif_file_sarL <- "AOI_SARL.tif" 
raster_data_sarL <- raster(tif_file_sarL)
raster_data_sarL
plot(raster_data_sarL)


output_dir<-"AOI_clipped_Patches_SARL2x2"
# Iterate through each tree and clip the optical image
for (i in 1:nrow(tree_coordinates)) {
  x <- tree_coordinates[i, 1]
  y <- tree_coordinates[i, 2]
  
  # Set the extent to a small window around the tree
  tree_extent <- extent(x - 16, x + 16, y - 16, y + 16)
  
  # Clip the optical image using the tree extent
  tree_clip <- crop(raster_data_sarL, tree_extent)
  
  # Create a unique filename for each clipped image
  output_filename <- file.path(output_dir, paste0("tree_", i, ".tif"))
  
  # Save the clipped image as a GeoTIFF file
  writeRaster(tree_clip, filename = output_filename, format = "GTiff")
}

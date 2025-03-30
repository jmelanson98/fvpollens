#Code for exploratory analysis of vegetation data and pollen availability

setwd('/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/pollencode')
#rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")
source("src/prepDF.R")
load.packages()

#load survey data
specimenData = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")
sampleEffort2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022sampledata.csv", sep = ","))
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023sampledata.csv", sep = ","))
vegData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022vegetationdata.csv", sep = ","))
vegData2022[is.na(vegData2022)] <- 0
vegData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023vegetationdata.csv", sep = ","))
vegData2023[is.na(vegData2023)] <- 0
plantcodes = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/plantcodes.csv", sep = ","))


# load landscape raster/shapefiles
landscape1 <- raster("/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/landscape/FValley_lc_1res.tif")
fv_points <- read_sf("/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/landscape/fvbombus/fvbombus_points.shp")
fv_points2022 = fv_points[fv_points$site_id %in% sampleEffort2022$sample_point,]
landcover = read_csv(paste0(here(), "/landscape/landcover.csv"))

#calculate landscape metrics
landscapeMetrics = calculateLandscapeMetrics(landscape1, fv_points2022, landcover)

#save and/or load landscape metrics
save(landscapeMetrics,
     file="saved/landscapeMetrics.Rdata")

load(file="saved/landscapeMetrics.Rdata")

#join vegetation data for both years
vegData2022[,6:ncol(vegData2022)] <- lapply(vegData2022[,6:ncol(vegData2022)], as.numeric)
vegData2023[,6:ncol(vegData2023)] <- lapply(vegData2023[,6:ncol(vegData2023)], as.numeric)
vegData = bind_rows(vegData2022, vegData2023)
vegData[is.na(vegData)] <- 0

#or leave veg data separate for years
vegData2022[is.na(vegData2022)] <- 0
vegData2023[is.na(vegData2023)] <- 0

#join sample data for both years and add julian dates
sampleEffort = bind_rows(sampleEffort2022, sampleEffort2023)
sampleEffort2022$date = paste(sampleEffort2022$day, sampleEffort2022$month, sampleEffort2022$year)
sampleEffort2022$date = gsub(" ", "", sampleEffort2022$date, fixed = TRUE)
sampleEffort2022$date <- as.POSIXlt(sampleEffort2022$date, format = "%d%b%y")
sampleEffort2022$julian_date = sampleEffort2022$date$yday
justdates = sampleEffort2022[,c("sample_id", "julian_date")]

#prep environmental matrix
environmental_df = left_join(sampleEffort2022, landscapeMetrics, by = c("sample_point" = "sample_pt"))
rownames(environmental_df) = environmental_df$sample_id
environmental_df = environmental_df[,names(environmental_df) %in% c("month", "sample_point", "julian_date", "prop_blueberry", "prop_edge", "landscape_shdi")]

#Create NMDS plots for vegetation surveys
#code by timepoint (julian date) and transect (or site, or landcover)

#prep veg matrix
#combined
rownames(vegData) = vegData$sample_id
vegMatrix = vegData[,!(names(vegData) %in% c("site", "round", "sample_point", "sample_id", "veg_observer"))]
vegPresence = vegMatrix
vegPresence[vegPresence > 0] = 1
vegPresence = vegPresence[rowSums(vegPresence) != 0,]

#separate for each year
rownames(vegData2022) = vegData2022$sample_id
rownames(vegData2023) = vegData2023$sample_id

vegMatrix2022 = vegData2022[,!(names(vegData2022) %in% c("site", "round", "sample_point", "sample_id", "veg_observer"))]
vegPresence2022 = vegMatrix2022
vegPresence2022[vegPresence2022 > 0] = 1
vegPresence2022 = vegPresence2022[rowSums(vegPresence2022) != 0,]

#remove weirdo outliers
vegPresence2022  = vegPresence2022[!(rownames(vegPresence2022) %in% c("1_PM14A", "2_W05A")),]

vegMatrix2023 = vegData2023[,!(names(vegData2023) %in% c("site", "round", "sample_point", "sample_id", "veg_observer"))]
vegPresence2023 = vegMatrix2023
vegPresence2023[vegPresence2023 > 0] = 1
vegPresence2023 = vegPresence2023[rowSums(vegPresence2023) != 0,]

#run NMDS with Jaccard distance for presence-absence (both years)
nmds_result <- metaMDS(vegPresence, distance = "jaccard", k = 2, trymax = 50, autotransform = FALSE)
#check stress value (should be < 0.2 for a good fit)
print(nmds_result$stress)

#run NMDS with jaccard distance of presence-absence (years separated)
nmds_result2022bray <- metaMDS(vegPresence2022, distance = "bray", k = 2, trymax = 100)
nmds_result2023 <- metaMDS(vegPresence2023, distance = "jaccard", k = 2, trymax = 100, autotransform = FALSE)
#check stress value (should be < 0.2 for a good fit)
print(nmds_result2022$stress)
print(nmds_result2023$stress)

#Add environmental fits (landscape complexity) with vegan envfit function
# Fit environmental variables
environmental_df = environmental_df[rownames(environmental_df) %in% rownames(vegPresence2022),colnames(environmental_df) %in% c("julian_date", "prop_blueberry", "prop_edge", "landscape_shdi")]
envfit_result <- envfit(nmds_result2022, environmental_df, permutations = 999, na.rm = TRUE)
# Display significant variables
print(envfit_result)

#Visualize
# Extract site scores
sample_scores <- as.data.frame(scores(nmds_result2022, display = "sites"))
environmental_df$sample_id = rownames(environmental_df)
sample_scores$sample_id = rownames(sample_scores)
sample_scores = left_join(sample_scores, environmental_df, by = "sample_id")

# Extract species scores
species_scores <- as.data.frame(scores(nmds_result2022, display = "species"))
species_scores$Species <- rownames(species_scores)
ggplot(data = sample_scores, aes(x = NMDS1, y = NMDS2, color = julian_date)) +
  geom_point(size = 3) +
  #geom_text(data = species_scores, aes(label = Species), color = "black", size = 3) +
  theme_minimal() +
  ggtitle("NMDS Plot")
# Extract significant vectors
env_vectors <- as.data.frame(envfit_result$vectors$arrows * sqrt(envfit_result$vectors$r))
env_vectors$Variable <- rownames(env_vectors)

# Add vectors to the plot
ggplot(data = sample_scores, aes(x = NMDS1, y = NMDS2, color = julian_date)) +
  geom_point(size = 3) +
  #geom_text(data = species_scores, aes(label = Species), color = "black", size = 3) +
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue") +
  scale_x_continuous() +
  scale_y_continuous() +
  geom_text(data = env_vectors, aes(x = NMDS1, y = NMDS2, label = Variable),
            color = "blue", vjust = -0.5) +
  theme_minimal() +
  ggtitle("NMDS Plot with Environmental Variables")


#make NMDS plot of vegetation surveys at site level

#get plant codes
colnames(plantcodes) = c("scientific_name", "common_name", "plant_code", "null")
plantcodes = plantcodes[,colnames(plantcodes) != "null"]
justcodes = plantcodes$plant_code
justcodes = gsub(" ", ".", justcodes)
extracodes = c("CRDU", "DISI", "POAV", "LInaria.spp.", "SYVU", "DRNE", "ALPO")
allcodes = c(justcodes, extracodes)
codes2022 = allcodes[allcodes %in% colnames(vegData2022)]


#exponentiate abundance values
#make sure NA's aren't 0's before running!!
vegData2022 %>%
  mutate(vegData2022, across(codes2022, as.numeric)) %>%
  mutate(across(codes2022, function(x) 10^(x-1))) %>%
  mutate(across(codes2022, ~replace(., is.na(.), 0))) -> vegData2022exp

#sum a sample of 25 veg surveys per site x round combo
vegData2022exp %>%
  group_by(site, round) %>%
  slice_sample(n=25) %>%
  summarise(across(all_of(codes2022), sum, na.rm = TRUE), .groups = "drop") -> site_level_veg

#log everything again
site_level_veg %>%
  mutate(across(codes2022, function(x) log(x+1))) -> log_site_level_veg

#prep for NMDS
log_site_level_veg = as.data.frame(log_site_level_veg)
rownames(log_site_level_veg) = paste(log_site_level_veg$site, log_site_level_veg$round, sep = "")
log_site_level_veg_red = log_site_level_veg[,colnames(log_site_level_veg) %in% codes2022]

#run NMDS
nmds_site2022 <- metaMDS(log_site_level_veg_red, distance = "bray", k = 2, trymax = 50)
#check stress value (should be < 0.2 for a good fit)
print(nmds_site2022$stress)

#add environmental fits (round and site) with vegan envfit function
# fit environmental variables
round_fit = log_site_level_veg[,colnames(log_site_level_veg) %in% c("round", "site")]
envfit_result <- envfit(nmds_site2022, round_fit, permutations = 999, na.rm = TRUE)
# Display significant variables
print(envfit_result)

# visualize
# extract site scores
sample_scores <- as.data.frame(scores(nmds_site2022, display = "sites"))
round_fit$siteround = rownames(round_fit)
sample_scores$siteround = rownames(sample_scores)
sample_scores = left_join(sample_scores, round_fit, by = "siteround")

# Extract species scores
species_scores <- as.data.frame(scores(nmds_site2022, display = "species"))
species_scores$Species <- rownames(species_scores)

ggplot(data = sample_scores, aes(x = NMDS1, y = NMDS2, color = round)) +
  geom_point(size = 3) +
  #geom_text(data = species_scores, aes(label = Species), color = "black", size = 3) +
  theme_minimal() +
  ggtitle("Site-level floral composition across time")
# Extract significant vectors
env_vectors <- as.data.frame(envfit_result$vectors$arrows * sqrt(envfit_result$vectors$r))
env_vectors$Variable <- rownames(env_vectors)

# Add vectors to the plot
ggplot(data = sample_scores, aes(x = NMDS1, y = NMDS2, color = site)) +
  geom_point(size = 3) +
  #geom_text(data = species_scores, aes(label = Species), color = "black", size = 3) +
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "blue") +
  stat_ellipse(level = 0.95) +
  scale_x_continuous() +
  scale_y_continuous() +
  geom_text(data = env_vectors, aes(x = NMDS1, y = NMDS2, label = Variable),
            color = "blue", vjust = -0.5) +
  theme_minimal() +
  ggtitle("Site level floral composition across time")


specimenData %>%
  group_by(final_id) %>%
  summarize(numbees = n(),
            unique_flowers = n_distinct(active_flower),
            flowerperbee = unique_flowers/numbees) -> speciessum

specimenData %>%
  filter(final_id != "") %>%
  filter(active_flower %in% codes2022) %>%
  group_by(final_id) %>%
  summarise(shannon_diversity = -sum((prop.table(table(active_flower))) * log(prop.table(table(active_flower)))),
            numbees = n(),
            unique_flowers = n_distinct(active_flower)) -> speciesdiv
plot(x = speciesdiv$numbees, y = speciesdiv$unique_flowers, pch = 16, col = "red")
text(speciesdiv$numbees, speciesdiv$unique_flowers, labels = speciesdiv$final_id, pos = 4)

#Model floral richness/diversity as a function of landscape interacted with time (constancy over time)
#Plot richness over time, color lines by landscape diversity
#could do this both at transect level and at subsite level (potentially more informative)

#NMDS of visitation data at transect or subsite level

#Determine how many pollen samples we have available (e.g., individuals per transect or subsite)
#do this next

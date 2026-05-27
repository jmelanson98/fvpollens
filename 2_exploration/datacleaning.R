##### Data cleaning script
### April 15, 2026
### J Melanson


### Load in packages
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)
library(sf)
library(terra)

###################################################
### Load in data
###################################################
bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"
#bombus_path = "/Users/jenna1/Documents/bombus_project/"
mix2022 = read.csv(paste0(bombus_path, "raw_data/siblingships/mixtus_sibships_2022.csv"))
mix2023 = read.csv(paste0(bombus_path, "raw_data/siblingships/mixtus_sibships_2023.csv"))
imp2022 = read.csv(paste0(bombus_path, "raw_data/siblingships/impatiens_sibships_2022.csv"))
imp2023 = read.csv(paste0(bombus_path, "raw_data/siblingships/impatiens_sibships_2023.csv"))
specimenData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022specimendata.csv"), sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023specimendata.csv"), sep = ",", header = T))
vegData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022vegetationdata.csv"), sep = ",", header = T))
vegData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023vegetationdata.csv"), sep = ",", header = T))
sampleData2022 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2022sampledata.csv"), sep = ",", header = T))
sampleData2023 = as.data.frame(read.csv(paste0(bombus_path, "raw_data/2023sampledata.csv"), sep = ",", header = T))
samplepoints = as.data.frame(read.csv(paste0(bombus_path, "raw_data/allsamplepoints.csv"), sep = ",", header = F))
fv_points = st_read(paste0(bombus_path, "landscape/fvbombus/fvbombus_points.shp"))
fv_points = st_transform(fv_points, 32610)


###################################################
### Combine close together sampling location
###################################################

# Combine sampling locations which are < 25m apart (e.g., transect was renamed in datasheet)
# targeting PM, where some transects were re-named after round 1
traps_m = data.frame(sample_point = fv_points$site_id,
                     trap_x = st_coordinates(fv_points)[,1],
                     trap_y = st_coordinates(fv_points)[,2])
dist_matrix = as.matrix(dist(traps_m[, c("trap_x", "trap_y")]))
pairs = which(dist_matrix < 20 & lower.tri(dist_matrix), arr.ind = TRUE)

close_pairs = data.frame(
  name1 = fv_points$site_id[pairs[, 1]],
  name2 = fv_points$site_id[pairs[, 2]],
  distance_m = dist_matrix[pairs])

# Rename sample points in ALL data sets
for (row in 1:nrow(close_pairs)){
  mix2022$sample_pt[mix2022$sample_pt == close_pairs$name2[row]] = close_pairs$name1[row]
  mix2023$sample_pt[mix2023$sample_pt == close_pairs$name2[row]] = close_pairs$name1[row]
  imp2022$sample_pt[imp2022$sample_pt == close_pairs$name2[row]] = close_pairs$name1[row]
  imp2023$sample_pt[imp2023$sample_pt == close_pairs$name2[row]] = close_pairs$name1[row]
  specimenData2022$sample_pt[specimenData2022$sample_pt == close_pairs$name2[row]] = close_pairs$name1[row]
  specimenData2023$sample_pt[specimenData2023$sample_pt == close_pairs$name2[row]] = close_pairs$name1[row]
  vegData2022$sample_point[vegData2022$sample_point == close_pairs$name2[row]] = close_pairs$name1[row]
  vegData2023$sample_point[vegData2023$sample_point == close_pairs$name2[row]] = close_pairs$name1[row]
  sampleData2022$sample_point[sampleData2022$sample_point == close_pairs$name2[row]] = close_pairs$name1[row]
  sampleData2023$sample_point[sampleData2023$sample_point == close_pairs$name2[row]] = close_pairs$name1[row]
}


###################################################
### Clean a couple of incorrect values
###################################################
# clean two floral values which were mis-entered from datasheets
vegData2023$TRRE[vegData2023$sample_id == "18_PM51"] = 2
vegData2023$VICR[vegData2023$sample_id == "18_SD28"] = 2

# in 2022, wild sweet peas (Lathyrus spp.) were incorrectly recorded as field peas (Pisum sativum)
# the only sites where field peas were truly observed in that year were W24, W25, W28 (beside
# Reynelda Farms large field where they were growing peas)
# this was noted in the plant list sheet in 2022 --> see OneDrive veg data sheets, plant list tab

vegData2022$`Lathyrus spp.` = vegData2022$PISA
vegData2022$`Lathyrus spp.`[vegData2022$sample_point %in% c("W24", "W25", "W28")] = NA
vegData2022$PISA[!vegData2022$sample_point %in% c("W24", "W25", "W28")] = NA


# fix one mixtus final_id
mix2022$final_id[mix2022$barcode_id == "W14_19_01"] = "B. mixtus"



################################################################
### Remove observations that occurred without veg surveys
################################################################

# first, get sample ids where we surveyed bees but not plants
missed2022 = sampleData2022$sample_id[!sampleData2022$sample_id %in% vegData2022$sample_id] # 4 surveys
missed2023 = sampleData2023$sample_id[!sampleData2023$sample_id %in% vegData2023$sample_id] # 20 surveys ugh that's so bad

# remove those surveys from cleaned specimen and sample data sheets
specimenData2022 = specimenData2022[!specimenData2022$sample_id %in% missed2022,]
specimenData2023 = specimenData2023[!specimenData2023$sample_id %in% missed2023,]
sampleData2022 = sampleData2022[!sampleData2022$sample_id %in% missed2022,]
sampleData2023 = sampleData2023[!sampleData2023$sample_id %in% missed2023,]

# also remove them from siblingships data sheets
mix2022 = mix2022[!mix2022$sample_id %in% missed2022,]
mix2023 = mix2023[!mix2023$sample_id %in% missed2023,]

imp2022 = imp2022[!imp2022$sample_id %in% missed2022,]
imp2023 = imp2023[!imp2023$sample_id %in% missed2023,]

###################################################
### Write new files
###################################################
write.csv(mix2022, "3_data/cleandata/siblingships/mixtus_sibships_2022.csv")
write.csv(mix2023, "3_data/cleandata/siblingships/mixtus_sibships_2023.csv")
write.csv(imp2022, "3_data/cleandata/siblingships/impatiens_sibships_2022.csv")
write.csv(imp2023, "3_data/cleandata/siblingships/impatiens_sibships_2023.csv")

write.csv(specimenData2022, "3_data/cleandata/fielddata/2022specimendata.csv")
write.csv(specimenData2023, "3_data/cleandata/fielddata/2023specimendata.csv")
write.csv(sampleData2022, "3_data/cleandata/fielddata/2022sampledata.csv")
write.csv(sampleData2023, "3_data/cleandata/fielddata/2023sampledata.csv")
write.csv(vegData2022, "3_data/cleandata/fielddata/2022vegetationdata.csv")
write.csv(vegData2023, "3_data/cleandata/fielddata/2023vegetationdata.csv")



###################################################
### Clean plant list
###################################################
plantlist = read.csv(paste0(bombus_path, "raw_data/plantcodes.csv"), header = FALSE)
colnames(plantlist) = c("scientific_name", "common_name", "plant_code", "notes")

# Clean plant list
plantlist$scientific_name[plantlist$scientific_name == "Aescules x carnea"] = "Aescules carnea"
plantlist$scientific_name[plantlist$scientific_name == "Convolvulus vulgaris?"] = "Convolvulus vulgaris"
plantlist$scientific_name[plantlist$scientific_name == "Crocosmia × crocosmiiflora"] = "Crocosmia crocosmiiflora"
plantlist$scientific_name[plantlist$scientific_name == "Geranium x oxonianum"] = "Geranium oxonianum"
plantlist$scientific_name[plantlist$scientific_name == "Glandularia x hybrida"] = "Glandularia hybrida"
plantlist$scientific_name[plantlist$scientific_name == "Hemerocallis x hybridis"] = "Hemerocallis hybridis"
plantlist$scientific_name[plantlist$scientific_name == "Iris x germanica"] = "Iris germanica"
plantlist$scientific_name[plantlist$scientific_name == "Lavandula x intermedia"] = "Lavandula intermedia"
plantlist$scientific_name[plantlist$scientific_name == "Leucanthemum x superbum"] = "Leucanthemum superbum"
plantlist$scientific_name[plantlist$scientific_name == "Magnolia × soulangeana"] = "Magnolia soulangeana"
plantlist$scientific_name[plantlist$scientific_name == "Mentha × piperita"] = "Mentha piperita"
plantlist$scientific_name[plantlist$scientific_name == "Petunia x atkinsiana"] = "Petunia atkinsiana"
plantlist$scientific_name[plantlist$scientific_name == "Pelargonium x hortorum"] = "Pelargonium hortorum"
plantlist$scientific_name[plantlist$scientific_name == "Pelargonium x hybridum"] = "Pelargonium hybridum"
plantlist$scientific_name[plantlist$scientific_name == "Photinia x fraseri"] = "Photinia fraseri"
plantlist$scientific_name[plantlist$scientific_name == "Rosa × damascena"] = "Rosa damascena"
plantlist$scientific_name[plantlist$scientific_name == "Aescules carnea"] = "Aesculus carnea"
plantlist$scientific_name[plantlist$scientific_name == "Calibrchoa hybrida"] = "Calibrachoa hybrida"
plantlist$scientific_name[plantlist$scientific_name == "Salomus valerandi"] = "Samolus valerandi"
plantlist$scientific_name[plantlist$scientific_name == "Brassica napus va pabularia"] = "Brassica napus"
plantlist$scientific_name[plantlist$scientific_name == "Convolvulus vulgaris"] = "Convolvulus arvensis" # synonym
plantlist$scientific_name[plantlist$scientific_name == "China Rose"] = "Rosa chinensis"
plantlist$scientific_name[plantlist$scientific_name == "Rosa chinensis"] = "China Rose"
plantlist$scientific_name[plantlist$scientific_name == "Curcurbita spp."] = "Cucurbita spp."
plantlist$scientific_name[plantlist$scientific_name == "Tripleospermum spp."] = "Tripleurospermum spp."
plantlist$scientific_name[plantlist$scientific_name == "Allioideae spp."] = "Alliaceae" # GBIF chokes on Allioideae, but finds the old name (Alliaceae)
plantlist$scientific_name[plantlist$scientific_name == "Maleae spp."] = "Rosaceae" # GBIF chokes on Maleae tribe, so go up to family


plantlist = plantlist[!plantlist$plant_code %in% c("dead", "flying", "Monocotyledon", "ALMA", "CRDU"),]
plantlist = plantlist[!plantlist$scientific_name == "",]

write.csv(plantlist, "3_data/cleandata/plantdata/observedplants.csv", row.names = FALSE)
plants = read.csv("3_data/cleandata/plantdata/observedplants.csv")


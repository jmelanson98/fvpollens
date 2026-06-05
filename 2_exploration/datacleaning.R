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


##################################################################
### Make sure visitation records and veg surveys are consistent
###################################################################

# Get specimen metadata data
specs2022 = read.csv("3_data/cleandata/fielddata/2022specimendata.csv")
specs2023 = read.csv("3_data/cleandata/fielddata/2023specimendata.csv")
sample2022 = read.csv("3_data/cleandata/fielddata/2022sampledata.csv")
sample2023 = read.csv("3_data/cleandata/fielddata/2023sampledata.csv")


# Join data frames
specs2022filt = specs2022[c("barcode_id", "active_flower", "final_id", "notes", "site", "round", "sample_pt", "sample_id", "year", "pollen")]
specs2023filt = specs2023[c("barcode_id", "active_flower", "final_id", "notes", "site", "round", "sample_pt", "sample_id", "year", "pollen")]
specs = rbind(specs2022filt, specs2023filt)

# Filter to two species
twospp = specs %>% filter(final_id == "B. mixtus" | final_id == "B. impatiens")


# filter to just flower visits
twospp = twospp %>% filter(active_flower != "flying" & 
                             active_flower != "leaf" & 
                             active_flower != "Jenna's elbow" &
                             active_flower != "road (dead)" &
                             active_flower != "nest?" &
                             active_flower != "N/A (leaf)" &
                             active_flower != "N/A" &
                             active_flower != "W01 nest" &
                             active_flower != "ground" &
                             active_flower != "dead" &
                             active_flower != "nest")
twospp$active_flower[twospp$active_flower == "dipu"] = "DIPU"

# Get floral vegetation surveys to match
veg2022 = read.csv("3_data/cleandata/fielddata/2022vegetationdata.csv")
veg2023 = read.csv("3_data/cleandata/fielddata/2023vegetationdata.csv")

veg2022 = veg2022 %>%
  select(-X, -site, -round, -sample_point, -veg_observer) %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), ~ sum(10^(.x - 1), na.rm = TRUE))) %>%
  column_to_rownames("sample_id")
veg2023 = veg2023 %>%
  select(-X, -site, -round, -sample_point, -veg_observer) %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), ~ sum(10^(.x - 1), na.rm = TRUE))) %>%
  column_to_rownames("sample_id")

veg = bind_rows(veg2022, veg2023) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  select(where(~ sum(.x) > 0))


# loop through sample events, check for mismatch
mismatched=c()

for (i in unique(twospp$sample_id)){
  # get plants
  plant_sub = veg[rownames(veg) == i,] %>%
    select(where(~ sum(.x) > 0))
  
  # get bees
  beesub = twospp %>% 
    filter(sample_id == i)
  
  # get full plant lists
  net = t(table(beesub$active_flower, beesub$final_id))
  
  #throw error if visited plant is not in plant abundance
  if(any(!colnames(net) %in% colnames(plant_sub))){
    mismatched = c(mismatched, i)
  }
}
mismatched


# Cases where the plant is just missing -- check before and after, then take average
veg2022$RULA[veg2022$sample_id == "9_W10"] = 1
veg2022$VICR[veg2022$sample_id == "7_SD28"] = 2


veg2023$LYSA[veg2023$sample_id == "25_W05"] = 2
veg2023$RUAR[veg2023$sample_id == "22_W09"] = 1

# If it's not in either, add the minimum amount
veg2023$FAES[veg2023$sample_id == "22_W17"] = 1
veg2023$SYAL[veg2023$sample_id == "18_W33"] = 1


# In cases for ROAC/RONU or TRRE/TRHY or PEMA/PELA
# --> go with veg survey, because Hazel's id was typically better than mine
specs2022$active_flower[specs2022$barcode_id == "W2_31_02"] = "RONU"
specs2022$active_flower[specs2022$barcode_id == "SD7_26_08"] = "TRHY"
specs2022$active_flower[specs2022$barcode_id == "NR10_15A_01"] = "PELA"

# CECY and CYSE denote the same plant, they are just synonyms (centaurea cynanus vs cyanus segetum)
# GBIF recognizes CECY
colnames(veg2022)[colnames(veg2022) == "CYSE"] = "CECY"
colnames(veg2023)[colnames(veg2023) == "CYSE"] = "CECY"
specs2022$active_flower[specs2022$active_flower == "CYSE"] = "CECY"
specs2023$active_flower[specs2023$active_flower == "CYSE"] = "CECY"



# pull them up and check manually
i = "19_W18"
plant_sub = veg[rownames(veg) == i,] %>%
  select(where(~ sum(.x) > 0))

beesub = twospp %>% 
  filter(sample_id == i)
beesub
plant_sub

# write to files
write.csv(specs2022, "3_data/cleandata/fielddata/2022specimendata.csv")
write.csv(specs2023, "3_data/cleandata/fielddata/2023specimendata.csv")
write.csv(veg2022, "3_data/cleandata/fielddata/2022vegetationdata.csv")
write.csv(veg2023, "3_data/cleandata/fielddata/2023vegetationdata.csv")

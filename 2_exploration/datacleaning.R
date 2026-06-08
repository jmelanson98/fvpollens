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

# somehow a character got added in the middle of the df:
vegData2023$TAPA[vegData2023$TAPA == "18_W22"] = NA

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

################### Put corrections here #########################

# First, check if it's entered incorrectly from veg datasheet
vegData2022$RUAR[vegData2022$sample_id == "4_SD27"] = 2
vegData2022$`Philadelphus spp.`[vegData2022$sample_id == "4_W24"] = 3
vegData2022$SYAL[vegData2022$sample_id == "4_W05"] = 3
vegData2022$SYAL[vegData2022$sample_id == "7_PM70"] = 3
vegData2022$SYAL[vegData2022$sample_id == "7_PM49"] = 2
vegData2022$SYAL[vegData2022$sample_id == "6_PM70"] = 3
vegData2022$SYAL[vegData2022$sample_id == "6_PM49"] = 2
vegData2022$SYOF[vegData2022$sample_id == "6_PM70"] = NA
vegData2022$SYOF[vegData2022$sample_id == "6_PM49"] = NA
vegData2022$RHDAP[vegData2022$sample_id == "3_W05"] = 4
vegData2022$SOAR[vegData2022$sample_id == "10_SD09"] = 1
vegData2022$ROSY[vegData2022$sample_id == "10_SD09"] = 4
vegData2022$BRRA[vegData2022$sample_id == "10_SD09"] = 2
vegData2022$PHVU[vegData2022$sample_id == "10_SD09"] = 1
vegData2022$ERHU[vegData2022$sample_id == "9_HR20A"] = 3
vegData2022$PHPA[vegData2022$sample_id == "9_HR20A"] = 3

vegData2023$ARDI[vegData2023$sample_id == "18_W18"] = 4
vegData2023$DIWI[vegData2023$sample_id == "18_W18"] = 2
vegData2023$AQAT[vegData2023$sample_id == "18_W18"] = 2
vegData2023$SYAL[vegData2023$sample_id == "22_PM49"] = 4
vegData2023$SYAL[vegData2023$sample_id == "22_PM48"] = NA
vegData2023$VICR[vegData2023$sample_id == "20_PM54A"] = 2
vegData2023$PIJA[vegData2023$sample_id == "12_HR20A"] = 3
vegData2023$Salix.spp.[vegData2023$sample_id == "12_HR24"] = 4
vegData2023$PIJA[vegData2023$sample_id == "12_HR18"] = 3
vegData2023$Solanum.spp.[vegData2023$sample_id == "24_W16"] = 2
vegData2023$LASA[vegData2023$sample_id == "24_W16"] = 1
vegData2023$ZEMA[vegData2023$sample_id == "24_W10"] = 5
vegData2023$RUHI[vegData2023$sample_id == "27_HR23"] = 3
vegData2023$HYMA[vegData2023$sample_id == "27_HR23"] = 4
vegData2023$RONU[vegData2023$sample_id == "18_W14"] = 2

# Names of sites go mixed up
specimenData2023$sample_pt[specimenData2023$sample_id == "18_W34"] = "W33"
specimenData2023$sample_id[specimenData2023$sample_id == "18_W34"] = "18_W33"
vegData2022$sample_point[vegData2022$sample_id == "3_W35"] = "W34"
vegData2022$sample_id[vegData2022$sample_id == "3_W35"] = "3_W34"
vegData2022$sample_point[vegData2022$sample_id == "W35A"] = "W35"
vegData2022$sample_id[vegData2022$sample_id == "3_W35A"] = "3_W35"


# Typos in plant names
colnames(vegData2022)[colnames(vegData2022) == "Rhododendron.spp."] = "Rhododendron spp."
colnames(vegData2023)[colnames(vegData2023) == "Rhododendron.spp."] = "Rhododendron spp."
colnames(vegData2022)[colnames(vegData2022) == "Prunus.spp."] = "Prunus spp."
colnames(vegData2023)[colnames(vegData2023) == "Prunus.spp."] = "Prunus spp."
colnames(vegData2022)[colnames(vegData2022) == "Acer.spp."] = "Acer spp."
colnames(vegData2023)[colnames(vegData2023) == "Acer.spp."] = "Acer spp."
colnames(vegData2023)[colnames(vegData2023) == "Salix.spp."] = "Salix spp."
colnames(vegData2023)[colnames(vegData2023) == "Solanum.spp."] = "Solanum spp."
colnames(vegData2022)[colnames(vegData2022) == "Crepis.spp."] = "Crepis spp."
colnames(vegData2023)[colnames(vegData2023) == "Crepis.spp."] = "Crepis spp."
colnames(vegData2022)[colnames(vegData2022) == "Persicaria.spp."] = "Persicaria spp."
colnames(vegData2023)[colnames(vegData2023) == "Persicaria.spp."] = "Persicaria spp."
colnames(vegData2022)[colnames(vegData2022) == "Hosta.spp."] = "Hosta spp."
colnames(vegData2023)[colnames(vegData2023) == "Hosta.spp."] = "Hosta spp."
specimenData2023$active_flower[specimenData2023$barcode_id == "HR22_01_01"] = "HYRA"
specimenData2023$active_flower[specimenData2023$barcode_id == "W18_10_01"] = "RUAR"
vegData2023$CETH[vegData2023$sample_id == "18_SD40"] = 4
vegData2023$Ceanothus.spp.[vegData2023$sample_id == "18_SD40"] = NA
vegData2022$LAPU[vegData2022$sample_id == "3_SD40"] = NA
vegData2022$LUPO[vegData2022$sample_id == "3_SD40"] = 3


# Cases where the plant is just missing -- check before and after, then take average
vegData2022$RULA[vegData2022$sample_id == "9_W10"] = 1
vegData2022$VICR[vegData2022$sample_id == "7_SD28"] = 2
vegData2022$HYPE[vegData2022$sample_id == "4_HR27"] = 1
vegData2022$VICR[vegData2022$sample_id == "8_SD28"] = 1
vegData2022$SYAL[vegData2022$sample_id == "8_W18"] = 1
vegData2022$LASY[vegData2022$sample_id == "5_HR10"] = 3
vegData2022$VICR[vegData2022$sample_id == "4_HR19A"] = 4
vegData2022$BRRA[vegData2022$sample_id == "1_SD24"] = 2
vegData2022$EPCI[vegData2022$sample_id == "6_ED07"] = 2
vegData2022$CASI[vegData2022$sample_id == "7_PM46"] = 1


vegData2023$LYSA[vegData2023$sample_id == "25_W05"] = 2
vegData2023$RUAR[vegData2023$sample_id == "22_W09"] = 1
vegData2023$RUAR[vegData2023$sample_id == "24_NR08"] = 2
vegData2023$RUPA[vegData2023$sample_id == "19_NR43"] = 2
vegData2023$ROPA[vegData2023$sample_id == "19_NR43"] = NA
vegData2023$TRRE[vegData2023$sample_id == "26_HR12"] = 3
vegData2023$RUAR[vegData2023$sample_id == "26_ED17"] = 1
vegData2023$RUAR[vegData2023$sample_id == "25_ED17"] = 1
vegData2023$LUPO[vegData2023$sample_id == "19_PM70"] = 2
vegData2023$VICR[vegData2023$sample_id == "18_HR05A"] = 3
vegData2023$RULA[vegData2023$sample_id == "22_NR22"] = 1
vegData2023$SYOF[vegData2023$sample_id == "18_ED01"] = 3
vegData2023$TRRE[vegData2023$sample_id == "27_HR12"] = 2
vegData2023$RONU[vegData2023$sample_id == "18_W30"] = 2



# If it's not in either, add the minimum value
vegData2023$SYAL[vegData2023$sample_id == "18_W33"] = 1
vegData2023$AQAT[vegData2023$sample_id == "18_W15"] = 1
vegData2023$LUPO[vegData2023$sample_id == "18_SD27"] = 2
vegData2023$VACO[vegData2023$sample_id == "13_NR05A"] = 1
vegData2023$RUSP[vegData2023$sample_id == "13_ED17A"] = 1
vegData2023$TAOF[vegData2023$sample_id == "26_SD41B"] = 1
vegData2022$CASI[vegData2022$sample_id == "8_PM46"] = 1
vegData2022$LASY[vegData2022$sample_id == "10_SD08"] = 1
vegData2022$LYSA[vegData2022$sample_id == "10_ED17"] = 1
vegData2022$VACO[vegData2022$sample_id == "3_PM45"] = 1
vegData2022$RUSP[vegData2022$sample_id == "1_W29"] = 1



# In cases for ROAC/RONU or TRRE/TRHY or PEMA/PELA
# --> go with veg survey, because Hazel's id was typically better than mine
specimenData2022$active_flower[specimenData2022$barcode_id == "W2_31_02"] = "RONU"
specimenData2022$active_flower[specimenData2022$barcode_id == "SD7_26_08"] = "TRHY"
specimenData2022$active_flower[specimenData2022$barcode_id == "NR10_15A_01"] = "PELA"
specimenData2022$active_flower[specimenData2022$barcode_id == "ED7_21A_02"] = "RULA"
specimenData2022$active_flower[specimenData2022$barcode_id == "HR1_14A_02"] = "CAHI"
specimenData2022$active_flower[specimenData2022$barcode_id == "HR1_14A_03"] = "CAHI"
specimenData2022$active_flower[specimenData2022$barcode_id == "ED8_28A_01"] = "TRHY"
specimenData2022$active_flower[specimenData2022$barcode_id == "PM2_65_01"] = "RARE"
specimenData2022$active_flower[specimenData2022$barcode_id == "ED7_24A_01"] = "TRRE"
specimenData2022$active_flower[specimenData2022$barcode_id == "ED7_24A_02"] = "TRRE"
specimenData2022$active_flower[specimenData2022$barcode_id == "PM1_20_01"] = "RONU"


# ID changed on veg sheet but not tracked on sample sheet
specimenData2023$active_flower[specimenData2023$barcode_id == "W27_16_01"] = "TRHY"
specimenData2023$active_flower[specimenData2023$barcode_id == "W27_16_02"] = "TRHY"
specimenData2023$active_flower[specimenData2023$barcode_id == "W27_16_03"] = "TRHY"
specimenData2023$active_flower[specimenData2023$barcode_id == "SD13_40_01"] = "Acer spp."
specimenData2023$active_flower[specimenData2023$barcode_id == "PM20_57_01"] = "CRTE"
specimenData2023$active_flower[specimenData2023$barcode_id == "HR27_26_01"] = "LASE"
specimenData2023$active_flower[specimenData2023$barcode_id == "ED26_23A_02"] = "RARA"
specimenData2023$active_flower[specimenData2023$barcode_id == "ED21_07_01"] = "SIORE"
specimenData2023$active_flower[specimenData2023$barcode_id == "W18_30_01"] = "ROAC"
specimenData2023$active_flower[specimenData2023$barcode_id == "W18_30_02"] = "ROAC"
specimenData2023$active_flower[specimenData2023$barcode_id == "W20_10_01"] = "RULA"
specimenData2023$active_flower[specimenData2023$barcode_id == "W18_30_01"] = "RONU"
specimenData2023$active_flower[specimenData2023$barcode_id == "W18_30_02"] = "RONU"
vegData2022$Brassicaceae[vegData2022$sample_id == "1_HR03"] = NA
vegData2022$`Rorippa spp.`[vegData2022$sample_id == "1_HR03"] = 4
vegData2022$Brassicaceae[vegData2022$sample_id == "1_SD09"] = NA
vegData2022$BRRA[vegData2022$sample_id == "1_SD09"] = 3
vegData2022$Lathyrus.spp.[vegData2022$sample_id == "5_HR06A"] = NA
vegData2022$LASY[vegData2022$sample_id == "5_HR06A"] = 2


# CECY and CYSE denote the same plant, they are just synonyms (centaurea cynanus vs cyanus segetum)
# GBIF recognizes CECY
colnames(vegData2022)[colnames(vegData2022) == "CYSE"] = "CECY"
colnames(vegData2023)[colnames(vegData2023) == "CYSE"] = "CECY"
specimenData2022$active_flower[specimenData2022$active_flower == "CYSE"] = "CECY"
specimenData2023$active_flower[specimenData2023$active_flower == "CYSE"] = "CECY"


#################### Check for mismatches ###############################
# Join data frames
specimenData2022filt = specimenData2022[c("barcode_id", "active_flower", "final_id", "notes", "site", "round", "sample_pt", "sample_id", "year", "pollen")]
specimenData2023filt = specimenData2023[c("barcode_id", "active_flower", "final_id", "notes", "site", "round", "sample_pt", "sample_id", "year", "pollen")]
specs = rbind(specimenData2022filt, specimenData2023filt)

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
                             active_flower != "nest" &
                             active_flower != "nest searching" &
                             active_flower != "grass" &
                             active_flower != "" &
                             active_flower != "road (alive)")
twospp$active_flower[twospp$active_flower == "dipu"] = "DIPU"

vegData2022_red = vegData2022 %>%
  select(-site, -round, -sample_point, -veg_observer) %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), ~ sum(10^(.x - 1), na.rm = TRUE))) %>%
  column_to_rownames("sample_id")
vegData2023_red = vegData2023 %>%
  select(-site, -round, -sample_point, -veg_observer) %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), ~ sum(10^(.x - 1), na.rm = TRUE))) %>%
  column_to_rownames("sample_id")

veg = bind_rows(vegData2022_red, vegData2023_red) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  select(where(~ sum(.x) > 0))



################### Loop to check for remaining mismatches ##################
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

# pull them up and check manually
i = "3_PM45"
plant_sub = veg[rownames(veg) == i,] %>%
  select(where(~ sum(.x) > 0))

beesub = twospp %>% 
  filter(sample_id == i)
beesub
plant_sub

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

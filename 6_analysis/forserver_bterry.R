##### Analysis script -- to get B Terry model
### June 2, 2026
### J Melanson


# Load in packages
library(dplyr)
library(brms)
library(rstan)
library(stringr)
library(ggplot2)
library(tibble)
library(vegan)
library(marginaleffects)
library(ordbetareg)
library(tidyr)


# Set color scheme for plots
faded_pale = "#D2E4D4"
faded_light = "#B6D3B8"
faded_medium = "#609F65"
faded_strong = "#4E8353"
faded_green = "#355938"
faded_dark = "#1B2D1C"


light_gold = "#F7EAC0"
lm_gold = "#F2DC97"
medium_gold = "#ECCF6F"
gold = "#E2B41D"
dark_gold = "#B99318"
darker_gold = "#907313"


# Get specimen metadata data
specs2022 = read.csv("3_data/cleandata/fielddata/2022specimendata.csv")
specs2023 = read.csv("3_data/cleandata/fielddata/2023specimendata.csv")
sample2022 = read.csv("3_data/cleandata/fielddata/2022sampledata.csv")
sample2023 = read.csv("3_data/cleandata/fielddata/2023sampledata.csv")


# Join data frames
specs2022filt = specs2022[c("barcode_id", "active_flower", "final_id", "notes", "site", "round", "sample_pt", "sample_id", "year", "pollen")]
specs2023filt = specs2023[c("barcode_id", "active_flower", "final_id", "notes", "site", "round", "sample_pt", "sample_id", "year", "pollen")]
specs = rbind(specs2022filt, specs2023filt)

# Get pollen data and filter to two species
specs$pollen_indicator = NA
specs$pollen_indicator[specs$pollen == "n"] = 0
specs$pollen_indicator[specs$pollen == "0"] = 0
specs$pollen_indicator[specs$pollen == "y"] = 1
specs$pollen_indicator[specs$pollen == "1"] = 1
specs$pollen_indicator[specs$pollen == "2"] = 1
specs$pollen_indicator[specs$pollen == "yy"] = 1
twospp = specs %>% filter(final_id == "B. mixtus" | final_id == "B. impatiens")

# Get day of year
#add julian date to sample effort data frame
sample2022$date = paste(sample2022$day, sample2022$month, sample2022$year, sep = "")
sample2022$date = as.POSIXlt(sample2022$date, format = "%d%b%y")
sample2022$doy = sample2022$date$yday
dates2022 = sample2022[c("round", "doy", "site", "sample_point")]

sample2023$date = paste(sample2023$day, sample2023$month, sample2023$year, sep = "")
sample2023$date = as.POSIXlt(sample2023$date, format = "%d%b%y")
sample2023$doy = sample2023$date$yday
dates2023 = sample2023[c("round", "doy", "site", "sample_point")]

dates = rbind(dates2022, dates2023)
twospp = left_join(twospp, dates, by = c("round", "site", "sample_pt" = "sample_point")) %>%
  filter(str_detect(notes, "queen", negate = TRUE)) %>%
  filter(str_detect(notes, "male", negate = TRUE)) %>%
  filter(!is.na(pollen_indicator))


###################################################
## Get landscape metrics
###################################################

landmets = read.csv("3_data/cleandata/landscapemetrics/calculated_metrics_exponential.csv")
iji = landmets %>%
  filter(metric == "iji") %>%
  filter(species == "B. mixtus" | species == "B. impatiens") %>%
  group_by(species, sample_pt) %>%
  summarize(iji = mean(value)/10)
ber = landmets %>%
  filter(metric == "idwBER") %>%
  group_by(species, sample_pt) %>%
  summarize(berry = sum(value)/10000)
semi = landmets %>%
  filter(metric == "idwSN" | metric == "idwEDG") %>%
  group_by(species, sample_pt) %>%
  summarize(semi = sum(value)/10000)

specs_withlandscape = twospp %>%
  left_join(iji, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(ber, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(semi, by = c("sample_pt", "final_id" = "species"))
specs_withlandscape$doy_centred = (specs_withlandscape$doy - median(specs_withlandscape$doy))/10
specs_withlandscape$year = as.factor(specs_withlandscape$year)

###################################################
## Get metabarcoding data
###################################################

# read in metabarcoding data
asvtable = read.csv("3_data/cleandata/metabarcoding/cleaned_asvtable.csv", row.names = 1)
taxa = read.csv("3_data/cleandata/metabarcoding/qiime_classifications.csv", row.names = 1)

# get long format table of ASVs with classification
asvtable$barcode_id = rownames(asvtable)
asv_long = pivot_longer(asvtable, cols = -barcode_id, 
                        names_to = "Feature.ID", 
                        values_to = "ReadCount") %>%
  filter(ReadCount > 0) %>%
  left_join(taxa, by = "Feature.ID") %>%
  separate(Taxon, 
           into = c("k", "kingdom", "phylum", "subphylum", "order", "family", "genus", "species"),
           sep = "__",
           fill = "right") %>%
  mutate(
    family = ifelse(genus == "Rhododendron_catawbiense", "Ericaceae", family),
    genus  = ifelse(genus == "Rhododendron_catawbiense", "Rhododendron", genus)) %>%
  mutate(across(c(family, genus), ~ gsub(";.*", "", .x))) %>%
  mutate(across(where(is.character), ~na_if(trimws(.x), "NA"))) %>%
  group_by(barcode_id, family, genus) %>%
  summarize(ReadCountTaxa = sum(ReadCount))

#################################################################################
##### What shapes pollen collection TURNOVER between individuals?
#################################################################################


comparison_composition = left_join(asv_long, specs_withlandscape) %>%
  filter(round %in% c(6, 7, 8, 9))

comparison_matrix = pivot_wider(comparison_composition, id_cols = barcode_id, names_from = genus, values_from = ReadCountTaxa) %>%
  column_to_rownames("barcode_id") %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

braydist = as.matrix(vegdist(comparison_matrix, method = "bray"))
jaccarddist = as.matrix(vegdist(comparison_matrix, method ="jaccard"))


braydf = as.matrix(braydist) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "bee1") %>% 
  pivot_longer(cols = -bee1, names_to = "bee2", values_to = "braydist") %>% 
  filter(bee1 > bee2)
jaccarddf = as.matrix(jaccarddist) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "bee1") %>% 
  pivot_longer(cols = -bee1, names_to = "bee2", values_to = "jaccarddist") %>% 
  filter(bee1 > bee2)

# get some metadata we want to join
metadata1 = specs_withlandscape[c("barcode_id", "site", "sample_pt", "year", "doy", "final_id")]
colnames(metadata1) = c("bee1", "site1", "sample_pt1", "year1", "doy1", "species1")
metadata2 = specs_withlandscape[c("barcode_id", "site", "sample_pt", "year", "doy", "final_id")]
colnames(metadata2) = c("bee2", "site2", "sample_pt2", "year2", "doy2", "species2")

braydf_meta = braydf %>%
  left_join(metadata1) %>%
  left_join(metadata2) %>%
  mutate(doydiff = abs(doy2-doy1)) %>%
  mutate(sametransect = as.factor(ifelse(sample_pt1 == sample_pt2, "same", "different"))) %>%
  mutate(samesite = as.factor(ifelse(site1 == site2, "same", "different"))) %>%
  mutate(species = as.factor(paste(pmin(species1, species2), pmax(species1, species2), sep = "/")))


fit = ordbetareg(
  braydist ~  0 + sametransect + samesite + species + (1|bee1) + (1|bee2),
  data = braydf_meta,
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))
saveRDS(mixtusfit, "5_models/brmsfits/bterry_comparison.rds")
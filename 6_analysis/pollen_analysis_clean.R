##### Analysis script 2
### June 9, 2026
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
library(tidyr)
library(taxize)
library(readxl)
library(rgbif)


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
  filter(str_detect(notes, "male", negate = TRUE))


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



#################################################################################
##### Get pollen metabarcoding data
#################################################################################

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

# add sample metadata
pollen_meta = left_join(specs_withlandscape, asv_long, by = "barcode_id", relationship = "one-to-many") %>%
  filter(!is.na(family)) %>%
  filter(final_id == "B. mixtus") %>%
  group_by(barcode_id, active_flower, site, round, sample_pt, sample_id, year, doy, iji, berry, semi) %>%
  summarize(major_prop = max(ReadCountTaxa)/sum(ReadCountTaxa),
            major_flower = genus[which.max(ReadCountTaxa)],
            pollen_richness = n())


#################################################################################
##### Get full plant classifications (for veg data)
#################################################################################

observedplants = read.csv("3_data/cleandata/plantdata/observedplants.csv")
observedplants = observedplants[c("scientific_name", "common_name", "plant_code")]

# GNA verifier function (safe version)
gna_verifier_safe = function(names, batch_size = 50) {
  do.call(rbind, lapply(seq(1, length(names), by = batch_size), function(i) {
    gna_verifier(names[i:min(i + batch_size - 1, length(names))])
  }))
}


# Get canonical names
resolved = gna_verifier_safe(observedplants$scientific_name)

resolved_df = tibble(
  input_plants     = observedplants$scientific_name,
  canonicalname    = resolved$matchedCanonicalFull)
resolved_df = resolved_df[!is.na(resolved_df$canonicalname), ]

# Get full classification from GBIF
get_classification = function(input_name, canonical_name) {
  tryCatch({
    
    match = name_backbone(name = canonical_name)
    
    # handle NONE match by falling back to name_suggest
    if (match$matchType == "NONE") {
      print(paste0(canonical_name, " not found. Trying name lookup."))
      suggestion = name_lookup(q = canonical_name, 
                               higherTaxonKey = 6, limit = 1)$data
      if (is.null(suggestion) || nrow(suggestion) == 0) stop("No match")
      key = suggestion$key[1]
      rank_out = "GENUS"
      canonical_out = canonical_name  # just use what we passed in
      rank_out     = suggestion$rank[1]
      status_out   = NA_character_
      authority_out = NA_character_
    } else {
      key           = match$usageKey
      canonical_out = match$canonicalName
      rank_out      = match$rank
      status_out    = match$status
      authority_out = match$authorship
    }
    
    if (is.null(key) || is.na(key)) stop("No usable key")
    
    cls = classification(key, db = "gbif")[[1]]
    if (is.null(cls) || nrow(cls) == 0) stop("No classification returned")
    
    ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    tax = setNames(rep(NA_character_, length(ranks)), ranks)
    for (r in ranks) {
      val = cls$name[cls$rank == r]
      if (length(val) > 0) tax[r] = val
    }
    
    tibble(
      input_plants  = input_name,
      canonicalname = canonical_out,
      taxa_rank     = rank_out,
      taxa_status   = status_out,
      synonym_name  = ifelse(!is.na(status_out) & status_out == "SYNONYM", 
                             canonical_name, NA_character_),
      authority     = authority_out,
      kingdom       = tax["kingdom"],
      phylum        = tax["phylum"],
      class         = tax["class"],
      order         = tax["order"],
      family        = tax["family"],
      genus         = tax["genus"],
      species       = tax["species"]
    )
    
  }, error = function(e) {
    message("Failed for: ", canonical_name, " — ", e$message)
    tibble(
      input_plants  = input_name,
      canonicalname = NA_character_,
      taxa_rank     = NA_character_,
      taxa_status   = NA_character_,
      synonym_name  = NA_character_,
      authority     = NA_character_,
      kingdom       = NA_character_,
      phylum        = NA_character_,
      class         = NA_character_,
      order         = NA_character_,
      family        = NA_character_,
      genus         = NA_character_,
      species       = NA_character_
    )
  })
}

# map to input plants
taxonomy_df = bind_rows(
  Map(get_classification, resolved_df$input_plants, resolved_df$canonicalname))

# add tayberry because it's a weirdo
tayberry_row = tibble(
  input_plants  = "Rubus fruticosus x Rubus idaeus",
  canonicalname = "Rubus fruticosus x Rubus idaeus",  # hybrid cross, not really recognized
  taxa_rank     = "UNRANKED",
  taxa_status   = "UNRESOLVED",
  synonym_name  = NA,
  kingdom       = "Plantae",
  phylum        = "Tracheophyta",
  class         = "Magnoliopsida",
  order         = "Rosales",
  family        = "Rosaceae",
  genus         = "Rubus",
  species       = NA,
  authority = NA)   # no species epithet --- hybrid
taxonomy_df = rbind(taxonomy_df, tayberry_row)
taxonomy_df = distinct(taxonomy_df)

# join back to original data
observedplants = observedplants %>%
  left_join(taxonomy_df, by = c("scientific_name" = "input_plants"))

# Do a little bit of cleaning up
observedplants$genus[observedplants$scientific_name == "Hosta spp."] = "Hosta"
observedplants$genus[observedplants$scientific_name == "Silene armeria"] = "Silene"
observedplants$genus[observedplants$scientific_name == "Rubus ursinus x Rubus idaeus"] = "Rubus"
observedplants$genus[observedplants$scientific_name == "Montia spp."] = "Montia"

# save and read out classified plants
write.csv(observedplants, "3_data/cleandata/plantdata/classifiedplants.csv")
observedplants = read.csv("3_data/cleandata/plantdata/classifiedplants.csv")



#################################################################################
##### Get MAJOR GENUS abundance for each bee
#################################################################################

# join vegetation data sheets together
veg2022 = read.csv("3_data/cleandata/fielddata/2022vegetationdata.csv")
veg2023 = read.csv("3_data/cleandata/fielddata/2023vegetationdata.csv")

floral_long1 = veg2022 %>%
  select(-veg_observer, -X) %>%
  pivot_longer(!c(site, round, sample_point, sample_id), names_to = "flower", values_to = "abundance") %>%
  mutate(across(-c(site, round, sample_point, sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(site, round, sample_point, sample_id, flower),
                ~ 10^(.-1)))
floral_long2 = veg2023 %>%
  select(-veg_observer, -X.5) %>%
  pivot_longer(!c(site, round, sample_point, sample_id), names_to = "flower", values_to = "abundance") %>%
  mutate(across(-c(site, round, sample_point, sample_id, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(site, round, sample_point, sample_id, flower),
                ~ 10^(.-1)))
floral_long = rbind(floral_long1, floral_long2)

# make sure observed plant name matches column names
observedplants$plant_code = gsub(" spp\\.", ".spp.", observedplants$plant_code)

# initiate columns
pollen_meta$landscape_major_abundance = NA
pollen_meta$local_major_abundance = NA
pollen_meta$local_richness = NA
pollen_meta$landscape_richness = NA
pollen_meta$landscape_abundance = NA
pollen_meta$local_abundance = NA

# loop and compute for all bees!
for (bee in 1:nrow(pollen_meta)){
  # get observed plants in the major genus
  active_codes = observedplants$plant_code[observedplants$genus == pollen_meta$major_flower[bee] & 
                                             !is.na(observedplants$genus)]
  
  # filter veg observations to the major genus / site / round
  veg_sub_site = floral_long %>%
    filter(site == pollen_meta$site[bee]) %>%
    filter(round == pollen_meta$round[bee]) %>%
    filter(flower %in% active_codes)
  
  # get landscape abundance
  pollen_meta$landscape_major_abundance[bee] = log10(sum(veg_sub_site$abundance, na.rm = TRUE) + 1)
  
  # filter veg observations to the major genus / sample_id
  veg_sub_sample = floral_long %>%
    filter(sample_id == pollen_meta$sample_id[bee]) %>%
    filter(flower %in% active_codes)
  
  # get landscape abundance
  pollen_meta$local_major_abundance[bee] = log10(sum(veg_sub_sample$abundance, na.rm = TRUE) + 1)
  
  # get (average) landscape floral richness (i.e., normalize by the number of surveys to prevent sampling effects)
  veg_sub_lrich = floral_long %>%
    filter(site == pollen_meta$site[bee]) %>%
    filter(round == pollen_meta$round[bee])
  pollen_meta$landscape_richness[bee] = length(unique(veg_sub_lrich$flower[!is.na(veg_sub_lrich$abundance)]))/length(unique(veg_sub_lrich$sample_id))
  pollen_meta$landscape_abundance[bee] = log10(sum(veg_sub_lrich$abundance, na.rm = TRUE) + 1)
  
  # get local floral richness
  veg_sub_locrich = floral_long %>%
    filter(sample_id == pollen_meta$sample_id[bee]) %>%
    filter(!is.na(abundance))
  pollen_meta$local_richness[bee] = length(unique(veg_sub_locrich$flower))
  pollen_meta$local_abundance[bee] = log10(sum(veg_sub_locrich$abundance, na.rm = TRUE) + 1)
}


ggplot(pollen_meta, aes(x = local_abundance, y = pollen_richness)) +
  geom_jitter(alpha = 0.5, height = 0.1, width = 0.1) +
  theme_bw()


pollen_meta$year = as.factor(pollen_meta$year)
pollen_meta$pollen_richness = pollen_meta$pollen_richness - 1

fit = brm(
  pollen_richness ~ local_richness + landscape_richness + local_major_abundance + landscape_major_abundance + iji + berry + semi + year + (1|sample_pt) + (1|site),
  data = pollen_meta,
  family = negbinomial(),
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.995))

fit_nonspecific = brm(
  pollen_richness ~ local_richness + landscape_richness + local_abundance + landscape_abundance + iji + berry + semi + year + (1|sample_pt) + (1|site),
  data = pollen_meta,
  family = negbinomial(),
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.995))

fit <- add_criterion(fit, "loo")
fit_nonspecific <- add_criterion(fit_nonspecific, "loo")
loo_compare(fit, fit_nonspecific)

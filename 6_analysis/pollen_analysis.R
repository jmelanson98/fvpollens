##### Analysis script 1
### May 28, 2026
### J Melanson


# Load in packages
library(dplyr)
library(brms)
library(rstan)
library(stringr)
library(ggplot2)
library(tibble)
library(vegan)


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


####################################################################
##### How does pollen collection (0/1) change through time??
####################################################################

ggplot(twospp, aes(x = doy, y = pollen_indicator, colour = final_id)) +
  geom_jitter(height = 0.03, alpha = 0.5) +
  xlab("Day of year") +
  facet_grid(year~final_id) +
  ylab("Collecting pollen") +
  theme_bw()

# center doy (by median day of sampling)
twospp$doy_centred = (twospp$doy - median(unique(twospp$doy)))/10
twospp$year = as.factor(twospp$year)


mixfit = brm(
  pollen_indicator ~ final_id*year + final_id*doy_centred + final_id*I(doy_centred^2) + (1|site) + (1|sample_pt),
  data = twospp,
  family = bernoulli(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))
saveRDS(mixfit, "analysis/brmsfits/sibcount_mixtus.rds")


####################################################################
##### How is pollen collection spp richness different between spp?
####################################################################

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
           into = c("phylum", "na_col", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right") %>%
  mutate(across(where(is.character), ~na_if(trimws(.x), "NA"))) %>%
  group_by(barcode_id, family, genus, species) %>%
  summarize(ReadCountTaxa = sum(ReadCount))


# add sample metadata
individual_pollens = left_join(asv_long, specs, by = "barcode_id")

# get richness per individual
pollen_richness = individual_pollens %>%
  group_by(barcode_id, final_id, site, round, sample_pt, year) %>%
  summarize(pollrich = n())

# for species comparisons, filter to phenological window!
unique(pollen_richness$round[pollen_richness$final_id == "B. impatiens"])
pollen_richness_comparison = pollen_richness %>% filter(round %in% c(6, 7, 8, 9))
pollen_richness_comparison$final_id = as.factor(pollen_richness_comparison$final_id)
pollen_richness_comparison$site = as.factor(pollen_richness_comparison$site)

# which species has richer pollen?
fit = brm(
  pollrich ~ final_id + (1|site) + (1|sample_pt),
  data = pollen_richness_comparison,
  family = negbinomial(),
  chains = 4, cores = 4,
  warmup = 3000, iter = 6000,
  control = list(adapt_delta = 0.99))



####################################################################
##### How is pollen collection DIVERSITY different between spp?
####################################################################

pollens_wide = individual_pollens %>%
  ungroup() %>%
  filter(round %in% c(6, 7, 8, 9)) %>%
  mutate(taxa = paste(family, genus, species)) %>%
  pivot_wider(id_cols = barcode_id,
              names_from = taxa,
              values_from = ReadCountTaxa) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))
  
pollen_mat = pollens_wide %>%
  column_to_rownames("barcode_id") %>%  
  as.matrix()                          

pollens_wide$shannon_div = diversity(pollen_mat, index = "shannon")

pollens_wide = left_join(pollens_wide, specs, by = "barcode_id")
pollens_wide$site = as.factor(pollens_wide$site)
pollens_wide$final_id = as.factor(pollens_wide$final_id)
pollens_wide$round = as.factor(pollens_wide$round)
pollens_wide$shannon_div_zo = pollens_wide$shannon_div/max(pollens_wide$shannon_div)

# which species has more diverse pollen?
# which species has richer pollen?
fitdiv = brm(
  shannon_div_zo ~ final_id + round + (1|site) + (1|sample_pt),
  data = pollens_wide,
  family = zero_one_inflated_beta(),
  chains = 4, cores = 4)

####################################################################
##### How about plot an NMDS for composition?
####################################################################

# this time group only to genus
asv_long = pivot_longer(asvtable, cols = -barcode_id, 
                        names_to = "Feature.ID", 
                        values_to = "ReadCount") %>%
  filter(ReadCount > 0) %>%
  left_join(taxa, by = "Feature.ID") %>%
  separate(Taxon, 
           into = c("phylum", "na_col", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right") %>%
  mutate(across(where(is.character), ~na_if(trimws(.x), "NA"))) %>%
  group_by(barcode_id, family, genus) %>%
  summarize(ReadCountTaxa = sum(ReadCount))


# add sample metadata
individual_pollens = left_join(asv_long, specs, by = "barcode_id")
pollens_wide = individual_pollens %>%
  ungroup() %>%
  filter(round %in% c(6, 7, 8, 9)) %>%
  mutate(taxa = paste(family, genus)) %>%
  pivot_wider(id_cols = barcode_id,
              names_from = taxa,
              values_from = ReadCountTaxa) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

pollen_mat = pollens_wide %>%
  column_to_rownames("barcode_id") %>%  
  as.matrix()   


meta = specs %>%
  filter(barcode_id %in% rownames(pollen_mat)) %>%
  select(c(barcode_id, final_id, site, round)) %>%
  arrange(match(barcode_id, rownames(pollen_mat)))


dist_mat = vegdist(pollen_mat, method = "jaccard", binary = TRUE)

adonis2(dist_mat ~ final_id + round + site,
        data = meta,
        permutations = 999, 
        by = "margin") 

adonis2(dist_mat ~ final_id, 
        data = meta,
        permutations = 999,
        strata = meta$site)

h = how(blocks = interaction(meta$round, meta$site),
         nperm = 999)

bd <- betadisper(dist_mat, meta$final_id)
permutest(bd, permutations = h)


meta_sub <- meta %>% filter(round != 9)
dist_mat_sub <- as.dist(as.matrix(dist_mat)[meta_sub$barcode_id, meta_sub$barcode_id])

adonis2(dist_mat_sub ~ final_id,
        blocks = interaction(meta_sub$round, meta_sub$site),
        data = meta_sub,
        permutations = h)



meta_sub <- meta %>% filter(round != 9)
dist_mat_sub <- as.dist(as.matrix(dist_mat)[rownames(as.matrix(dist_mat)) %in% meta_sub$barcode_id,
                                            rownames(as.matrix(dist_mat)) %in% meta_sub$barcode_id])

# rebuild h on the subset
h <- how(blocks = interaction(meta_sub$round, meta_sub$site),
         nperm = 999)

adonis2(dist_mat_sub ~ final_id,
        data = meta_sub,
        permutations = h)

#################
# Try some visitation stuff


twospp = left_join(twospp, dates, by = c("round", "site", "sample_pt" = "sample_point")) %>%
  filter(str_detect(notes, "queen", negate = TRUE)) %>%
  filter(str_detect(notes, "male", negate = TRUE))

ggplot(twospp, aes(x = doy)) +
  geom_histogram() +
  facet_grid(final_id ~ year) + 
  theme_bw()






###################################################
## How do flower visitation networks vary by round
###################################################

library(bipartite)

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
                             active_flower != "dead")
twospp$active_flower[twospp$active_flower == "dipu"] = "DIPU"

# network per round x site
specializations = expand.grid(round = c(2:10, 17:27),
                              site = unique(twospp$site))
specializations$specialization_mix = NA
specializations$specialization_imp = NA

for (i in 1:nrow(specializations)){
  sub = twospp %>% 
    filter(site == specializations$site[i]) %>%
    filter(round == specializations$round[i])
  net = table(sub$active_flower, sub$final_id)
  specializations$specialization_mix[i] = specieslevel(net, index = "d")$`higher level`$d[2]
  specializations$specialization_imp[i] = specieslevel(net, index = "d")$`higher level`$d[1]
}

hprimes = specializations %>%
  filter(!is.na(specialization_mix) & !is.na(specialization_imp)) %>%
  pivot_longer(cols = c("specialization_mix", "specialization_imp"), 
               names_to = "species", 
               values_to = "hprime")

sitedates = dates %>%
  group_by(site, round) %>%
  summarize(doy=mean(doy))
roundyears = twospp %>% distinct(round, year)


hprimes = hprimes %>%
  left_join(sitedates) %>%
  left_join(roundyears)
hprimes$doy_centred = (hprimes$doy - median(hprimes$doy))/10
hprimes$year = as.factor(hprimes$year)
hprimes$hprime[near(hprimes$hprime, 1)] = 1

fit = brm(
  hprime ~ species*year + species*doy_centred + species*I(doy_centred^2) + (1|site),
  data = hprimes,
  family = zero_one_inflated_beta(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))




# is visitation composition different?
# community matrix: rows = bee species × round × site, cols = plant species visited
visit_mat <- twospp %>%
  group_by(final_id, round, year, site, active_flower) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = active_flower, values_from = n, values_fill = 0)
plant_cols = unique(twospp$active_flower)

visit_community <- visit_mat %>% 
  ungroup() %>%
  select(-final_id, -year, -round, -site) %>% 
  as.matrix()
visit_meta <- visit_mat %>% select(final_id, year, round, site)

visit_hell <- decostand(visit_community, method = "hellinger")
dist_visit <- vegdist(visit_hell, method = "bray")


h <- how(blocks = interaction(visit_meta$round, visit_meta$site, visit_meta$year), nperm = 999)

adonis2(dist_visit ~ final_id, data = visit_meta, permutations = h)


adonis2(dist_visit ~ final_id + round + site + year,
        data = visit_meta,
        permutations = 999, 
        by = "margin")      



###################################################
## Now try some stuff with landscape metrics
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
  left_join(semi, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(pollen_richness_comparison)
forpoll = specs_withlandscape %>% filter(!is.na(pollrich))
forpoll$doy_centred = (forpoll$doy - median(forpoll$doy))


fit = brm(
  pollrich ~ final_id*berry + final_id*semi + final_id*iji + final_id*doy_centred + (1|site) + (1|sample_pt),
  data = forpoll,
  family = negbinomial(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))



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
           into = c("phylum", "na_col", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right") %>%
  mutate(across(where(is.character), ~na_if(trimws(.x), "NA"))) %>%
  group_by(barcode_id, family, genus, species) %>%
  summarize(ReadCountTaxa = sum(ReadCount))


# add sample metadata
individual_pollens = left_join(asv_long, specs, by = "barcode_id")

# get richness per individual
pollen_richness = individual_pollens %>%
  group_by(barcode_id, final_id, site, round, sample_pt, year) %>%
  summarize(pollrich = n())
pollen_richness$final_id = as.factor(pollen_richness$final_id)
pollen_richness$site = as.factor(pollen_richness$site)


specs_withlandscape = twospp %>%
  left_join(iji, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(ber, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(semi, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(pollen_richness) %>%
  filter(final_id == "B. mixtus") %>%
  filter(!is.na(pollrich))
specs_withlandscape$doy_centred = (specs_withlandscape$doy - median(specs_withlandscape$doy))/10
specs_withlandscape$year = as.factor(specs_withlandscape$year)


fit = brm(
  pollrich ~ berry + semi + iji + year + doy_centred + I(doy_centred^2) + (1|site) + (1|sample_pt),
  data = specs_withlandscape,
  family = negbinomial(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))


specs_withlandscape = twospp %>%
  left_join(iji, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(ber, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(semi, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(pollen_richness) %>%
  filter(final_id == "B. impatiens") %>%
  filter(!is.na(pollrich))
specs_withlandscape$doy_centred = (specs_withlandscape$doy - median(specs_withlandscape$doy))/10


fit = brm(
  pollrich ~ berry + semi + iji + doy_centred + I(doy_centred^2) + (1|site) + (1|sample_pt),
  data = specs_withlandscape,
  family = negbinomial(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))



impcollection = twospp %>%
  left_join(iji, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(ber, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(semi, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(pollen_richness) %>%
  filter(final_id == "B. impatiens")
impcollection$doy_centred = (impcollection$doy - median(impcollection$doy))/10

pollencollection = brm(
  pollen_indicator ~ year + doy_centred + I(doy_centred^2) + berry + iji + semi + (1|site) + (1|sample_pt),
  data = impcollection,
  family = bernoulli(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))


mixcollection = twospp %>%
  left_join(iji, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(ber, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(semi, by = c("sample_pt", "final_id" = "species")) %>%
  left_join(pollen_richness) %>%
  filter(final_id == "B. mixtus")
mixcollection$doy_centred = (mixcollection$doy - median(mixcollection$doy))/10
mixcollection$year = as.factor((mixcollection$year))

pollencollection = brm(
  pollen_indicator ~ year + doy_centred + I(doy_centred^2) + berry + iji + semi + (1|site) + (1|sample_pt),
  data = mixcollection,
  family = bernoulli(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))

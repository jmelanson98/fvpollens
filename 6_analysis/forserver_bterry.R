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

setwd("~/projects/def-ckremen/melanson/fvpollens")


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

# add sample metadata
pollen_richness = left_join(specs_withlandscape, asv_long, by = "barcode_id", relationship = "one-to-many") %>%
  filter(!is.na(family)) %>%
  group_by(barcode_id, final_id, site, round, sample_pt, year, doy, doy_centred, iji, berry, semi) %>%
  summarize(pollrich = n()) %>%
  mutate(pollrich = pollrich -1)
pollen_richness$final_id = as.factor(pollen_richness$final_id)
pollen_richness$site = as.factor(pollen_richness$site)
pollen_richness$year = as.factor(pollen_richness$year)


comparisonrichness = pollen_richness %>% filter(round %in% c(6, 7, 8, 9))
mixtusrichness = pollen_richness %>% filter(final_id == "B. mixtus")

#################################################################################
##### What shapes pollen collection TURNOVER between individuals?
#################################################################################

# Get jaccard distance between pollen loads
mixtus_pollen = left_join(asv_long, specs_withlandscape) %>%
  filter(final_id == "B. mixtus")

mixtus_matrix = pivot_wider(mixtus_pollen, id_cols = barcode_id, names_from = genus, values_from = ReadCountTaxa) %>%
  column_to_rownames("barcode_id") %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

jaccarddist = as.matrix(vegdist(mixtus_matrix, method ="jaccard"))

jaccarddf = as.matrix(jaccarddist) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "bee1") %>% 
  pivot_longer(cols = -bee1, names_to = "bee2", values_to = "jaccarddist") %>% 
  filter(bee1 > bee2)

# get some metadata we want to join
metadata1 = specs_withlandscape[c("barcode_id", "site", "sample_pt", "sample_id", "year", "doy", "iji", "berry", "semi")]
colnames(metadata1) = c("bee1", "site1", "sample_pt1", "sample_id1", "year1", "doy1", "iji1", "berry1", "semi1")
metadata2 = specs_withlandscape[c("barcode_id", "site", "sample_pt", "sample_id", "year", "doy", "iji", "berry", "semi")]
colnames(metadata2) = c("bee2", "site2", "sample_pt2", "sample_id2", "year2", "doy2", "iji2", "berry2", "semi2")

jacdf_meta = jaccarddf %>%
  left_join(metadata1) %>%
  left_join(metadata2) %>%
  mutate(doydiff = abs(doy2-doy1)) %>%
  filter(sample_id1 == sample_id2)
jacdf_meta$doy_centered = (jacdf_meta$doy1 - 180)/10


fit = ordbetareg(
  jaccarddist ~  iji1 + berry1 + semi1 + year1*doy_centered + year1*I(doy_centered^2) + (1|bee1) + (1|bee2) + (1|site1),
  data = jacdf_meta,
  chains = 4, cores = 4,
  controls(list=c(adapt_delta = 0.99)))
saveRDS(fit, "5_models/brmsfits/ordbetareg_mixtus.rds")

pp_check(fit, type = "hist")
hypothesis(fit, "semi1 > 0")


# Get effect of blueberry
berry_seq = seq(min(jacdf_meta$berry1),
                max(jacdf_meta$berry1),
                length.out = 50)

berrypred = lapply(berry_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- jacdf_meta |> mutate(berry1 = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(berry1) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows()

#create axis values for standardized variables
labs.berry = pretty(jacdf_meta$berry1*10, n = 5)
axis.berry = labs.berry/10

# Plot
bberryplot = ggplot() +
  geom_line(data = berrypred, aes(x = berry1, y = estimate), colour = faded_green, linewidth = 1) +
  geom_ribbon(data = berrypred, aes(x = berry1, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_jitter(data = jacdf_meta, aes(x = berry1, y = jaccarddist), colour = "black", width = 0.5, alpha = 0.2) +
  xlab("Distance-weighted berry cultivation") +
  ylab("Jaccard distance") +
  scale_x_continuous(
    breaks = axis.berry,
    labels = labs.berry) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))



# Get effects of DOY
doy_seq = seq(min(jacdf_meta$doy_centered),
                max(jacdf_meta$doy_centered),
                length.out = 50)

doypred = lapply(doy_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- jacdf_meta |> mutate(doy_centered = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(year1, doy_centered) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows()

#create axis values for standardized variables
labs.doy = pretty(jacdf_meta$doy_centered*10 + 180, n = 5)
axis.doy = (labs.doy - 180)/10

# Plot
doyplot = ggplot() +
  geom_line(data = doypred[doypred$year1 == 2023,], aes(x = doy_centered, y = estimate), colour = faded_green, linewidth = 1) +
  geom_ribbon(data = doypred[doypred$year1 == 2023,], aes(x = doy_centered, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_jitter(data = jacdf_meta, aes(x = doy_centered, y = jaccarddist), colour = "black", width = 0.5, alpha = 0.2) +
  xlab("Day of year") +
  ylab("Jaccard distance") +
  scale_x_continuous(
    breaks = axis.doy,
    labels = labs.doy) +
  facet_grid(year1~1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

grid = bberryplot + doyplot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave("7_figures/manuscript_figures/mixtus_betadiv.jpg", grid,
       height = 2000, width = 3000, units = "px")

####################################################################
####################################################################
##### SPLIT BEES INTO EARLY -- MID -- LATE SEASON
####################################################################
####################################################################

q33 = quantile(mixtusrichness$doy[mixtusrichness$year == 2022], 0.33)
q67 = quantile(mixtusrichness$doy[mixtusrichness$year == 2022], 0.67)

byseason2022 = mixtusrichness %>%
  filter(year == 2022) %>%
  mutate(season = case_when(
    doy <= q33 ~ "early",
    doy <= q67 ~ "mid",
    TRUE ~ "late"))

q33 = quantile(mixtusrichness$doy[mixtusrichness$year == 2023], 0.33)
q67 = quantile(mixtusrichness$doy[mixtusrichness$year == 2023], 0.67)

byseason2023 = mixtusrichness %>%
  filter(year == 2023) %>%
  mutate(season = case_when(
    doy <= q33 ~ "early",
    doy <= q67 ~ "mid",
    TRUE ~ "late"))

byseason = rbind(byseason2023, byseason2022)

####################################################################
##### How does pollen turnover vary by landscape x season?
####################################################################

# Get jaccard distance between pollen loads
mixtus_pollen = left_join(asv_long, byseason) %>%
  filter(final_id == "B. mixtus")

mixtus_matrix = pivot_wider(mixtus_pollen, id_cols = barcode_id, names_from = genus, values_from = ReadCountTaxa) %>%
  column_to_rownames("barcode_id") %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

jaccarddist = as.matrix(vegdist(mixtus_matrix, method ="jaccard"))

jaccarddf = as.matrix(jaccarddist) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "bee1") %>% 
  pivot_longer(cols = -bee1, names_to = "bee2", values_to = "jaccarddist") %>% 
  filter(bee1 > bee2)

# get some metadata we want to join
metadata1 = byseason[c("barcode_id", "site", "sample_pt", "year", "doy", "iji", "berry", "semi", "season")]
colnames(metadata1) = c("bee1", "site1", "sample_pt1", "year1", "doy1", "iji1", "berry1", "semi1", "season1")
metadata2 = byseason[c("barcode_id", "site", "sample_pt", "year", "doy", "iji", "berry", "semi", "season")]
colnames(metadata2) = c("bee2", "site2", "sample_pt2", "year2", "doy2", "iji2", "berry2", "semi2", "season2")

jacdf_meta = jaccarddf %>%
  left_join(metadata1) %>%
  left_join(metadata2) %>%
  mutate(doydiff = abs(doy2-doy1)) %>%
  filter(sample_pt1 == sample_pt2) %>%
  filter(season1 == season2)


fit = ordbetareg(
  jaccarddist ~  season1*iji1 + season1*berry1 + season1*semi1 + year1 + doydiff + (1|bee1) + (1|bee2) + (1|site1),
  data = jacdf_meta,
  chains = 4, cores = 4,
  controls(list=c(adapt_delta = 0.99)))
saveRDS(fit, "5_models/brmsfits/ordbetareg_seasonxlandscape.rds")

summary(fit)

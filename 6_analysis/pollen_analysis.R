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
library(marginaleffects)


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


#############################################################################
##### How does pollen collection (0/1) change through time and wrt landscape?
#############################################################################

fit = brm(
  pollen_indicator ~ final_id*year + final_id*doy_centred +
                    final_id*iji + final_id*berry + final_id*semi + (1|site) + (1|sample_pt),
  data = specs_withlandscape,
  family = bernoulli(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))
saveRDS(fit, "5_models/brmsfits/pollencollection.rds")
fit = readRDS("5_models/brmsfits/pollencollection.rds")



doy_seq <- seq(min(specs_withlandscape$doy_centred),
               max(specs_withlandscape$doy_centred),
               length.out = 50)

doypred <- lapply(doy_seq, function(d) {
  
  # Swap in this doy value, keep all other covariates real
  tmp <- specs_withlandscape |> mutate(doy_centred = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp) |>
    as.data.frame() |>
    group_by(final_id, doy_centred) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows()

#create axis values for standardized variables
labs.doy = pretty(specs_withlandscape$doy, n = 6)
axis.doy = (labs.doy-median(specs_withlandscape$doy)) / 10

# Plot
pollencollection = ggplot() +
  geom_line(data = doypred[doypred$final_id == "B. mixtus",], aes(x = doy_centred, y = estimate, colour = final_id), linewidth = 1) +
  geom_ribbon(data = doypred[doypred$final_id == "B. mixtus",], aes(x = doy_centred, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.5) +
  geom_line(data = doypred[doypred$final_id == "B. impatiens",], aes(x = doy_centred, y = estimate, colour = final_id), linetype = "dashed", linewidth = 1) +
  geom_ribbon(data = doypred[doypred$final_id == "B. impatiens",], aes(x = doy_centred, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.5) +
  geom_jitter(data = specs_withlandscape, aes(x = doy_centred, y = pollen_indicator, colour = final_id), height = 0.03, alpha = 0.2) +
  scale_colour_manual(values=c(gold, faded_green)) +
  scale_fill_manual(values=c("grey", faded_green)) +
  scale_x_continuous(
    breaks = axis.doy,
    labels = labs.doy,
    limits = (c(120, 240) - median(specs_withlandscape$doy)) / 10) +
  xlab("Day of year") +
  facet_grid(~final_id) +
  ylab("Collecting pollen") +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))
ggsave("7_figures/manuscript_figures/pollen_collection.png", pollencollection,
       width = 2000, height = 1000, units = "px")



#################################################################################
##### How is pollen GENUS richness different between spp?
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



comparisonfit = brm(
  pollrich ~ final_id*berry + final_id*semi + final_id*iji + final_id*doy_centred + (1|site) + (1|sample_pt),
  data = comparisonrichness,
  family = negbinomial(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))
hypothesis(comparisonfit, "doy_centred + final_idB.mixtus:doy_centred = 0")
# no difference in any var except impatiens richness increases through time
pp_check(comparisonfit, type = "bars")
pp_check(comparisonfit)

spp_diff = avg_predictions(comparisonfit,
                variables = list(final_id = c("B. mixtus", "B. impatiens")),
                newdata = comparisonrichness,
                re_formula = NA) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

spp_diff_plot = ggplot(spp_diff) +
  geom_rect(aes(xmin = as.numeric(factor(final_id)) - 0.3, 
                xmax = as.numeric(factor(final_id)) + 0.3, 
                ymin = conf.low, ymax = conf.high,
                fill = final_id), alpha = 0.5) +
  geom_segment(aes(x = as.numeric(factor(final_id)) - 0.3, 
                   xend = as.numeric(factor(final_id)) + 0.3, 
                   y = estimate, yend = estimate, 
                   colour = final_id), linewidth = 1) +
  scale_x_continuous(breaks = 1:2, 
                     labels = levels(factor(spp_diff$final_id))) +
  scale_colour_manual(values = c(dark_gold, faded_dark)) +
  scale_fill_manual(values = c(medium_gold, faded_medium)) +
  labs(y = "Pollen richness (no. of genera)", x = "Bumblebee species") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, face = "italic"),
        axis.text.y = element_text(size = 11))


# DOY changes
doy_seq = seq(min(comparisonrichness$doy_centred),
               max(comparisonrichness$doy_centred),
               length.out = 10)

doypred <- lapply(doy_seq, function(d) {
  
  # Swap in this doy value, keep all other covariates real
  tmp <- comparisonrichness |> mutate(doy_centred = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(comparisonfit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(final_id, doy_centred) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.doy = pretty(specs_withlandscape$doy, n = 3)
axis.doy = (labs.doy-median(specs_withlandscape$doy)) / 10

# Plot
richnessdoycomparison = ggplot() +
  geom_line(data = doypred[doypred$final_id == "B. mixtus",], aes(x = doy_centred, y = estimate, colour = final_id), linewidth = 1, linetype = "dashed") +
  geom_ribbon(data = doypred[doypred$final_id == "B. mixtus",], aes(x = doy_centred, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.5) +
  geom_line(data = doypred[doypred$final_id == "B. impatiens",], aes(x = doy_centred, y = estimate, colour = final_id), linewidth = 1) +
  geom_ribbon(data = doypred[doypred$final_id == "B. impatiens",], aes(x = doy_centred, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.3) +
  geom_jitter(data = comparisonrichness, aes(x = doy_centred, y = pollrich, colour = final_id), height = 0.03, alpha = 0.5) +
  scale_colour_manual(values=c(faded_green, gold)) +
  scale_fill_manual(values=c("grey", gold)) +
  scale_x_continuous(
    breaks = axis.doy,
    labels = labs.doy,
    limits = (c(195, 230) - median(comparisonrichness$doy)) / 10) +
  xlab("Day of year") +
  facet_grid(~final_id) +
  ylab("Pollen richness (no. of genera)") +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

grid = spp_diff_plot + richnessdoycomparison +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave("7_figures/appendix_figures/pollen_richness_comparison.png", grid,
       width = 3000, height = 1500, units = "px")




#################################################################################
##### How is pollen GENUS richness different between landscapes/years?
#################################################################################

mixtusfit = brm(
  pollrich ~ berry + semi + iji + year*doy_centred + year*I(doy_centred^2) + (1|site) + (1|sample_pt),
  data = mixtusrichness,
  family = negbinomial(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))
saveRDS(mixtusfit, "5_models/brmsfits/mixtus_richness.rds")

pp_check(mixtusfit, type = "bars")
pp_check(mixtusfit)
pp_check(mixtusfit, type = "stat", stat = "max")

hypothesis(mixtusfit, "doy_centred + year2023:doy_centred = 0")
hypothesis(mixtusfit, "Idoy_centredE2 + year2023:Idoy_centredE2 = 0")
hypothesis(mixtusfit, "iji > 0")


# Get effect of DOY
doy_seq = seq(min(mixtusrichness$doy_centred),
              max(mixtusrichness$doy_centred),
              length.out = 50)

doypred = lapply(doy_seq, function(d) {
  
  # Swap in this doy value, keep all other covariates real
  tmp <- mixtusrichness |> mutate(doy_centred = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(mixtusfit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(year, doy_centred) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.doy = pretty(specs_withlandscape$doy, n = 5)
axis.doy = (labs.doy-median(specs_withlandscape$doy)) / 10

# Plot
mixtusdoyplot = ggplot() +
  geom_line(data = doypred[doypred$year == 2022,], aes(x = doy_centred, y = estimate), colour = faded_green, linewidth = 1, linetype = "dashed") +
  geom_ribbon(data = doypred[doypred$year == 2022,], aes(x = doy_centred, ymin = conf.low, ymax = conf.high), fill = "grey", alpha = 0.3) +
  geom_line(data = doypred[doypred$year == 2023,], aes(x = doy_centred, y = estimate), colour = faded_green, linewidth = 1) +
  geom_ribbon(data = doypred[doypred$year == 2023,], aes(x = doy_centred, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_point(data = mixtusrichness, aes(x = doy_centred, y = pollrich + 1), colour = "black", height = 0.03, alpha = 0.3) +
  scale_x_continuous(
    breaks = axis.doy,
    labels = labs.doy) +
  xlab("Day of year") +
  facet_grid(year~final_id) +
  ylab("") +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))


# Get effect of IJI
iji_seq = seq(min(mixtusrichness$iji),
              max(mixtusrichness$iji),
              length.out = 50)

ijipred = lapply(iji_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- mixtusrichness |> mutate(iji = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(mixtusfit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(iji) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.iji = pretty(mixtusrichness$iji*10, n = 5)
axis.iji = labs.iji/10

# Plot
mixtusijiplot = ggplot() +
  geom_line(data = ijipred, aes(x = iji, y = estimate), colour = faded_green, linewidth = 1, linetype = "dashed") +
  geom_ribbon(data = ijipred, aes(x = iji, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_point(data = mixtusrichness, aes(x = iji, y = pollrich + 1), colour = "black", height = 0.03, alpha = 0.3) +
  scale_x_continuous(
    breaks = axis.iji,
    labels = labs.iji) +
  xlab("Landscape interspersion") +
  ylab("Pollen richness (no. of genera)") +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))


# get grid
grid = mixtusijiplot + mixtusdoyplot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave("7_figures/manuscript_figures/pollen_richness_mixtus_landscape.png", grid,
       width = 3000, height = 2000, units = "px")




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




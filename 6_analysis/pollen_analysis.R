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
  geom_jitter(data = specs_withlandscape, aes(x = doy_centred, y = pollen_indicator), colour = "black", height = 0.03, alpha = 0.1) +
  scale_colour_manual(values=c(faded_green)) +
  scale_fill_manual(values=c(faded_green)) +
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


# Now just for mixtus
mixtuscollection = ggplot() +
  geom_line(data = doypred[doypred$final_id == "B. mixtus",], aes(x = doy_centred, y = estimate, colour = final_id), linewidth = 1) +
  geom_ribbon(data = doypred[doypred$final_id == "B. mixtus",], aes(x = doy_centred, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.5) +
  geom_jitter(data = specs_withlandscape[specs_withlandscape$final_id == "B. mixtus",], aes(x = doy_centred, y = pollen_indicator), colour = "black", height = 0.03, alpha = 0.1) +
  scale_colour_manual(values=c(faded_green)) +
  scale_fill_manual(values=c(faded_green)) +
  scale_x_continuous(
    breaks = axis.doy,
    labels = labs.doy,
    limits = (c(120, 240) - median(specs_withlandscape$doy)) / 10) +
  xlab("Day of year") +
  ylab("") +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

berry_seq <- seq(min(specs_withlandscape$berry),
                 max(specs_withlandscape$berry),
                 length.out = 50)

berrypred <- lapply(berry_seq, function(d) {
  
  # Swap in this doy value, keep all other covariates real
  tmp <- specs_withlandscape |> mutate(berry = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp) |>
    as.data.frame() |>
    group_by(final_id, berry) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows()

#create axis values for standardized variables
labs.berry = pretty(specs_withlandscape$berry*10, n = 6)
axis.berry = labs.berry/10

mixtuscollection_berry = ggplot() +
  geom_line(data = berrypred[berrypred$final_id == "B. mixtus",], aes(x = berry, y = estimate, colour = final_id), linewidth = 1) +
  geom_ribbon(data = berrypred[berrypred$final_id == "B. mixtus",], aes(x = berry, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.5) +
  geom_jitter(data = specs_withlandscape[specs_withlandscape$final_id == "B. mixtus",], aes(x = berry, y = pollen_indicator), colour = "black", height = 0.03, alpha = 0.1) +
  scale_colour_manual(values=c(faded_green)) +
  scale_fill_manual(values=c(faded_green)) +
  xlab("Distance-weighted blueberry area") +
  ylab("P(collecting pollen)") +
  scale_x_continuous(
    breaks = axis.berry,
    labels = labs.berry,
    limits = (c(0, 250) / 10)) +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

grid = mixtuscollection_berry + mixtuscollection +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))
ggsave("7_figures/manuscript_figures/pollen_collection_mixtus.png", grid,
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
##### Create some rarefaction curves of pollen collection
#################################################################################
pollen_richness = left_join(asv_long, specs_withlandscape, by = "barcode_id") %>%
  mutate(OccurrenceTaxa = ifelse(ReadCountTaxa>0, 1, 0)) %>%
  filter(final_id == "B. mixtus")

pollen_wide = pivot_wider(pollen_richness, id_cols = barcode_id, names_from = genus, values_from = OccurrenceTaxa) %>%
  column_to_rownames("barcode_id")
pollen_wide[is.na(pollen_wide)] = 0
meta = pollen_richness %>% 
  group_by(barcode_id, sample_pt, site) %>%
  summarize(n=n()) %>%
  filter(!is.na(site)) %>%
  arrange(match(barcode_id, rownames(pollen_wide)))

transects = unique(meta$sample_pt)
sites <- unique(meta$site)

# Site-level accumulation curves — cleaner
accum_list <- lapply(sites, function(s){
  idx <- which(meta$site == s)
  site = meta$site[idx[1]]
  sub_mat <- pollen_wide[idx, ]
  sub_mat <- sub_mat[, colSums(sub_mat) > 0]
  
  sp <- specaccum(sub_mat, method = "rarefaction")
  
  data.frame(
    sample_pt    = s,
    site = site,
    bees    = sp$sites,
    richness = sp$richness,
    sd      = sp$sd
  )
})

landscapes = specs_withlandscape[c("sample_pt", "iji", "berry", "semi", "final_id")] %>%
  filter(final_id == "B. mixtus") %>%
  distinct()
accum_df <- bind_rows(accum_list)
accum_df = left_join(accum_df, landscapes, by = "sample_pt")

rarefaction = ggplot(accum_df, aes(x = bees, y = richness, colour = site, fill = site, group = sample_pt)) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), 
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Number of bees sampled",
       y = "Pollen taxa detected",
       title = "Sample-based rarefaction by site") +
  theme_bw() 
theme(legend.position= "none")
ggsave("7_figures/appendix_figures/site_level_rarefaction.png", rarefaction,
       height = 1000, width = 1500, units = "px")

#################################################################################
##### How is pollen GAMMA diversity related to landscape?
#################################################################################

pollen_points = pollen_richness %>%
  filter(final_id == "B. mixtus") %>%
  group_by(site, sample_pt, year) %>%
  summarize(num_genera = length(unique(genus)),
            num_bees = length(unique(barcode_id)),
            num_rounds = length(unique(round))) %>%
  left_join(landscapes)
goodcovg = pollen_points %>% filter(num_bees >= 5)
rarefactionbees = pollen_richness %>%
  ungroup() %>%
  filter(final_id == "B. mixtus") %>%
  distinct(barcode_id, sample_pt, year) %>%
  filter(paste0(sample_pt, year) %in% paste0(goodcovg$sample_pt, goodcovg$year))

# Get rarefied pollen richness for sites with 5+ bees
# (Rarify to 5 bees)

library(iNEXT)

transect_list <- rarefactionbees %>%
  group_by(sample_pt, year) %>%
  group_map(~ {
    sub_mat <- pollen_wide[.x$barcode_id, , drop = FALSE]
    
    freqs <- colSums(sub_mat)
    freqs <- freqs[freqs > 0]
    
    c(nrow(.x), freqs)
  }, .keep = TRUE) %>%
  setNames(
    rarefactionbees %>%
      group_by(sample_pt, year) %>%
      group_keys() %>%
      mutate(name = paste(sample_pt, year, sep = "_")) %>%
      pull(name)
  )
transect_list <- lapply(transect_list, unname)

# Run iNEXT
est_5 = estimateD(transect_list,
                  datatype = "incidence_freq",
                  base     = "size",
                  level = 5,
                  q = 0)

# Get standard error on each estimate
est_5 = est_5 %>%
  mutate(sample_pt = sub("_\\d{4}$", "", Assemblage)) %>%
  mutate(year = sub(".*_(\\d{4})$", "\\1", Assemblage)) %>%
  left_join(goodcovg)

# Fit model
fit = brm(
  qD  ~ berry + semi + iji + num_rounds + year + (1|site),
  data = est_5,
  family = lognormal(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))
fitdraws = as_draws_df(fit)
mean(fitdraws$b_berry > 0)


# Get effect of blueberry
berry_seq = seq(min(est_5$berry),
                max(est_5$berry),
                length.out = 50)

berrypred = lapply(berry_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- goodcovg |> mutate(berry = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(berry) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.berry = pretty(est_5$berry*10, n = 5)
axis.berry = labs.berry/10

# Plot
bberryplot = ggplot() +
  geom_line(data = berrypred, aes(x = berry, y = estimate), colour = faded_green, linewidth = 1) +
  geom_ribbon(data = berrypred, aes(x = berry, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_point(data = est_5, aes(x = berry, y = qD), colour = "black", alpha = 0.5) +
  xlab("Distance-weighted berry cultivation") +
  ylab("Rarified pollen richness (no. of genera)") +
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


year_diff = avg_predictions(fit,
                            variables = list(year = c(2022, 2023)),
                            newdata = est_5,
                            re_formula = NA)

year_diff_plot = ggplot(year_diff) +
  geom_rect(aes(xmin = as.numeric(factor(year)) - 0.3, 
                xmax = as.numeric(factor(year)) + 0.3, 
                ymin = conf.low, ymax = conf.high,
                fill = as.factor(year)), alpha = 0.4)+
  geom_segment(aes(x = as.numeric(factor(year)) - 0.3, 
                   xend = as.numeric(factor(year)) + 0.3, 
                   y = estimate, yend = estimate, 
                   colour = as.factor(year)), linewidth = 1) +
  geom_jitter(data = est_5, aes(x = as.numeric(factor(year)), y = qD, colour = as.factor(year)), 
              width = 0.1, alpha = 0.6, colour = "black") +
  scale_x_continuous(breaks = 1:2, 
                     labels = levels(factor(year_diff$year))) +
  scale_colour_manual(values = c(faded_green, faded_green)) +
  scale_fill_manual(values = c(faded_green, faded_green)) +
  labs(y = "Rarified pollen richness (no. of genera)", x = "Year") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11))

grid = bberryplot + year_diff_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave("7_figures/manuscript_figures/rarified_pollen_richness_mixtus.png", grid,
       width = 3000, height = 1500, units = "px")



#################################################################################
##### How is pollen GAMMA diversity related to time?
#################################################################################

pollen_timepoints = pollen_richness %>%
  filter(final_id == "B. mixtus") %>%
  group_by(site, round, year) %>%
  summarize(num_genera = length(unique(genus)),
            num_bees = length(unique(barcode_id)),
            num_points = length(unique(sample_pt)))
goodcovg = pollen_timepoints %>% filter(num_bees >= 5)
rarefactionbees = pollen_richness %>%
  ungroup() %>%
  filter(final_id == "B. mixtus") %>%
  distinct(barcode_id, round, site, year) %>%
  filter(paste0(round, site) %in% paste0(goodcovg$round, goodcovg$site))

# Get rarefied pollen richness for sites with 5+ bees
# (Rarify to 5 bees)

library(iNEXT)

SR_list <- rarefactionbees %>%
  group_by(round, site) %>%
  group_map(~ {
    sub_mat <- pollen_wide[.x$barcode_id, , drop = FALSE]
    
    freqs <- colSums(sub_mat)
    freqs <- freqs[freqs > 0]
    
    c(nrow(.x), freqs)
  }, .keep = TRUE) %>%
  setNames(
    rarefactionbees %>%
      group_by(round, site) %>%
      group_keys() %>%
      mutate(name = paste(round, site, sep = "_")) %>%
      pull(name)
  )
SR_list <- lapply(SR_list, unname)

# Run iNEXT
SR_est_5 = estimateD(SR_list,
                     datatype = "incidence_freq",
                     base     = "size",
                     level = 5,
                     q = 0)

# Join to coverage dataframe
datesframe = goodcovg %>% 
  left_join(dates) %>%
  group_by(site, round, year, num_points) %>%
  summarize(doy = mean(doy))

SR_est_5 = SR_est_5 %>%
  separate(Assemblage, into = c("round", "site"), sep = "_") %>%
  mutate(round = as.numeric(round),
         site = as.factor(site)) %>%
  left_join(datesframe)
SR_est_5$doy_centered = (SR_est_5$doy - 180)/10

# Fit model

fit2 = brm(
  qD  ~ year*doy_centered + year*I(doy_centered^2) + num_points + (1|site),
  data = SR_est_5,
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))
hypothesis(fit2, "Idoy_centeredE2 > 0")


# Get effect of blueberry
doy_seq = seq(min(SR_est_5$doy_centered),
              max(SR_est_5$doy_centered),
              length.out = 50)

doypred = lapply(doy_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- goodcovg |> mutate(doy_centered = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit2, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(doy_centered, year) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.doy = pretty(SR_est_5$doy_centered*10 + 180, n = 5)
axis.doy = (labs.doy-180)/10

# Plot
rarified_doy_plot = ggplot() +
  geom_line(data = doypred[doypred$year == 2022,], aes(x = doy_centered, y = estimate), colour = faded_green, linewidth = 1, linetype = "dashed") +
  geom_ribbon(data = doypred[doypred$year == 2022,], aes(x = doy_centered, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_line(data = doypred[doypred$year == 2023,], aes(x = doy_centered, y = estimate), colour = faded_green, linewidth = 1) +
  geom_ribbon(data = doypred[doypred$year == 2023,], aes(x = doy_centered, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_point(data = SR_est_5, aes(x = doy_centered, y = qD), colour = "black", alpha = 0.5) +
  xlab("Day of year") +
  ylab("") +
  scale_x_continuous(
    breaks = axis.doy,
    labels = labs.doy) +
  facet_grid(year~1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

grid = (bberryplot / (year_diff_plot + rarified_doy_plot)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave("7_figures/manuscript_figures/rarified_pollen_mixtus.png", grid,
       width = 2000, height = 3000, units = "px")




#################################################################################
##### What shapes pollen collection TURNOVER between individuals?
#################################################################################

# see forserver_bterry.R

###################################################################
## How does bee specialization change through time?
###################################################################

library(bipartite)

# Get bees from sites x timepoints with enough abundance
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
filteredspp = twospp %>%
  group_by(final_id, site, round) %>%
  filter(n() > 10)
# fix the active plant names to match the format of column names
filteredspp$active_flower = gsub(" spp\\.", ".spp.", filteredspp$active_flower)

# Get floral vegetation surveys to match
veg2022 = read.csv("3_data/cleandata/fielddata/2022vegetationdata.csv")
veg2023 = read.csv("3_data/cleandata/fielddata/2023vegetationdata.csv")

SRveg2022 = veg2022 %>%
  select(-X, -sample_id, -veg_observer) %>%
  group_by(site, round) %>%
  summarise(across(where(is.numeric), ~ sum(10^(.x - 1), na.rm = TRUE))) %>%
  unite("rowname", site, round, sep = "_") %>%
  column_to_rownames("rowname")
SRveg2023 = veg2023 %>%
  select(-X.5, -sample_id, -veg_observer) %>%
  group_by(site, round) %>%
  summarise(across(where(is.numeric), ~ sum(10^(.x - 1), na.rm = TRUE))) %>%
  unite("rowname", site, round, sep = "_") %>%
  column_to_rownames("rowname")
SRveg = bind_rows(SRveg2022, SRveg2023) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  select(where(~ sum(.x) > 0))


# network per round x site
specializations = distinct(filteredspp, final_id, round, site)
specializations$dprime = NA

for (i in 1:nrow(specializations)){
  # get plants
  plant_sub = SRveg[rownames(SRveg) == paste(specializations$site[i], specializations$round[i], sep = "_"),] %>%
    select(where(~ sum(.x) > 0))
  
  # get bees
  beesub = filteredspp %>% 
    filter(site == specializations$site[i]) %>%
    filter(round == specializations$round[i]) %>%
    filter(final_id == specializations$final_id[i])
    
  # get full plant lists
  net = t(table(beesub$active_flower, beesub$final_id))
  all_plants = union(colnames(net), colnames(plant_sub))
  
  #throw error if visited plant is not in plant abundance
  if(any(!colnames(net) %in% colnames(plant_sub))){
    print("ERROR: visited plant missing from vegetation surveys.")
    print(i)
  }
  
  # add missing plants to net table with zeros
  missing = setdiff(all_plants, colnames(net))
  net = cbind(net, matrix(0, nrow = nrow(net), 
                               ncol = length(missing),
                               dimnames = list(NULL, missing)))
  
  # reorder to match plant_sub column order
  net = net[, colnames(plant_sub), drop = FALSE]
  
  # make right format
  net_matrix = as.matrix(net)
  class(net_matrix) = "numeric"
  plant_vector = as.numeric(unlist(plant_sub))
  
  # now give them to bipartite for d'!
  specializations$dprime[i] = dfun(net_matrix, abuns = plant_vector)$dprime

}

# quickly visualize the specialization of the two species
ggplot(specializations, aes(x = dprime, colour = final_id, fill = final_id)) +
  geom_histogram() +
  scale_colour_manual(values = c(dark_gold, faded_green)) +
  scale_fill_manual(values = c(medium_gold, faded_medium))+
  facet_grid(final_id~1) +
  ylab("Count") +
  xlab("d'") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 11, face = "italic"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 11))

# get dates of reach round
sitedates = dates %>%
  group_by(site, round) %>%
  summarize(doy=mean(doy))
roundyears = twospp %>% distinct(round, year)


specializations = specializations %>%
  left_join(sitedates) %>%
  left_join(roundyears)
specializations$doy_centred = (specializations$doy - median(specializations$doy))/10
specializations$year = as.factor(specializations$year)


fit = brm(
  dprime ~ final_id*year*doy_centred + final_id*year*I(doy_centred^2) + (1|site),
  data = specializations,
  family = Beta(),
  chains = 4, cores = 4,
  init = 0,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99))

# test hypotheses for mixtus
hypothesis(fit, "doy_centred + final_idB.mixtus:doy_centred > 0") # mixtus doy 2022
hypothesis(fit, "doy_centred + final_idB.mixtus:doy_centred + year2023:doy_centred + final_idB.mixtus:year2023:doy_centred > 0") # mixtus doy 2023
hypothesis(fit, "Idoy_centredE2 + final_idB.mixtus:Idoy_centredE2 > 0") # mixtus doy^2 2022
hypothesis(fit, "Idoy_centredE2 + final_idB.mixtus:Idoy_centredE2 + year2023:Idoy_centredE2 + final_idB.mixtus:year2023:Idoy_centredE2 > 0") # mixtus doy^2 2023


# test hypotheses for impatiens
hypothesis(fit, "doy_centred > 0") # impatiens doy 2022
hypothesis(fit, "doy_centred + year2023:doy_centred > 0") # impatiens doy 2023
hypothesis(fit, "Idoy_centredE2 > 0") # impatiens doy^2 2022
hypothesis(fit, "Idoy_centredE2 + year2023:Idoy_centredE2 < 0") # impatiens doy^2 2023



# make some plots!
doy_seq <- seq(min(specializations$doy_centred),
               max(specializations$doy_centred),
               length.out = 50)

doypred <- lapply(doy_seq, function(d) {
  
  # Swap in this doy value, keep all other covariates real
  tmp <- specializations |> mutate(doy_centred = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp) |>
    as.data.frame() |>
    group_by(final_id, year, doy_centred) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows()

mixmin = min(specializations$doy_centred[specializations$final_id == "B. mixtus"])
mixmax = max(specializations$doy_centred[specializations$final_id == "B. mixtus" & specializations$year == 2023])
impmin = min(specializations$doy_centred[specializations$final_id == "B. impatiens"])
impmin22 = min(specializations$doy_centred[specializations$final_id == "B. impatiens" & specializations$year == 2022])
impmax = max(specializations$doy_centred[specializations$final_id == "B. impatiens"])

doypredfilt = doypred %>%
  filter((final_id == "B. mixtus" & doy_centred >= mixmin -0.5 & doy_centred <= mixmax +0.5) |
           (final_id == "B. impatiens" & doy_centred >= impmin -0.5 & doy_centred <= impmax + 0.5)) %>%
  filter(!(final_id == "B. impatiens" & year == 2022 & doy_centred < impmin22 - 0.5))


#create axis values for standardized variables
labs.doy = pretty(specializations$doy, n = 6)
axis.doy = (labs.doy-median(specializations$doy)) / 10

# Plot
visitation_dprime_plot = ggplot() +
  # plot impatiens predictions
  geom_line(data = doypredfilt[doypredfilt$final_id == "B. impatiens",], aes(x = doy_centred, y = estimate, colour = final_id), linewidth = 1) +
  geom_ribbon(data = doypredfilt[doypredfilt$final_id == "B. impatiens",], aes(x = doy_centred, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.5) +
  
  #plot mixtus predictions
  geom_line(data = doypredfilt[doypredfilt$final_id == "B. mixtus" & doypredfilt$year == 2023,], aes(x = doy_centred, y = estimate, colour = final_id), linewidth = 1) +
  geom_ribbon(data = doypredfilt[doypredfilt$final_id == "B. mixtus" & doypredfilt$year == 2023,], aes(x = doy_centred, ymin = conf.low, ymax = conf.high, fill = final_id), alpha = 0.5) +
  
  # plot data
  geom_point(data = specializations, aes(x = doy_centred, y = dprime), colour = "black", alpha = 0.6) +
  scale_colour_manual(values=c(gold, faded_medium)) +
  scale_fill_manual(values=c(gold, faded_medium)) +
  scale_x_continuous(
    breaks = axis.doy,
    labels = labs.doy) +
  xlab("Day of year") +
  facet_grid(year~final_id, scales = "free_x") +
  ylab("Specialization (d')") +
  theme_minimal() + 
  theme(legend.position = "none",
        strip.text.y = element_text(size = 14),
        strip.text.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

ggsave("7_figures/manuscript_figures/temporal_visitation_specialization.jpg", visitation_dprime_plot,
       width = 2000, height = 2000, units = "px")



####################################################################
##### How does plant abundance shift through time?
####################################################################

# Get sample effort
SReffort2022 = sample2022 %>% 
  group_by(site, round) %>%
  summarize(n = n())
SReffort2023 = sample2023 %>% 
  group_by(site, round) %>%
  summarize(n = n())

# Get site level richness, abundance, evenness
beeflowers = unique(c(specs2022$active_flower, specs2023$active_flower))
beeflowers = beeflowers[str_detect(beeflowers, "nest|N/A|flying|dipu|road|dead|grass|leaf|ground", negate = T)]

columnlist = c("site", "round", beeflowers)
floral_df_wide = veg2022[,colnames(veg2022) %in% columnlist]

floral_abundance1 = floral_df_wide %>%
  pivot_longer(!c(site, round), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(site, round, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(site, round, flower),
                ~ 10^(.-1))) %>%
  group_by(site, round) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE) + 1) %>%
  ungroup() %>%
  mutate(across(-c(site, round),
                ~ log10(.)))

floral_richness1 = floral_df_wide %>%
  pivot_longer(!c(site, round), names_to = "plant_code", values_to = "floral_abundance") %>%
  filter(!is.na(floral_abundance)) %>%
  group_by(site, round) %>%
  summarize(floral_richness = n())

floral_exponential1 = floral_df_wide %>%
  mutate(across(-c(site, round), as.numeric)) %>%
  mutate(across(-c(site, round), function(x) 10^(x-1))) %>%
  mutate(across(-c(site, round), ~replace(., is.na(.), 0))) %>%
  group_by(site, round) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  unite("rowname", site, round, sep = "_") %>%
  column_to_rownames("rowname")
  
floral_exponential1$floral_diversity = vegan::diversity(floral_exponential1)
floral_exponential1 = floral_exponential1 %>%
  rownames_to_column("rowname") %>%
  separate(rowname, into = c("site", "round"), sep = "_") %>%
  mutate(round = as.numeric(round))
floral_evenness1 = floral_exponential1[,colnames(floral_exponential1) %in% c("site", "round", "floral_diversity")] %>%
  left_join(floral_richness1, by = c("site", "round")) %>%
  mutate(even = floral_diversity/log(floral_richness))


floral_df1 = floral_abundance1 %>%
  left_join(floral_evenness1, by = c("site", "round")) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  pivot_longer(!c(site, round), names_to = "metric", values_to = "value") %>%
  left_join(sitedates) %>%
  mutate(year = 2022)

ggplot(floral_df1, aes(x = doy, y = value, group = metric)) +
  geom_point() +
  facet_grid(metric~1, scales = "free_y") +
  xlim(c(140,240)) +
  theme_bw()


# now for 2023
floral_df_wide = veg2023[,colnames(veg2023) %in% columnlist]

floral_abundance2 = floral_df_wide %>%
  pivot_longer(!c(site, round), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(site, round, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(site, round, flower),
                ~ 10^(.-1))) %>%
  group_by(site, round) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE) + 1) %>%
  ungroup() %>%
  mutate(across(-c(site, round),
                ~ log10(.)))

floral_richness2 = floral_df_wide %>%
  pivot_longer(!c(site, round), names_to = "plant_code", values_to = "floral_abundance") %>%
  filter(!is.na(floral_abundance)) %>%
  group_by(site, round) %>%
  summarize(floral_richness = n())

floral_exponential2 = floral_df_wide %>%
  mutate(across(-c(site, round), as.numeric)) %>%
  mutate(across(-c(site, round), function(x) 10^(x-1))) %>%
  mutate(across(-c(site, round), ~replace(., is.na(.), 0))) %>%
  group_by(site, round) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  unite("rowname", site, round, sep = "_") %>%
  column_to_rownames("rowname")

floral_exponential2$floral_diversity = vegan::diversity(floral_exponential2)
floral_exponential2 = floral_exponential2 %>%
  rownames_to_column("rowname") %>%
  separate(rowname, into = c("site", "round"), sep = "_") %>%
  mutate(round = as.numeric(round))
floral_evenness2 = floral_exponential2[,colnames(floral_exponential2) %in% c("site", "round", "floral_diversity")] %>%
  left_join(floral_richness2, by = c("site", "round")) %>%
  mutate(even = floral_diversity/log(floral_richness))


floral_df2 = floral_abundance2 %>%
  left_join(floral_evenness2, by = c("site", "round")) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  pivot_longer(!c(site, round), names_to = "metric", values_to = "value") %>%
  left_join(sitedates) %>%
  mutate(year = 2023)

floral_df = rbind(floral_df1, floral_df2)

floral_availability_plot = ggplot(floral_df, aes(x = doy, y = value, colour = site)) +
  geom_point() +
  facet_grid(metric~year, scales = "free_y",
             labeller = labeller(metric = c("even" = "Evenness", 
                                            "floral_abundance" = "Abundance", 
                                            "floral_diversity" = "Diversity", 
                                            "floral_richness" = "Richness"))) +
  xlab("Day of year") +
  ylab("Value") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", colour = "white"))

ggsave("7_figures/appendix_figures/season_floral_availability.png", floral_availability_plot,
       height = 3000, width = 2000, units = "px")


####################################################################
#### How does Rubus pollen collection change through time?
####################################################################

asvtable$barcode_id = rownames(asvtable)
asv_long = asvtable %>%
  pivot_longer(cols = -barcode_id, names_to = "Feature.ID", values_to = "ReadCount") %>%
  filter(ReadCount > 0) %>%
  group_by(barcode_id) %>%
  mutate(ReadProp = ReadCount/sum(ReadCount))%>%
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
  left_join(specs, by = "barcode_id") %>%
  left_join(dates, by = c("round", "sample_pt" = "sample_point", "site"), relationship = "many-to-one") %>%
  filter(!is.na(year)) %>%
  mutate(indicator = ifelse(genus == "Rubus", 1, 0)) %>%
  mutate(rubusspp = ifelse(indicator ==1, species, "other")) %>%
  filter(final_id == "B. mixtus")
asv_long$indicator[is.na(asv_long$indicator)] = 0

grouped_genus = asv_long %>%
  group_by(round, doy, year, barcode_id, site) %>%
  summarize(indicator = max(indicator)) %>%
  ungroup() %>%
  group_by(round, year, site) %>%
  summarize(doy = mean(doy),
            proprubus = sum(indicator)/n(),
            numbees = n())

rubus_poll_fig = ggplot(grouped_genus) +
  #geom_bar(aes(x = doy, y = proprubus*numbees, colour = site), stat = "identity") +
  geom_point(aes(x = doy, y = proprubus, colour = site, size = numbees)) +
  geom_line(aes(x = doy, y = proprubus, colour = site), linewidth = 1) +
  labs(x = "Day of year", y = "Proportion bees collecting Rubus spp. pollen", colour = "Site", size = "No. of bees", fill = NULL) +
  facet_grid(year~1) +
  scale_colour_viridis_d(option = "cividis") +
  theme_bw() +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_blank())
ggsave("7_figures/manuscript_figures/mixtus_rubus_consumption.png", rubus_poll_fig, 
       width = 18, height = 8, dpi = 400)


####################################################################
#### How does Rubus flower availability change through time?
####################################################################

ruarla_abundance1 = veg2022 %>%
  select(site, round, RULA, RUAR) %>%
  pivot_longer(!c(site, round), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(site, round, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(site, round, flower),
                ~ 10^(.-1))) %>%
  group_by(site, round) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE) + 1) %>%
  ungroup() %>%
  left_join(SReffort2022) %>%
  mutate(floral_abundance = floral_abundance/n) %>%
  mutate(across(-c(site, round),
                ~ log10(.))) %>%
  left_join(sitedates) %>%
  mutate(year = 2022)

ruarla_abundance2 = veg2023 %>%
  select(site, round, RULA, RUAR) %>%
  pivot_longer(!c(site, round), names_to = "flower", values_to = "floral_abundance") %>%
  mutate(across(-c(site, round, flower),
                ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(-c(site, round, flower),
                ~ 10^(.-1))) %>%
  group_by(site, round) %>%
  summarize(floral_abundance = sum(floral_abundance, na.rm = TRUE) + 1) %>%
  ungroup() %>%
  left_join(SReffort2023) %>%
  mutate(floral_abundance = floral_abundance/n) %>%
  mutate(across(-c(site, round),
                ~ log10(.))) %>%
  left_join(sitedates) %>%
  mutate(year = 2023)

ruarla_abundance = rbind(ruarla_abundance1, ruarla_abundance2)

ggplot() +
  # plot plant abundance data
  geom_point(data = ruarla_abundance, aes(x = doy, y = floral_abundance, colour = site)) +
  geom_line(data = ruarla_abundance, aes(x = doy, y = floral_abundance, colour = site), linewidth = 1) +
  labs(x = "Day of year", y = "Log-scaled abundance of Rubus armeniacus and Rubus lacerifolia", colour = "Site", size = "No. of bees", fill = NULL) +
  facet_grid(year~1) +
  scale_colour_viridis_d(option = "cividis") +
  theme_bw() +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text.x = element_blank())




####################################################################
#### Now get it all together!
####################################################################

SReffort = rbind(SReffort2022, SReffort2023)

# get mixtus phenology
mix = specs %>%
  filter(final_id == "B. mixtus") %>%
  group_by(site, round, year) %>%
  summarize(number_mixtus = n()) %>%
  left_join(SReffort) %>%
  mutate(number_mixtus = number_mixtus/n) %>%
  select(-n)

# get pollen richness
rich = mixtusrichness %>%
  group_by(site, round, year) %>%
  summarize(avg_richness = mean(pollrich)) %>%
  mutate(year = as.numeric(as.character((year))))

# get beta diversity
# Get jaccard distance between pollen loads
asv_long = asv_long %>%
  group_by(barcode_id, family, genus) %>%
  summarize(ReadCountTaxa = sum(ReadCount))

mixtus_matrix = pivot_wider(asv_long, id_cols = barcode_id, names_from = genus, values_from = ReadCountTaxa) %>%
  column_to_rownames("barcode_id") %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

jaccarddist = as.matrix(vegdist(mixtus_matrix, method ="jaccard"))

jaccarddf = as.matrix(jaccarddist) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "bee1") %>% 
  pivot_longer(cols = -bee1, names_to = "bee2", values_to = "jaccarddist") %>% 
  filter(bee1 > bee2)

# get some metadata we want to join
metadata1 = specs_withlandscape[c("round", "barcode_id", "site", "sample_pt", "sample_id", "year", "doy", "iji", "berry", "semi")]
colnames(metadata1) = c("round1", "bee1", "site1", "sample_pt1", "sample_id1", "year1", "doy1", "iji1", "berry1", "semi1")
metadata2 = specs_withlandscape[c("round", "barcode_id", "site", "sample_pt", "sample_id", "year", "doy", "iji", "berry", "semi")]
colnames(metadata2) = c("round2", "bee2", "site2", "sample_pt2", "sample_id2", "year2", "doy2", "iji2", "berry2", "semi2")

beta = jaccarddf %>%
  left_join(metadata1) %>%
  left_join(metadata2) %>%
  mutate(doydiff = abs(doy2-doy1)) %>%
  filter(sample_id1 == sample_id2) %>%
  group_by(site1, round1, year1) %>%
  summarize(mean_jacdist = mean(jaccarddist))
colnames(beta) = c("site", "round", "year", "mean_jacdist")
beta$year = as.numeric(as.character((beta$year)))

# get specialization
special = specializations %>%
  filter(final_id == "B. mixtus") %>%
  select(-doy_centred, -final_id) %>%
  mutate(year = as.numeric(as.character((year)))) %>%
  select(-doy)

# get rarefaction
rarefaction = SR_est_5 %>%
  distinct(site, round, year, qD) %>%
  mutate(year = as.numeric(as.character((year))))

# get all raw data
rawdata = ruarla_abundance %>%
  select(-doy, -n) %>%
  mutate(site = as.factor(site)) %>%
  left_join(grouped_genus) %>%
  left_join(mix) %>%
  left_join(rich) %>%
  left_join(beta) %>%
  left_join(special) %>%
  left_join(rarefaction)

rawdata_long = rawdata %>%
  pivot_longer(cols = c(floral_abundance, proprubus, number_mixtus, avg_richness, mean_jacdist, dprime, qD),
                            names_to = "metric",
                            values_to = "value")
  

rawdata_long$metric = factor(rawdata_long$metric, levels = c("number_mixtus", 
                                                             "floral_abundance",
                                                             "proprubus",
                                                             "avg_richness",
                                                             "mean_jacdist", 
                                                             "qD",
                                                             "dprime"
                                                             ))

comboplot = ggplot(rawdata_long) +
  #plot bee data
  geom_bar(aes(x = doy, y = value, colour = site, group = metric), 
           fill = "white", stat = "identity") +
  scale_colour_viridis_d(option = "cividis") +
  labs(x = "Day of year", y = "", colour = "Landscape", fill = "Landscape") +
  facet_grid(metric~year, scales = "free_y",
             labeller = labeller(metric = 
                                   c("number_mixtus" = "No. workers", 
                                     "floral_abundance" = "Rubus abundance",
                                     "proprubus" = "Rubus pollen use",
                                     "avg_richness" = "Pollen richness (ind)",
                                     "mean_jacdist" = "Pollen turnover", 
                                     "qD" = "Rarefied pollen richness",
                                     "dprime" = "Network specialization"
                                   ))) +
  geom_vline(xintercept = 170) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = 11),
        axis.title = element_text(size = 16))
  
ggsave("7_figures/manuscript_figures/temporal_plot.jpg", comboplot,
       width = 3000, height = 4000, units = "px")
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

splitbyseason = ggplot(byseason, aes(x = doy, y = pollrich, colour = season)) +
  geom_point(alpha = 0.4) +
  scale_colour_manual(values = c("lightgreen", "orange", "darkgreen")) +
  facet_grid(year~1) +
  xlab("Day of year") +
  ylab("Pollen richness (no. of genera per individual)") +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.background.y = element_rect(colour = "white", fill = "white")
  )

ggsave("7_figures/appendix_figures/samples_by_season.png", splitbyseason,
       height = 1500, width = 1500, units = "px")


#############################################################################
##### How does landscape impact pollen alpha diversity in different seasons?
#############################################################################

fit = brm(
  pollrich ~ year + season*iji + season*berry + season*semi + (1|site) + (1|sample_pt),
  data = byseason,
  family = negbinomial(),
  chains = 4, cores = 4,
  init = 0,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99)
)

         
# Get effect of IJI
iji_seq = seq(min(byseason$iji),
              max(byseason$iji),
              length.out = 50)

ijipred = lapply(iji_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- byseason |> mutate(iji = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(season, iji) |>
    summarise(estimate  = mean(estimate) + 1,
              conf.low  = mean(conf.low) + 1,
              conf.high = mean(conf.high) + 1,
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.iji = pretty(byseason$iji*10, n = 5)
axis.iji = labs.iji/10

# reorder factor levels for plotting
byseason = byseason %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))
ijipred = ijipred %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))

# Plot
ijiseasonplot = ggplot() +
  geom_line(data = ijipred[ijipred$season == "mid",], aes(x = iji, y = estimate), colour = "darkgreen", linewidth = 1) +
  geom_ribbon(data = ijipred[ijipred$season == "mid",], aes(x = iji, ymin = conf.low, ymax = conf.high), fill = "darkgreen", alpha = 0.3) +
  geom_point(data = byseason, aes(x = iji, y = pollrich + 1, colour = season), alpha = 0.5) +
  scale_x_continuous(
    breaks = axis.iji,
    labels = labs.iji) +
  scale_colour_manual(values = c("lightgreen", "darkgreen", "orange")) +
  scale_fill_manual(values = c("lightgreen", "darkgreen", "orange")) +
  xlab("Landscape interspersion") +
  ylab("") +
  facet_grid(season~1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))


# Get effect of IJI
semi_seq = seq(min(byseason$semi),
              max(byseason$semi),
              length.out = 50)

semipred = lapply(semi_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- byseason |> mutate(semi = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(season, semi) |>
    summarise(estimate  = mean(estimate) + 1,
              conf.low  = mean(conf.low) + 1,
              conf.high = mean(conf.high) + 1,
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.semi = pretty(byseason$semi*10, n = 5)
axis.semi = labs.semi/10

# reorder factor levels for plotting
byseason = byseason %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))
semipred = semipred %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))

# Plot
semiseasonplot = ggplot() +
  geom_line(data = semipred[semipred$season == "early",], aes(x = semi, y = estimate), colour = "darkgreen", linewidth = 1, linetype = "dashed") +
  geom_ribbon(data = semipred[semipred$season == "early",], aes(x = semi, ymin = conf.low, ymax = conf.high), fill = "lightgreen", alpha = 0.6) +
  geom_point(data = byseason, aes(x = semi, y = pollrich + 1, colour = season), alpha = 0.5) +
  scale_x_continuous(
    breaks = axis.semi,
    labels = labs.semi) +
  scale_colour_manual(values = c("lightgreen", "darkgreen", "orange")) +
  scale_fill_manual(values = c("lightgreen", "darkgreen", "orange")) +
  xlab("Distance-weighted seminatural area") +
  ylab("Pollen richness (no. of genera)") +
  facet_grid(season~1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))

grid = semiseasonplot + ijiseasonplot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave("7_figures/manuscript_figures/landscapexseason_alpha.png", grid,
       height = 3000, width = 2000, units = "px")


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
  geom_point(data = jacdf_meta, aes(x = berry1, y = jaccarddist), colour = "black", alpha = 0.5) +
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
ggplot() +
  geom_line(data = doypred[doypred$year1 == 2023,], aes(x = doy_centered, y = estimate), colour = faded_green, linewidth = 1) +
  geom_ribbon(data = doypred[doypred$year1 == 2023,], aes(x = doy_centered, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_point(data = jacdf_meta, aes(x = doy_centered, y = jaccarddist), colour = "black", alpha = 0.5) +
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




#################################################################################
##### How is pollen GAMMA diversity related to season x landscape?
#################################################################################

pollen_points = byseason %>%
  filter(final_id == "B. mixtus") %>%
  group_by(site, sample_pt, year, season) %>%
  summarize(num_bees = length(unique(barcode_id)),
            num_rounds = length(unique(round))) %>%
  left_join(landscapes)
goodcovg = pollen_points %>% filter(num_bees >= 5)
rarefactionbees = pollen_richness %>%
  ungroup() %>%
  left_join(byseason) %>%
  filter(final_id == "B. mixtus") %>%
  distinct(barcode_id, sample_pt, season, year) %>%
  filter(paste0(sample_pt, season, year) %in% paste0(goodcovg$sample_pt, goodcovg$season, goodcovg$year))

# Get rarefied pollen richness for sites with 5+ bees
# (Rarify to 5 bees)

library(iNEXT)

transect_list <- rarefactionbees %>%
  group_by(sample_pt, season, year) %>%
  group_map(~ {
    sub_mat <- pollen_wide[.x$barcode_id, , drop = FALSE]
    
    freqs <- colSums(sub_mat)
    freqs <- freqs[freqs > 0]
    
    c(nrow(.x), freqs)
  }, .keep = TRUE) %>%
  setNames(
    rarefactionbees %>%
      group_by(sample_pt, season, year) %>%
      group_keys() %>%
      mutate(name = paste(sample_pt, season, year, sep = "_")) %>%
      pull(name)
  )
transect_list <- lapply(transect_list, unname)

# Run iNEXT
est_5 = estimateD(transect_list,
                  datatype = "incidence_freq",
                  base     = "size",
                  level = 5,
                  q = 0)


est_5 = est_5 %>%
  separate(Assemblage, into = c("sample_pt", "season", "year"), sep = "_") %>%
  left_join(goodcovg)



# Fit model
fit = brm(
  qD  ~ season*berry + season*semi + season*iji + year + (1|site),
  data = est_5,
  family = Gamma(link="log"),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  init = 0,
  control = list(adapt_delta = 0.99))


fitdraws = as_draws_df(fit)
mean(fitdraws$b_berry > 0)
mean(fitdraws$b_berry + fitdraws$`b_seasonlate:berry` > 0)
mean(fitdraws$b_berry + fitdraws$`b_seasonmid:berry` > 0)

mean(fitdraws$b_semi > 0)
mean(fitdraws$b_semi + fitdraws$`b_seasonlate:semi` > 0)
mean(fitdraws$b_semi + fitdraws$`b_seasonmid:semi` > 0)

mean(fitdraws$b_iji > 0)
mean(fitdraws$b_iji + fitdraws$`b_seasonlate:iji` > 0)
mean(fitdraws$b_iji + fitdraws$`b_seasonmid:iji` > 0)


# Get effect of blueberry
berry_seq = seq(min(est_5$berry),
                max(est_5$berry),
                length.out = 50)

berrypred = lapply(berry_seq, function(d) {
  
  # Swap in this iji value, keep all other covariates real
  tmp <- goodcovg |> mutate(berry = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(berry, season) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.berry = pretty(est_5$berry*10, n = 5)
axis.berry = labs.berry/10

est_5 = est_5 %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))
berrypred = berrypred %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))

# Plot
berryplot = ggplot() +
  geom_line(data = berrypred[berrypred$season == "early",], aes(x = berry, y = estimate), colour = faded_green, linewidth = 1) +
  geom_ribbon(data = berrypred[berrypred$season == "early",], aes(x = berry, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_point(data = est_5, aes(x = berry, y = qD, colour = season)) +
  scale_colour_manual(values = c("lightgreen", "darkgreen", "orange")) +
  xlab("Distance-weighted berry cultivation") +
  ylab("Rarified pollen richness (no. of genera)") +
  scale_x_continuous(
    breaks = axis.berry,
    labels = labs.berry) +
  facet_grid(season~1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))


semi_seq = seq(min(est_5$semi),
                max(est_5$semi),
                length.out = 50)

semipred = lapply(semi_seq, function(d) {
  
  # Swap in this semi value, keep all other covariates real
  tmp <- goodcovg |> mutate(semi = d)
  
  # Predict, then immediately summarise — never accumulate full predictions
  predictions(fit, newdata = tmp, re_formula = NA) |>
    as.data.frame() |>
    group_by(semi, season) |>
    summarise(estimate  = mean(estimate),
              conf.low  = mean(conf.low),
              conf.high = mean(conf.high),
              .groups   = "drop")
}) |> bind_rows() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ .x + 1))

#create axis values for standardized variables
labs.semi = pretty(est_5$semi*10, n = 5)
axis.semi = labs.semi/10

est_5 = est_5 %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))
semipred = semipred %>%
  mutate(season = factor(season, levels = c("early", "mid", "late")))

# Plot
semiplot = ggplot() +
  geom_line(data = semipred[semipred$season == "early",], aes(x = semi, y = estimate), colour = faded_green, linewidth = 1, linetype = "dashed") +
  geom_ribbon(data = semipred[semipred$season == "early",], aes(x = semi, ymin = conf.low, ymax = conf.high), fill = faded_green, alpha = 0.3) +
  geom_point(data = est_5, aes(x = semi, y = qD, colour = season)) +
  scale_colour_manual(values = c("lightgreen", "darkgreen", "orange")) +
  xlab("Distance-weighted seminatural area") +
  ylab("") +
  scale_x_continuous(
    breaks = axis.berry,
    labels = labs.berry) +
  facet_grid(season~1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11))



grid = berryplot + semiplot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16))

ggsave("7_figures/manuscript_figures/rarified_richness_byseason_mixtus.png", grid,
       width = 2000, height = 2000, units = "px")






####################################################################
##### How is pollen collection COMPOSITION different between spp?
####################################################################

# Get bees x plant species matrix
y = as.matrix(comparison_matrix)
y[y>0] = 1
storage.mode(y) <- "integer"

# Get bee metadata
x_covars = specs_withlandscape %>%
  filter(barcode_id %in% rownames(y)) %>%
  select(barcode_id, final_id, round, iji, berry, semi, sample_pt, site)  %>%
  arrange(match(barcode_id, rownames(y)))
site_covars = x_covars[c("sample_pt", "site")]
site_covars$barcode_id = rownames(site_covars)

# HMSC needs the random effects as factors
XData <- x_covars %>%
  select(final_id, round) %>%
  mutate(
    final_id = as.factor(final_id),
    round = as.factor(round))

studyDesign <- data.frame(
  site      = as.factor(x_covars$site),
  sample_pt  = as.factor(x_covars$sample_pt),
  barcode_id       = as.factor(x_covars$barcode_id))

# Site-level random effect
rL.site <- HmscRandomLevel(units = levels(studyDesign$site))

# Transect-level random effect (nested within site)
rL.transect <- HmscRandomLevel(units = levels(studyDesign$sample_pt))

# Observation-level random effect
rL.bee <- HmscRandomLevel(units = levels(studyDesign$barcode_id))


# Prep for model
m <- Hmsc(
  Y             = y,
  XData         = XData,
  XFormula      = ~ final_id + round,   # fixed effects
  studyDesign   = studyDesign,
  ranLevels     = list(
    site      = rL.site,
    sample_pt  = rL.transect,
    barcode_id       = rL.bee
  ),
  distr = "probit"
)

# Fit with MCMC
# Start with short chains to check it runs, then scale up
thin    <- 10      # thinning interval
samples <- 1000     # posterior samples per chain
transient <- 500    # burn-in (in thinned units, so 500*thin actual steps)
nChains <- 4

m_fit <- sampleMcmc(
  m,
  thin       = thin,
  samples    = samples,
  transient  = transient,
  nChains    = nChains,
  nParallel  = 4,     # set to number of cores available
  verbose = TRUE)


# Convert to coda for standard MCMC diagnostics
mpost <- convertToCodaObject(m_fit)
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
gelman.diag(mpost$Lambda[[3]], multivariate = FALSE)$psrf  # [[3]] = bee level
effectiveSize(mpost$Beta)

# Posterior support for fixed effects
postBeta <- getPostEstimate(m_fit, parName = "Beta")

plotBeta(
  m_fit,
  post        = postBeta,
  param       = "Support",    # shows posterior support, not just mean
  supportLevel = 0.95,
  split        = 0.4,
  spNamesNumbers = c(FALSE, TRUE)
)

species_support <- postBeta$support[2, ]      # P(effect > 0)
species_support_neg <- postBeta$supportNeg[2, ] # P(effect < 0)

positive <- species_support[species_support > 0.95]      # overrepresented in sp2
negative <- species_support[species_support_neg > 0.95]  # overrepresented in sp1

positive
negative


# Extract latent factor scores for each bee -- the bee coordinates in "pollen composition space"
eta_samples <- lapply(m_fit$postList, function(chain) {
  lapply(chain, function(sample) sample$Eta[[3]])
})

# Average across all samples and chains
eta_mean <- Reduce("+", unlist(eta_samples, recursive = FALSE)) / 
  (length(m_fit$postList) * length(m_fit$postList[[1]]))

eta_df <- as.data.frame(eta_mean)
colnames(eta_df) <- paste0("factor", 1:5)
eta_df$bee_species <- XData$final_id

# Ordination plot
ggplot(eta_df, aes(x = factor1, y = factor2, colour = bee_species)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.95) +
  labs(x = "Latent factor 1", y = "Latent factor 2") +
  theme_bw()


####################################################################
##### PERMANOVA
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
  select(c(barcode_id, sample_pt, final_id, site, round)) %>%
  arrange(match(barcode_id, rownames(pollen_mat)))


dist_mat = vegdist(pollen_mat, method = "jaccard", binary = TRUE)

adonis2(dist_mat ~ round + site + final_id,
        data = meta,
        permutations = 999,
        by = "margin")

# PERMDISP: test dispersion differences between bee species
disp <- betadisper(dist_mat, group = meta$final_id)

# Quick omnibus permutation test
permutest(disp, permutations = 999)

# Visualise
## Option 1: PCoA ordination with hulls by species
plot(disp, hull = TRUE, label = FALSE,
     col = c("#E69F00", "#56B4E9"))  # colourblind-friendly

## Option 2: Boxplot of distances to group centroid (more interpretable)
disp_df <- data.frame(
  distance_to_centroid = disp$distances,
  bee_species          = meta$final_id,
  site                 = meta$site,
  round                = meta$round
)

ggplot(disp_df, aes(x = bee_species, y = distance_to_centroid, fill = bee_species)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  facet_wrap(~ site) +   # lets you see if pattern is consistent across sites
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  labs(y = "Distance to group centroid (Jaccard)",
       x = "Bee species",
       title = "Within-group dispersion by bee species") +
  theme_bw() +
  theme(legend.position = "none")

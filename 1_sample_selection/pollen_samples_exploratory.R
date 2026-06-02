setwd('/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/pollencode')
#rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")
load.packages()

#load survey data
specs2022 = read.csv("3_data/cleandata/fielddata/2022specimendata.csv")
specs2023 = read.csv("3_data/cleandata/fielddata/2023specimendata.csv")
sample2022 = read.csv("3_data/cleandata/fielddata/2022sampledata.csv")
sample2023 = read.csv("3_data/cleandata/fielddata/2023sampledata.csv")
vegData2022 = read.csv("3_data/cleandata/fielddata/2022vegetationdata.csv")
vegData2022[is.na(vegData2022)] <- 0
vegData2023 = read.csv("3_data/cleandata/fielddata/2023vegetationdata.csv")
vegData2023[is.na(vegData2023)] <- 0


#clean up: remove males / bees without final_id
allworkers2022 = specs2022 %>% 
  filter(final_id != "")%>%
  filter(str_detect(notes, "queen|male", negate = T))

allworkers2023 = specs2023 %>% 
  filter(final_id != "")%>%
  filter(str_detect(notes, "queen|male", negate = T))

#add date info to specimen dataframes
sample2022$date = paste(sample2022$day, sample2022$month, sample2022$year)
sample2022$date = gsub(" ", "", sample2022$date, fixed = TRUE)
sample2022$date <- as.POSIXlt(sample2022$date, format = "%d%b%y")
sample2022$julian_date = sample2022$date$yday

sample2023$date = paste(sample2023$day, sample2023$month, sample2023$year)
sample2023$date = gsub(" ", "", sample2023$date, fixed = TRUE)
sample2023$date <- as.POSIXlt(sample2023$date, format = "%d%b%y")
sample2023$julian_date = sample2023$date$yday

#join sample effort df to specimen df
allData2022 = left_join(allworkers2022, sample2022, by = c("sample_id", "site", "round", "year"))
allData2023 = left_join(allworkers2023, sample2023, by = c("sample_id", "site", "round", "year"))

#filter down to just mixtus and impatiens with pollen
allData2022 %>% 
  filter(pollen == "1" | pollen == "2" | pollen == "y" | pollen == "yy") %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") -> pollen2022
allData2023 %>% 
  filter(pollen == "1" | pollen == "2" | pollen == "y" | pollen == "yy") %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") -> pollen2023

#create summary df by subsite, month for each year
pollen2022 %>%
  group_by(final_id, site, month) %>%
  summarize(n = n()) %>%
  ungroup() -> summary2022
summary2022 %>% complete(final_id, site, month, fill = list(total = 0)) -> summarycomplete2022
summarycomplete2022$n <- ifelse(is.na(summarycomplete2022$n), 0, 
                                ifelse(summarycomplete2022$n > 10, 10,
                                       ifelse(summarycomplete2022$n < 3, 0, summarycomplete2022$n)))

summarycomplete2022$year = 2022

pollen2023 %>%
  filter(site != "") %>%
  group_by(final_id, site, month) %>%
  summarize(n = n()) %>%
  ungroup() -> summary2023
summary2023 %>% complete(final_id, site, month, fill = list(total = 0)) -> summarycomplete2023
summarycomplete2023$n <- ifelse(is.na(summarycomplete2023$n), 0, 
                                ifelse(summarycomplete2023$n > 10, 10,
                                       ifelse(summarycomplete2023$n < 3, 0, summarycomplete2023$n)))
summarycomplete2023$year = 2023

fullsummary = rbind(summarycomplete2022, summarycomplete2023)
fullsummary$monthyear = paste(fullsummary$month, fullsummary$year, sep="")
fullsummary = fullsummary %>% filter(round != 11 & round != 12)
fullsummary$monthyear <- factor(fullsummary$monthyear, levels = c("May2022", "June2022", "July2022", "August2022",
                                                                  "May2023", "June2023", "July2023", "August2023"))

fullmix = fullsummary[fullsummary$final_id == "B. mixtus",]
fullimp = fullsummary[fullsummary$final_id == "B. impatiens",]

#make heatmaps of which combos we have enough pollen for
ggplot(fullmix, aes(x = monthyear, y = site, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "orange") + # Color scale
  geom_vline(xintercept = 4.5, color = "black", linewidth = 1) +
  labs(title = "B. mixtus availability by subsite", x = "Timepoints", y = "subsite", fill = "Enough pollen?") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 8))

ggplot(fullimp, aes(x = monthyear, y = site, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") + # Color scale
  labs(title = "B. impatiens availability by subsite", x = "Timepoints", y = "subsite", fill = "Enough pollen?") +
  geom_vline(xintercept = 4.5, color = "black", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))


#plot overlap of mixtus and impatiens in each year


# Get capture rates per day
phenology2022 = allData2022 %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") %>%
  group_by(final_id, site, julian_date) %>%
  summarize(num = n())
phenology2023 = allData2023 %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") %>%
  group_by(final_id, site, julian_date) %>%
  summarize(num = n())

effort2022 = sample2022 %>%
  group_by(site, julian_date) %>%
  summarize(num_surveys = n())
effort2023 = sample2023 %>%
  group_by(site, julian_date) %>%
  summarize(num_surveys = n())

phenology2022 = phenology2022 %>%
  left_join(effort2022) %>%
  mutate(rate = num/num_surveys) %>%
  mutate(year = 2022)
phenology2023 = phenology2023 %>%
  left_join(effort2023) %>%
  mutate(rate = num/num_surveys) %>%
  mutate(year = 2023)

phenology = rbind(phenology2022, phenology2023)


phenologyfig = ggplot(phenology, aes(x = julian_date, y = rate, fill = site), colour = "black") +
  geom_col(position = "identity", alpha = 0.6) +
  facet_grid(final_id~year) +
  scale_fill_viridis_d(option = "turbo") +
  ylab("Captures per 5 minutes") +
  xlab ("Day of year") +
  labs(fill ="Landscape") +
  xlim(c(120, 235)) +
  theme_minimal()  +
  theme(strip.text.y = element_text(face = "italic"))
ggsave("7_figures/appendix_figures/speciesphenology.png",
       width = 3000, height = 1500, units = "px")




#plot mixtus vs impatiens abundance for 2022
specimenData2023_joined = left_join(specimenData2023, sample2023, by = "sample_id")
id_withpollen = specimenData2023 %>% filter(final_id == "B. mixtus") %>%
  filter(pollen == "1" | pollen == "2" | pollen =="y") %>%
  group_by(sample_id) %>%
  summarize(n=n()) %>%
  filter(n>=3)
impid_withpollen = specimenData2023 %>% filter(final_id == "B. impatiens") %>%
  filter(pollen == "1" | pollen == "2" | pollen =="y") %>%
  group_by(sample_id) %>%
  summarize(n=n()) %>%
  filter(n>=3)

specimenData2023_joined %>%
  filter(final_id != "") %>%
  count(sample_id, final_id) %>% # Count the number of individuals for each species in each sampling event
  pivot_wider(names_from = final_id, values_from = n, values_fill = 0) -> summary
  filter(subsite != "") -> summary # Ensure 0 counts for missing combinations

summary$round = str_split(summary$sample_id, "_", simplify = TRUE)[,1]
summary = summary[summary$`B. mixtus` > 0,]
summary = summary[summary$sample_id %in% id_withpollen$sample_id,]
summary = summary[summary$sample_id %in% impid_withpollen$sample_id,]

ggplot(summary, aes(x = `B. impatiens`, y = `B. mixtus`, color = round)) +
  geom_jitter() +
  scale_x_continuous() +
  theme_minimal()


################################################
#### run models for pollen presence / absence
################################################
specimenData2022 %>%
  filter(final_id == "B. mixtus" & pollen != "") -> mixtus2022
mixtus2022workers = mixtus2022[str_detect(mixtus2022$notes, "released|male|lost|queen", negate = T),]
mixtus2022workers$pollPres = ifelse(mixtus2022workers$pollen == "y" | 
                                      mixtus2022workers$pollen == "1" |
                                      mixtus2022workers$pollen =="2", 1, 0)

#standardize data in sampleEffort frame
sample2022$julian_date = scale(sample2022$julian_date, center = TRUE, scale = TRUE)
landscapeMetrics$prop_blueberry = scale(landscapeMetrics$prop_blueberry, center = TRUE, scale = TRUE)
landscapeMetrics$prop_edge = scale(landscapeMetrics$prop_edge, center = TRUE, scale = TRUE)
landscapeMetrics$landscape_shdi = scale(landscapeMetrics$landscape_shdi, center = TRUE, scale = TRUE)

mixtus2022frame = mixtus2022workers %>%
  left_join(sample2022 %>% select(-c("site", "round", "sample_pt", "year")), by = "sample_id")
mixtus2022frame = merge(mixtus2022frame, landscapeMetrics, by = "sample_pt")
mixtus2022frame %>% group_by(sample_pt) %>% summarize(n=n()) -> out

#%>% filter(n > 4) -> out
mixtus2022thinned = mixtus2022frame[mixtus2022frame$sample_pt %in% out$sample_pt,]

formula.pollpres <- formula(pollPres ~
                               month*prop_blueberry +
                              month*prop_edge +
                              month*landscape_shdi +
                               (1|sample_pt) + (1|site)
)

bf.poll <- bf(formula.pollpres, family="bernoulli")

ncores = 4
## run model without interactions
fit.poll.pres <- brm(bf.poll, mixtus2022thinned,
                           cores=ncores,
                           iter = (10^4),
                           chains =1,
                           thin=1,
                           init=0,
                           save_pars = save_pars(all = TRUE),
                           open_progress = FALSE,
                           control = list(adapt_delta = 0.999,
                                          stepsize = 0.001,
                                          max_treedepth = 20)
)

summary(fit.poll.pres)
#####
#check model diagnostics
######

run_plot_freq_model_diagnostics(formula.pollpres,
                                this_data=mixtus2022thinned,
                                this_family="bernoulli", site.lat="site")


########
# make some plots
#########
all.cond.effects <- conditional_effects(fit.poll.pres)

pollblue <-
  all.cond.effects[["landscape_shdi:month"]]

#ggplot
poll.blue <- 
  
  #plot model prediction with credible interval
  ggplot(pollblue, aes(x = landscape_shdi,
                   y = estimate__), color = month) +
  geom_line(aes(x = landscape_shdi, y=estimate__, color = month)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  alpha=0.3, fill = month)) +
  
  #plot raw data
  geom_jitter(data=mixtus2022thinned,
             aes(x=landscape_shdi, y=pollPres, color = month), cex=2, alpha = 0.5, height = 0.02) +
  labs(x = "Landscape diversity", y = "Pollen presence",
       fill = "") +
  theme_ms() +
  theme(legend.position = "bottom") +
  guides(color = "none", alpha = "none") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))
poll.blue

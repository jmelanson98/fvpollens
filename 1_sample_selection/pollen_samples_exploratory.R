setwd('/Users/jenna1/Documents/UBC/Bombus Project/FVBombus_code/pollencode')
#rm(list=ls())
source("src/ggplotThemes.R")
source("src/init.R")
source("src/misc.R")
load.packages()

#load survey data
specimenData = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022specimendata.csv", sep = ",", header = T))
specimenData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023specimendata.csv", sep = ",", header = T, fill = TRUE))
samplePoints = as.data.frame(read.table("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/allsamplepoints.csv", sep = ","))
colnames(samplePoints) = c("sample_pt", "gps","landowner","subsite")
sampleEffort2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022sampledata.csv", sep = ","))
sampleEffort2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023sampledata.csv", sep = ","))
vegData2022 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2022vegetationdata.csv", sep = ","))
vegData2022[is.na(vegData2022)] <- 0
vegData2023 = as.data.frame(read.csv("/Users/jenna1/Documents/UBC/Bombus Project/Raw Data/2023vegetationdata.csv", sep = ","))
vegData2023[is.na(vegData2023)] <- 0


#clean up: remove males / bees without final_id
specimenData %>% filter(final_id != "") -> allbees
allworkers = allbees[str_detect(allbees$notes, "released|male|lost", negate = T),]

specimenData2023 %>% filter(final_id != "") -> allbees2023
allworkers2023 = allbees2023[str_detect(allbees2023$notes, "released|male|lost", negate = T),]

#add subsite + date info to specimen dataframes
#wrangle gps coordinates into shape
samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
samplePoints = samplePoints[,c("sample_pt", "subsite", "lat", "long")]

#join sample points to sample effort data
colnames(sampleEffort2022)[colnames(sampleEffort2022) == "sample_point"] = "sample_pt"
sampleEffort2022 = left_join(sampleEffort2022, samplePoints, by = "sample_pt")

colnames(sampleEffort2023)[colnames(sampleEffort2023) == "sample_point"] = "sample_pt"
sampleEffort2023 = left_join(sampleEffort2023, samplePoints, by = "sample_pt")

#add julian date to sample effort data frame
sampleEffort2022$date = paste(sampleEffort2022$day, sampleEffort2022$month, sampleEffort2022$year)
sampleEffort2022$date = gsub(" ", "", sampleEffort2022$date, fixed = TRUE)
sampleEffort2022$date <- as.POSIXlt(sampleEffort2022$date, format = "%d%b%y")
sampleEffort2022$julian_date = sampleEffort2022$date$yday

sampleEffort2023$date = paste(sampleEffort2023$day, sampleEffort2023$month, sampleEffort2023$year)
sampleEffort2023$date = gsub(" ", "", sampleEffort2023$date, fixed = TRUE)
sampleEffort2023$date <- as.POSIXlt(sampleEffort2023$date, format = "%d%b%y")
sampleEffort2023$julian_date = sampleEffort2023$date$yday

#join sample effort df to specimen df
allData2022 = left_join(allworkers, sampleEffort2022, by = "sample_id")
allData2023 = left_join(allworkers2023, sampleEffort2023, by = "sample_id")

#filter down to just mixtus and impatiens with pollen
allData2022 %>% 
  filter(pollen == "1" | pollen == "2" | pollen == "y") %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") -> pollen2022
allData2023 %>% 
  filter(pollen == "1" | pollen == "2" | pollen == "y") %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") -> pollen2023

#create summary df by subsite, month for each year
pollen2022 %>%
  filter(subsite != "") %>%
  group_by(final_id, subsite, month) %>%
  summarize(n = n()) %>%
  ungroup() -> summary2022
summary2022 %>% complete(final_id, subsite, month, fill = list(total = 0)) -> summarycomplete2022
summarycomplete2022$n <- ifelse(is.na(summarycomplete2022$n), 0, 
                                ifelse(summarycomplete2022$n > 10, 10,
                                       ifelse(summarycomplete2022$n < 3, 0, summarycomplete2022$n)))

summarycomplete2022$year = 2022

pollen2023 %>%
  filter(subsite != "") %>%
  group_by(final_id, subsite, month) %>%
  summarize(n = n()) %>%
  ungroup() -> summary2023
summary2023 %>% complete(final_id, subsite, month, fill = list(total = 0)) -> summarycomplete2023
summarycomplete2023$n <- ifelse(is.na(summarycomplete2023$n), 0, 
                                ifelse(summarycomplete2023$n > 10, 10,
                                       ifelse(summarycomplete2023$n < 3, 0, summarycomplete2023$n)))
summarycomplete2023$year = 2023

fullsummary = rbind(summarycomplete2022, summarycomplete2023)
fullsummary$monthyear = paste(fullsummary$month, fullsummary$year, sep="")
fullsummary = filter(fullsummary, !(round.x %in% c(11, 12)))
fullsummary$monthyear <- factor(fullsummary$monthyear, levels = c("May2022", "June2022", "July2022", "August2022",
                                                                  "May2023", "June2023", "July2023", "August2023"))

fullmix = fullsummary[fullsummary$final_id == "B. mixtus",]
fullimp = fullsummary[fullsummary$final_id == "B. impatiens",]

#make heatmaps of which combos we have enough pollen for
ggplot(fullmix, aes(x = monthyear, y = subsite, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "orange") + # Color scale
  geom_vline(xintercept = 10.5, color = "black", linewidth = 1) +
  labs(title = "B. mixtus availability by subsite", x = "Timepoints", y = "subsite", fill = "Enough pollen?") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 8))

ggplot(fullimp, aes(x = monthyear, y = subsite, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") + # Color scale
  labs(title = "B. impatiens availability by subsite", x = "Timepoints", y = "subsite", fill = "Enough pollen?") +
  geom_vline(xintercept = 4.5, color = "black", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))

#potential sample selection
pollen2022 %>% 
  filter(subsite != "") %>%
  group_by(subsite, month, final_id) %>%
  filter(n() >= 5) %>%
  slice_sample(n=5) -> subsiteselected2022
pollen2023 %>% 
  filter(subsite != "") %>%
  group_by(subsite, month, final_id) %>%
  filter(n() >= 5) %>%
  slice_sample(n=5) -> subsiteselected2023

pollen2022 %>%
  group_by(site.x, month, final_id) %>%
  filter(n() >= 5) %>%
  slice_sample(n=5) -> siteselected2022
pollen2023 %>%
  group_by(site.x, month, final_id) %>%
  filter(n() >= 5) %>%
  slice_sample(n=5) -> siteselected2023

# pollen2022 %>%
#   group_by(sample_pt.x, month, final_id) %>%
#   filter(n() >= 3) %>%
#   slice_sample(n=5) -> transectselected2022
# pollen2023 %>%
#   group_by(sample_pt.x, month, final_id) %>%
#   filter(n() >= 3) %>%
#   slice_sample(n=5) -> transectselected2023

all2022 = rbind(subsiteselected2022, siteselected2022)
all2023 = rbind(subsiteselected2023, siteselected2023)

unique2022 = all2022 %>%
  distinct() 
unique2023 = all2023 %>%
  distinct()


  group_by(subsite, month, final_id) %>%
  summarize(n = n_distinct(sample_pt.x)) -> subsitesum23
  
unique2022 %>%
  group_by(final_id, sample_pt.x) %>%
  summarize(n=n()) -> out

dim(unique2022 %>%
  filter(final_id == "B. impatiens"))











#repeat for mixtus 2023
#clean up: remove males / bees without final_id
specimenData2023 %>% filter(final_id != "" & barcode_id != "NA") %>%
  filter(final_id == "B. impatiens") -> cleanedmixtus
cleanedmixtus = cleanedmixtus[str_detect(cleanedmixtus$notes, "released|male|lost", negate = T),]

#join sample points to sample effort data
colnames(sampleEffort2023)[colnames(sampleEffort2023) == "sample_point"] = "sample_pt"
sampleEffort2023 = left_join(sampleEffort2023, samplePoints, by = "sample_pt")

#add julian date to sample effort data frame
sampleEffort2023$date = paste(sampleEffort2023$day, sampleEffort2023$month, sampleEffort2023$year)
sampleEffort2023$date = gsub(" ", "", sampleEffort2023$date, fixed = TRUE)
sampleEffort2023$date <- as.POSIXlt(sampleEffort2023$date, format = "%d%b%y")
sampleEffort2023$julian_date = sampleEffort2023$date$yday

#join sample effort df to specimen df
allData2023 = left_join(cleanedmixtus, sampleEffort2023, by = "sample_id")




allData2023 %>%
  filter(pollen != "n" & pollen != "0") %>%
  group_by(month, subsite) %>%
  summarize(n=n()) %>%
  filter(n >= 4 & subsite != "") -> summary

#4 specimens * 45 subsite x months = 180 specimens
# 4 specimens * 50 subsite x months = 200 specimens
# mixtus: 320 (for landscape stuff)
# impatiens: 360 (for landscape stuff)

#680 out of 760 or 1140 --> so we can use another 80 or 460, depending on # of runs


#plot overlap of mixtus and impatiens in each year
#clean up: remove males / bees without final_id
specimenData %>% filter(final_id != "") %>%
  filter(final_id == "B. impatiens" | final_id == "B. mixtus") -> cleanedspecs
cleanedspecs = cleanedspecs[str_detect(cleanedspecs$notes, "released|male|lost", negate = T),]
allData2022 = left_join(cleanedspecs, sampleEffort2022, by = "sample_id")

specimenData2023 %>% filter(final_id != "") %>%
  filter(final_id == "B. impatiens" | final_id == "B. mixtus") -> cleanedspecs2023
cleanedspecs2023 = cleanedspecs2023[str_detect(cleanedspecs2023$notes, "released|male|lost", negate = T),]
allData2023 = left_join(cleanedspecs2023, sampleEffort2023, by = "sample_id")


#make some plots
allData2022 %>%
  group_by(final_id, julian_date) %>%
  summarize(num = n()) -> phenology2022

mix2022all = allData2022[allData2022$final_id == "B. mixtus",]
imp2022all = allData2022[allData2022$final_id == "B. impatiens",]
mix2022 = phenology2022[phenology2022$final_id == "B. mixtus",]
imp2022 = phenology2022[phenology2022$final_id == "B. impatiens",]

allData2023 %>%
  group_by(final_id, julian_date) %>%
  summarize(num = n()) -> phenology2023

mix2023all = allData2023[allData2023$final_id == "B. mixtus",]
imp2023all = allData2023[allData2023$final_id == "B. impatiens",]
mix2023 = phenology2023[phenology2023$final_id == "B. mixtus",]
imp2023 = phenology2023[phenology2023$final_id == "B. impatiens",]

x_limits <- range(c(81, 235))

# Create the first histogram
hist1 <- ggplot(mix2022all, aes(x = julian_date, fill = site.x)) +
  geom_histogram(binwidth = 0.5, position = "stack", alpha = 0.7) +
  scale_x_continuous(limits = x_limits) +
  theme_minimal() +
  labs(title = "Mixtus 2022")

# Create the second histogram
hist2 <- ggplot(imp2022all, aes(x = julian_date, fill = site.x)) +
  geom_histogram(binwidth = 0.5, position = "stack", alpha = 0.7) +
  scale_x_continuous(limits = x_limits) +
  theme_minimal() +
  labs(title = "Impatiens 2022")

# Create the third histogram
hist3 <- ggplot(mix2023all, aes(x = julian_date, fill = site.x)) +
  geom_histogram(binwidth = 0.5, position = "stack", alpha = 0.7) +
  scale_x_continuous(limits = x_limits) +
  theme_minimal() +
  labs(title = "Mixtus 2023")

# Create the fourth histogram
hist4 <- ggplot(imp2023all, aes(x = julian_date, fill = site.x)) +
  geom_histogram(binwidth = 0.5, position = "stack", alpha = 0.7) +
  scale_x_continuous(limits = x_limits) +
  theme_minimal() +
  labs(title = "Impatiens 2023")

# Stack the plots vertically, ensuring alignment
final_plot <- (hist1 | hist3) / (hist2 | hist4)
print(final_plot)      


allData2022 %>%
  filter(pollen != "n" & pollen != "0") %>%
  group_by(sample_pt.x, final_id, round.x) %>%
  summarize(n=n()) %>%
  filter(n>=5) -> out



#plot mixtus vs impatiens abundance for 2022
specimenData2023_joined = left_join(specimenData2023, sampleEffort2023, by = "sample_id")
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
sampleEffort2022$julian_date = scale(sampleEffort2022$julian_date, center = TRUE, scale = TRUE)
landscapeMetrics$prop_blueberry = scale(landscapeMetrics$prop_blueberry, center = TRUE, scale = TRUE)
landscapeMetrics$prop_edge = scale(landscapeMetrics$prop_edge, center = TRUE, scale = TRUE)
landscapeMetrics$landscape_shdi = scale(landscapeMetrics$landscape_shdi, center = TRUE, scale = TRUE)

mixtus2022frame = mixtus2022workers %>%
  left_join(sampleEffort2022 %>% select(-c("site", "round", "sample_pt", "year")), by = "sample_id")
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




# Load in packages
library(biomformat)
library(qiime2R)


#################################################
### Load in data
#################################################
# Get specimen metadata data
specs2022 = read.csv("3_data/cleandata/fielddata/2022specimendata.csv")
specs2023 = read.csv("3_data/cleandata/fielddata/2023specimendata.csv")
sample2022 = read.csv("3_data/cleandata/fielddata/2022sampledata.csv")
sample2023 = read.csv("3_data/cleandata/fielddata/2023sampledata.csv")

specs2022filt = specs2022[c("barcode_id", "active_flower", "final_id", "site", "round", "sample_pt", "sample_id", "year")]
specs2023filt = specs2023[c("barcode_id", "active_flower", "final_id", "site", "round", "sample_pt", "sample_id", "year")]
specs = rbind(specs2022filt, specs2023filt)

# Get raw seqs
raw_seqs = as.character(read_qza("3_data/merged_ASVs.qza")$data)

# Read in merged OTU tables and classifications
otu = read_qza("3_data/merged_qiimeseqtab.qza")$data
classifications = read_qza("3_data/merged_taxonomy.qza")$data


#################################################
### Clean ASVs with low resolution
#################################################

# Get low resolution taxa
lowres = classifications %>% filter(Taxon == "Streptophyta")
otu_lowres = otu[rownames(otu) %in% lowres$Feature.ID,]
print(sum(otu_lowres)/sum(otu)) # low resolution sequences make up about 7% of reads
# blast searching says these are genuinely not identifiable (not just PlantITS gaps)
# potentially some chimeras that made it through dada2?

# Remove low resolution ASVs
classifications_highres = classifications %>% filter(Taxon != "Streptophyta")
otu_highres = otu[rownames(otu) %in% classifications_highres$Feature.ID,]

# Tranpose otu table and pivot longer for plotting
otut = as.data.frame(t(otu_highres))
otut$barcode_id = rownames(otut)
otu_long = pivot_longer(otut, cols = -barcode_id, names_to = "Feature.ID", values_to = "ReadCount")
otu_long = otu_long %>%
  group_by(barcode_id) %>%
  mutate(ReadProp = ReadCount/sum(ReadCount))
otu_taxa = left_join(otu_long, classifications_highres, by = "Feature.ID")
otu_taxa = otu_taxa %>%
  separate(Taxon, 
           into = c("phylum", "na_col", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right")

samples = otu_taxa %>% filter(ReadProp > 0)
samples$scientificname = paste(samples$genus, samples$species, sep = " ")
samples = left_join(samples, specs, by = "barcode_id")


#################################################
### Visualize negative controls
#################################################
negative_controls = samples %>% filter(barcode_id %in% c("PN01", "PN02", "PN03", "PN04",
                                                         "PN05", "PN06", "PN07", "PN08",
                                                         "PN09", "PN10", "PN11", "PN12")) %>%
  filter(ReadCount > 0) %>%
  mutate(indicator = ifelse(Feature.ID == "ASV1", 1, 0))

ggplot(negative_controls, aes(x = barcode_id, y = ReadCount, fill = genus)) +
  geom_col(width = 0.9, colour = "white", linewidth = 0.1) +
  #scale_fill_manual(values = shuffled_colors) +
  scale_fill_viridis_d(option = "turbo") +
  labs(x = "Sample", y = "Read count", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        #legend.key.size = unit(0.1, "cm"),
        #legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  guides(fill = guide_legend(ncol = 10))


#################################################
### Remove contaminant ASVs from samples
#################################################
otut = as.data.frame(t(otu_highres))

# remove reads with constitute < 0.5% of the sample (following https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13780)
rel_abund = otut / rowSums(otut)
otu_sper = otut
otu_sper[rel_abund < 0.005] = 0
  
# remove reads which are less than max contamination (following https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13780)
max_contam = negative_controls %>%
  group_by(Feature.ID) %>%
  summarize(maxReadCount = max(ReadCount))
otu_decon = otu_sper

for (ASV in max_contam$Feature.ID){
  otu_decon[ASV][otu_decon[ASV] < max_contam$maxReadCount[max_contam$Feature.ID == ASV]] = 0
}
sum(otu_decon)


# first remove samples with < 1000 high quality filtered reads
otu_decon = otu_decon[rowSums(otu_decon) > 1000,]


#################################################
### Visualize samples
#################################################

otu_decon$barcode_id = rownames(otu_decon)
otu_long = pivot_longer(otu_decon, cols = -barcode_id, names_to = "Feature.ID", values_to = "ReadCount")
otu_long = otu_long %>%
  group_by(barcode_id) %>%
  mutate(ReadProp = ReadCount/sum(ReadCount))
otu_final = left_join(otu_long, classifications_highres, by = "Feature.ID")
otu_final = otu_final %>%
  separate(Taxon, 
           into = c("phylum", "na_col", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right") %>%
  filter(ReadProp > 0)
otu_final$scientificname = paste(otu_final$genus, otu_final$species, sep = " ")
otu_final = left_join(otu_final, specs, by = "barcode_id")

samples2022 = otu_final %>% filter(year == 2022) %>%
  filter(final_id == "B. mixtus" | final_id == "B. impatiens") %>%
  filter(ReadCount > 0)

temp2022 = ggplot(samples2022, aes(x = barcode_id, y = ReadProp, fill = genus)) +
  geom_col(width = 0.9, colour = "white", linewidth = 0.1) +
  #scale_fill_manual(values = shuffled_colors) +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  facet_grid(final_id ~ site, scales = "free", space = "free_x") +
  labs(x = "Sample", y = "Relative abundance", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  guides(fill = guide_legend(ncol = 20))
ggsave("6_figures/2022_pollen_consumption.png", temp2022, width = 14, height = 8, dpi = 300)



barcode_order <- otu_final %>%
  filter(!is.na(year)) %>%
  filter(final_id == "B. mixtus") %>%
  filter(ReadCount > 0) %>%
  mutate(round = as.numeric(round)) %>%
  distinct(barcode_id, round) %>%
  arrange(round, barcode_id) %>%
  pull(barcode_id)

samplesmixtus = otu_final %>% 
  filter(!is.na(year)) %>%
  filter(final_id == "B. mixtus") %>%
  filter(ReadCount > 0) %>%
  mutate(indicator = ifelse(genus == "Vaccinium", 1, 0)) %>%
  mutate(barcode_id = factor(barcode_id, levels = barcode_order))
samplesmixtus$indicator[is.na(samplesmixtus$indicator)] = 0

blueberrymixtus = ggplot(samplesmixtus, aes(x = barcode_id, y = ReadProp, fill = as.factor(indicator))) +
  geom_col(width = 0.9, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("beige", "blue")) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  facet_grid(year ~ site, scales = "free", space = "free_x") +
  labs(x = "Sample", y = "Relative abundance", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  guides(fill = guide_legend(ncol = 20))
ggsave("6_figures/mixtus_blueberry_consumption.png", blueberrymixtus, width = 18, height = 8, dpi = 400)

tempmixtus = ggplot(samplesmixtus, aes(x = barcode_id, y = ReadProp, fill = scientificname)) +
  geom_col(width = 0.9, colour = "white", linewidth = 0.1) +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  facet_grid(year ~ site, scales = "free", space = "free_x") +
  labs(x = "Sample", y = "Relative abundance", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  guides(fill = guide_legend(ncol = 20))
ggsave("6_figures/mixtus_temp_consumption.png", tempmixtus, width = 18, height = 8, dpi = 400)

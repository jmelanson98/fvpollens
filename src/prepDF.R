prepSpecimenTable <- function(specimen.data, parasite.scores){
  #this function will:
  #### - add parasite scores to the specimen dataframe and clean unwanted individuals
  
  #join dataframes to get specimenData and parasite scores into a single dataframe
  parasiteScores = mutate(parasiteScores, across(-c(barcode_id), as.numeric))
  allspecimens = left_join(specimenData, parasiteScores, by = "barcode_id")
  
  #clean up: remove males / bees without final_id
  allspecimens %>% filter(final_id != "") %>%
    filter(final_id != "B. huntii") -> cleanedspecimens
  cleanedspecimens = cleanedspecimens[str_detect(cleanedspecimens$notes, "released|male|lost", negate = T),]
  
  return(cleanedspecimens)
  
}


calculateLandscapeMetrics <- function(landcover.raster, site.shapefile, landcover.classification){
  #this function will:
  #### - calculate area of different land classes in 500m buffer around each sample_pt
  #### - calculate landscape shdi in 500m buffer around each sample_pt
  #### - return a dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt
  
  #set and change CRS
  crs(landcover.raster) <- 900913
  st_crs(site.shapefile) <- 900913

  #area of each class in each buffer zone
  res1_lsm_c_ca <- sample_lsm(landcover.raster,
                              y = site.shapefile,
                              plot_id = site.shapefile$site_id,
                              shape = "circle",
                              size = 500,
                              return_raster = TRUE,
                              what = "lsm_c_ca")
  #set class IDs
  res1_lsm_c_ca = merge(res1_lsm_c_ca, landcover.classification[,c("lc_type", "class")], by = "class")
  colnames(res1_lsm_c_ca)[colnames(res1_lsm_c_ca) == "plot_id"] <- "sample_pt"
  
  #shannon diversity
  res1_lsm_l_shdi <- sample_lsm(landcover.raster,
                              y = site.shapefile,
                              plot_id = site.shapefile$site_id,
                              shape = "circle",
                              size = 500,
                              return_raster = TRUE,
                              what = "lsm_l_shdi")
  colnames(res1_lsm_l_shdi)[colnames(res1_lsm_l_shdi) == "plot_id"] <- "sample_pt"
  colnames(res1_lsm_l_shdi)[colnames(res1_lsm_l_shdi) == "value"] <- "landscape_shdi"
  
  #calculate % landscape blueberry
  res1_lsm_c_ca %>% 
    group_by(sample_pt) %>%
    summarise(buffer_area=sum(value)) -> buffer_areas
  res1_lsm_c_ca = merge(res1_lsm_c_ca, buffer_areas, by = "sample_pt")
  res1_blueberry = filter(res1_lsm_c_ca, lc_type == "blueberry")
  res1_blueberry$prop_blueberry = res1_blueberry$value/res1_blueberry$buffer_area
  
  #calculate landscape edge area (hedgerow + blackberry + low edge)
  res1_edge = filter(res1_lsm_c_ca, lc_type == "blackberry" | lc_type == "lowedge" | lc_type == "hedgerow")
  res1_edge$split_prop = res1_edge$value/res1_edge$buffer_area
  res1_edge %>%
    group_by(sample_pt) %>%
    summarize(prop_edge = sum(split_prop)) -> res1_edge_summary
  
  #create dataframe containing prop_blueberry, prop_edge, and shdi for each sample_pt
  landscape.dat = full_join(res1_blueberry[,c("sample_pt", "prop_blueberry")], res1_edge_summary, by = "sample_pt")
  landscape.dat$prop_blueberry[is.na(landscape.dat$prop_blueberry)] <- 0
  landscape.dat = full_join(landscape.dat, res1_lsm_l_shdi[,c("sample_pt", "landscape_shdi")], by = "sample_pt")
  landscape.dat$prop_edge[is.na(landscape.dat$prop_edge)] <- 0
  
  vars = c("prop_blueberry", "prop_edge", "landscape_shdi")
  for(var in vars){
    hist(landscape.dat[,var], main=var)
  }
  
  return(landscape.dat)
}




prepBRMSdf <- function(sampleEffort, samplePoints, cleanedSpecimens, flowerList, vegData, landscapeData){
  #this function will add the following to the sample effort dataframe:
  #### - add GPS coords
  #### - julian date
  #### - floral/bee diversity/abundance data
  #### - landscape metrics (prop_blueberry, prop_edge, shdi)
  # it will then merge the sample effort and specimen data frames, and remove columns that are not
  # used by brms
  # The output is a dataframe where every row is an observation of a single bee, or an "empty row" for a sampling
  #event where no bees were captured. Suitable for bernoulli family parasite models.
  
  #first, wrangle gps coordinates into shape
  samplePoints$lat = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) tail(x, 1))
  samplePoints$long = rapply(strsplit(gsub("\\(|\\).*", "", samplePoints$gps), split = " "), function(x) head(x, 3)[2])
  samplePoints = samplePoints[,c("sample_pt", "subsite", "lat", "long")]
  
  #join sample points to sample effort data
  colnames(sampleEffort)[colnames(sampleEffort) == "sample_point"] = "sample_pt"
  sampleEffort = left_join(sampleEffort, samplePoints, by = "sample_pt")
  
  #add julian date to sample effort data frame
  sampleEffort$date = paste(sampleEffort$day, sampleEffort$month, sampleEffort$year)
  sampleEffort$date = gsub(" ", "", sampleEffort$date, fixed = TRUE)
  sampleEffort$date <- as.POSIXlt(sampleEffort$date, format = "%d%b%y")
  sampleEffort$julian_date = sampleEffort$date$yday
  
  #add bumblebee diversity and abundance to sample effort dataframe
  community_df_long = group_by(cleanedSpecimens, sample_id) %>%
    count(final_id)
  community_df_wide = as.data.frame(pivot_wider(community_df_long, names_from = "final_id", values_from = "n"))
  community_df_wide[is.na(community_df_wide)] <- 0
  rownames(community_df_wide) = community_df_wide[,1]
  community_df = community_df_wide[,-1]
  community_df_wide$H <- vegan::diversity(community_df)
  colnames(community_df_wide) = c("sample_id", "impatiens_abundance", "rufocinctus_abundance", 
                                  "californicus_abundance", "vos_abundance", "flavifrons_abundance", 
                                  "mixtus_abundance", "melanopygus_abundance", "sitkensis_abundance", "bombus_shannon_diversity")
  community_df_wide$bombus_abundance = community_df_wide$impatiens_abundance + 
    community_df_wide$melanopygus_abundance + community_df_wide$mixtus_abundance + 
    community_df_wide$flavifrons_abundance + community_df_wide$rufocinctus_abundance + 
    community_df_wide$californicus_abundance + community_df_wide$vos_abundance
  
  sampleEffort = left_join(sampleEffort, community_df_wide, by = "sample_id")
  sampleEffort[is.na(sampleEffort)] <- 0
  
  #prepare veg data frame
  colnames(vegData)[3] = "sample_pt"
  
  #calculate floral diversity and add to sample effort dataframe
  columnlist = c("sample_id", flowerList)
  floral_df_wide = vegData[,colnames(vegData) %in% columnlist]
  floral_df_wide %>%
    mutate(floral_df_wide, across(-c(sample_id), as.numeric)) %>%
    mutate(across(-c(sample_id), function(x) 10^x)) %>%
    mutate(across(-c(sample_id), ~replace(., is.na(.), 0))) -> floral_df_wide_exponential 
  rownames(floral_df_wide_exponential) = floral_df_wide_exponential[,1]
  floral_df_wide_exponential = floral_df_wide_exponential[,-1]
  floral_df_wide_exponential$floral_diversity <- vegan::diversity(floral_df_wide_exponential)
  floral_df_wide_exponential$sample_id = rownames(floral_df_wide_exponential)
  floral_df = floral_df_wide_exponential[,colnames(floral_df_wide_exponential) %in% c("sample_id", "floral_diversity")]
  sampleEffort = left_join(sampleEffort, floral_df, by = "sample_id")
  
  #calculate floral abundance and add to sample effort dataframe
  floral_df_wide %>%
    pivot_longer(!sample_id, names_to = "flower", values_to = "floral_abundance") -> floral_df_long
  floral_df_long$floral_abundance = as.numeric(floral_df_long$floral_abundance)
  floral_df_long$floral_abundance = 10^(floral_df_long$floral_abundance - 1)
  abundance_aggregate = aggregate(floral_abundance ~ sample_id, floral_df_long, sum)
  abundance_aggregate$floral_abundance = log10(abundance_aggregate$floral_abundance + 1)
  sampleEffort = left_join(sampleEffort, abundance_aggregate, by = "sample_id")
  
  #join landscape data to sample effort dataframe
  sampleEffort = left_join(sampleEffort, landscapeData, by = "sample_pt")
  
  #join sampleEffort dataframe to specimens dataframe
  #check that empty rows (e.g., sampling effort where no bees were captured) are maintained
  brmsdf = full_join(sampleEffort, cleanedSpecimens, by = 'sample_id')
  brmsdf <- brmsdf %>%
    select(-site.y, -round.y, -sample_pt.y, -year.y) %>%
    rename(site = site.x,
           round = round.x,
           sample_pt = sample_pt.x,
           year = year.x,
           sample_notes = notes.x,
           specimen_notes = notes.y)
  
  #create a column for all crithidia
  brmsdf$hascrithidia <- brmsdf$cbombii + brmsdf$cexpoeki + brmsdf$crithidiaspp
  brmsdf$hascrithidia[brmsdf$hascrithidia>0] <- 1 
  
  #create column for nosema
  brmsdf$hasnosema <- brmsdf$nbombii + brmsdf$nceranae
  brmsdf$hasnosema[brmsdf$hasnosema>0] <- 1 
  
  #create a column for all parasites
  brmsdf$any_parasite = brmsdf$nbombii + brmsdf$nceranae +
    brmsdf$crithidiaspp + brmsdf$cbombii + brmsdf$cexpoeki +
    brmsdf$apicystis + brmsdf$ascosphaera
  brmsdf$any_parasite[brmsdf$any_parasite > 0] <- 1
  
  #create a column for native bee abundance
  brmsdf$native_bee_abundance = brmsdf$bombus_abundance - brmsdf$impatiens_abundance
  print(head(brmsdf))
  
  #remove columns that are not used in brms
  colsToKeep = c("site", "round", "sample_pt", "sample_id", "barcode_id", "caste", "sampledImp", "subsite", "lat", "long", "julian_date", "floral_abundance", "floral_diversity", "impatiens_abundance",
                 "bombus_shannon_diversity", "native_bee_abundance", "prop_blueberry", "prop_edge", "landscape_shdi",
                 "final_id", "apidae", "apicystis", "ascosphaera", "cbombii", "crithidiaspp", "cexpoeki", "hascrithidia",
                 "nbombii", "nceranae", "hasnosema", "any_parasite")
  brmsdf_reduced = brmsdf[,colnames(brmsdf) %in% colsToKeep]
  write.csv(brmsdf_reduced, "fvimp_brmsdf.csv")
  return(brmsdf_reduced )
  
}







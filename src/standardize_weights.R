standardize <- function(x)
  ## subtract mean and devide by sd to standardize data
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)


makeDataMultiLevel <- function(indiv.data, site.col, year.col){
  ## Sets weights to 0s and 1s so that a single dataset can be passed
  ## to brms but site-level data isn't duplicated by the number of
  ## specimens

  ## split data by year
  indiv.data.split <- split(indiv.data, indiv.data[, year.col])

  out.indiv.data <- lapply(indiv.data.split, addWeightCol,
                           site.col=site.col)

  out.indiv.data <- do.call(rbind, out.indiv.data)

  return(out.indiv.data)

}


addWeightCol <- function(each.year.dat, site.col){
  ## for all individuals from the same site, 
  site.ids <- unlist(tapply(each.year.dat[, site.col],
                            each.year.dat[, site.col],
                            function(x) 1:length(x)))


  names(site.ids) <- NULL
  each.year.dat$SiteIDs <- site.ids
  each.year.dat$Weights <- each.year.dat$SiteIDs
  each.year.dat$Weights[each.year.dat$Weights > 1] <- 0
  each.year.dat$Subset <- each.year.dat$Weights == 1
  return(each.year.dat)
}

standardizeVars <- function(spec.data, vars, key = NULL, by.site=TRUE){
  ##  center all of the x variables, need to use unique values to avoid
  ##  repetition by the number of specimens
  if(by.site){
    unique.site.vals <-  unique(spec.data[,c("sample_pt", key, vars)])
  } else {
    unique.site.vals <-  unique(spec.data[,c(key, vars)])
  }
  unique.site.vals[, vars] <- apply(as.matrix(unique.site.vals[, vars]), 2, standardize)
  print("Dimensions of the data before merging the standardize data")
  print(dim(spec.data))
  spec.data[, vars] <- NULL
  spec.data <- merge(spec.data, unique.site.vals, all.x=TRUE)
  print("Dimensions of the data after merging the standardize data")
  print(dim(spec.data))
  layout(matrix(1:(2*round(length(vars)/2)), nrow=2))
  for(var in vars){
    hist(unique.site.vals[, var], main=var)
  }
  
  return(spec.data)
}

prepParasiteWeights <- function(spec.data){
  ## create a dummy varaible "WeightPar" for the parasite data. The
  ## intention is to keep stan from dropping data for site-level models,
  ## but weight is 0 for parasite models.
  spec.data$WeightsPar <- 1
  spec.data$WeightsPar[spec.data$apidae == 0 |
                       is.na(spec.data$apidae)] <- 0
  
  ## stan drops all NA data, so can set AnyParasite to 0 with WeightsPar
  ## to keep it in the models
  spec.data$any_parasite[is.na(spec.data$any_parasite)] <- 0
  spec.data$cbombii[is.na(spec.data$cbombii)] <- 0
  spec.data$nbombii[is.na(spec.data$nbombii)] <- 0
  spec.data$hasnosema[is.na(spec.data$hasnosema)] <- 0
  spec.data$crithidiaspp[is.na(spec.data$crithidiaspp)] <- 0
  spec.data$hascrithidia[is.na(spec.data$hascrithidia)] <- 0
  spec.data$apicystis[is.na(spec.data$apicystis)] <- 0
  spec.data$final_id = as.character(spec.data$final_id)
  spec.data$final_id[is.na(spec.data$final_id)] <- "B. medius"
  spec.data$final_id = sapply(strsplit(spec.data$final_id, ". "), function(x){paste("Bombus", x[[2]], sep="_")})
  spec.data$caste[is.na(spec.data$caste)] <- "worker"
  
  return(spec.data)
}


prepDataSEM_negbin <-
  function(spec.data,#individual level specimen data
           variables.to.log = NULL, #variables to be logged
           variables.to.log.1 = NULL, #variables to be logged + 1
           vars_transect,#variables to standardize at transect level
           vars_transect_sr,#variables to standardize at transect x sampling round level
           vars_sp = NULL,  #variables to standardize at the species level
           vars_sp_yearsr= NULL) #variables to standardize at the species/year/sr level
{
  ## Function for making the SEM weights and standarizing variables.
  spec.data <- spec.data[order(spec.data$sample_pt), ]
  
  ## create a dummy variable "Weight" to deal with the data sets being at
  ## different levels to get around the issue of having to pass in one
  ## data set into brms

  ## species, year, SR key
  spec.data$round_species <-
    paste(spec.data$round, spec.data$final_id, sep = ";")
  
  print("Number of unique site, sampling round combinations")
  print(length(unique(paste(spec.data$sample_pt, spec.data$round))))
  spec.data <- makeDataMultiLevel(spec.data, "sample_pt", "round")
  print("Number of individuals with Weights == 1, should be the same as above")
  print(sum(spec.data$Weights))

  ## log variables before standardizing
  if(!is.null(variables.to.log)){
    spec.data[, variables.to.log] <-
      log(spec.data[, variables.to.log])
  }
  if(!is.null(variables.to.log.1)){
    spec.data[, variables.to.log.1] <-
      log(spec.data[, variables.to.log.1]+ 1)
  }
  
  ##  center all of the x variables, need to use unique values to avoid
  ##  repetition by the number of specimens
  
  if(!is.null(vars_transect_sr)){
    print("Standardizing variables with site, sampling round combinations")
    spec.data <- standardizeVars(spec.data, vars_transect_sr, "round")
  }
  if(!is.null(vars_sp)){
    print("Standardizing variables by species")
    spec.data <-
      standardizeVars(spec.data, vars_sp, "final_id", by.stand=FALSE)
  }
  if(!is.null(vars_transect)){
    print("Standardizing variables with site combinations")
    spec.data <- standardizeVars(spec.data, vars_transect)
  }

  if(!is.null(vars_sp_yearsr)){
    print("Standardizing variables with sampling round, site, species combinations")
    spec.data <- standardizeVars(spec.data, vars_sp_yearsr, "round_species")
  }
  
  ## create a dummy varaible "WeightPar" for the parasite data. The
  ## original intention was to keep stan from dropping data for
  ## site-level models, but weight is 0 for parasite models.
  
  #print("Number of successful parasite screenings")
  #print(sum(spec.data$apidae, na.rm = TRUE))
  #spec.data <- prepParasiteWeights(spec.data)
  #print("Number of of individuals with WeightsPar == 1, should be the same as above")
  #print(sum(spec.data$WeightsPar))
  #print("Final dim of data after adding WeightsPar")
  #print(dim(spec.data))
  
  ##moved the following code here from prepParasiteWeights because the prepParasiteWeights
  ##function is not necessary for the neg binomial dataset (I set weightsPar already)
  
  spec.data$spp_anyparasite[is.na(spec.data$spp_anyparasite)] <- 0
  spec.data$spp_cbombi[is.na(spec.data$spp_cbombi)] <- 0
  spec.data$spp_nbombi[is.na(spec.data$spp_nbombi)] <- 0
  spec.data$spp_anycrithidia[is.na(spec.data$spp_anycrithidia)] <- 0
  spec.data$spp_apicystis[is.na(spec.data$spp_apicystis)] <- 0
  spec.data$number_screened[is.na(spec.data$number_screened)] <- 0
  spec.data$weightsPar[is.na(spec.data$weightsPar)] <- 0
  spec.data$final_id[is.na(spec.data$final_id)] <- "none"
  #spec.data$final_id = as.factor(spec.data$final_id)
  spec.data$spp_crithidia[is.na(spec.data$spp_crithidia)] <- 0
  spec.data$subsetPar = spec.data$weightsPar > 0
  
  return(spec.data)
}

prepDataSEM_bernoulli <-
  function(spec.data,#individual level specimen data
           variables.to.log = NULL, #variables to be logged
           variables.to.log.1 = NULL, #variables to be logged + 1
           vars_transect,#variables to standardize at transect level
           vars_transect_sr,#variables to standardize at transect x sampling round level
           vars_sp = NULL,  #variables to standardize at the species level
           vars_sp_yearsr= NULL) #variables to standardize at the species/year/sr level
  {
    ## Function for making the SEM weights and standarizing variables.
    spec.data <- spec.data[order(spec.data$sample_pt), ]
    
    ## create a dummy variable "Weight" to deal with the data sets being at
    ## different levels to get around the issue of having to pass in one
    ## data set into brms
    
    ## species, year, SR key
    spec.data$round_species <-
      paste(spec.data$round, spec.data$final_id, sep = ";")
    
    print("Number of unique site, sampling round combinations")
    print(length(unique(paste(spec.data$sample_pt, spec.data$round))))
    spec.data <- makeDataMultiLevel(spec.data, "sample_pt", "round")
    print("Number of individuals with Weights == 1, should be the same as above")
    print(sum(spec.data$Weights))
    
    ## log variables before standardizing
    if(!is.null(variables.to.log)){
      spec.data[, variables.to.log] <-
        log(spec.data[, variables.to.log])
    }
    if(!is.null(variables.to.log.1)){
      spec.data[, variables.to.log.1] <-
        log(spec.data[, variables.to.log.1]+ 1)
    }
    
    ##  center all of the x variables, need to use unique values to avoid
    ##  repetition by the number of specimens
    
    if(!is.null(vars_transect_sr)){
      print("Standardizing variables with site, sampling round combinations")
      spec.data <- standardizeVars(spec.data, vars_transect_sr, "round")
    }
    if(!is.null(vars_sp)){
      print("Standardizing variables by species")
      spec.data <-
        standardizeVars(spec.data, vars_sp, "final_id", by.stand=FALSE)
    }
    if(!is.null(vars_transect)){
      print("Standardizing variables with site combinations")
      spec.data <- standardizeVars(spec.data, vars_transect)
    }
    
    if(!is.null(vars_sp_yearsr)){
      print("Standardizing variables with sampling round, site, species combinations")
      spec.data <- standardizeVars(spec.data, vars_sp_yearsr, "round_species")
    }
    
    ## create a dummy varaible "WeightPar" for the parasite data. The
    ## original intention was to keep stan from dropping data for
    ## site-level models, but weight is 0 for parasite models.
    
    print("Number of successful parasite screenings")
    print(sum(spec.data$apidae, na.rm = TRUE))
    spec.data <- prepParasiteWeights(spec.data)
    print("Number of of individuals with WeightsPar == 1, should be the same as above")
    print(sum(spec.data$WeightsPar))
    print("Final dim of data after adding WeightsPar")
    print(dim(spec.data))
    
    #convert Weights/WeightPar columns into Subset/SubsetPar columns
    spec.data$subsetPar = spec.data$WeightsPar > 0
    spec.data$Subset = spec.data$Weights > 0
    spec.data$impSubset = spec.data$Subset & spec.data$sampledImp #same as subset, but rows where we did not sample for impatiens are now FALSE
    
    #add column for native/non-native status
    spec.data$status = ifelse(spec.data$final_id == "Bombus_impatiens" | spec.data$final_id == "Bombus_medius", "nonnative", "native")
    
    #drop all rows where there are NA's
    spec.data$floral_abundance[is.na(spec.data$floral_abundance)] = 0
    spec.data = filter(spec.data, !is.na(floral_diversity))
    return(spec.data)
  }
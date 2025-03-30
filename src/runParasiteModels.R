
runParasiteModels <- function(spec.data,
                              parasite,
                              xvars){

  formula.parasite  <- as.formula(paste(
    paste(parasite, "| subset(subsetPar)"),  #+ trials(1)
    paste(xvars,
          collapse=" + "),
    sep=" ~ "))

  return(formula.parasite)
}



runParasiteModelsNoSubset <- function(spec.data,
                              parasite,
                              xvars){
  
  formula.parasite  <- as.formula(paste(
    paste(parasite),  #+ trials(1)
    paste(xvars,
          collapse=" + "),
    sep=" ~ "))
  
  return(formula.parasite)
}
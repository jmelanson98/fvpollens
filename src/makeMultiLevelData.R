
makeDataMultiLevel <- function(indiv.data, site.col, yr.col){
    ## split data by year
    indiv.data.split <- split(indiv.data, indiv.data[, yr.col])

    ## maybe in the future this will need to be an sapply
    out.indiv.data <- lapply(indiv.data.split, addWeightCol,
                             site.col=site.col)

    out.indiv.data <- do.call(rbind, out.indiv.data)

    return(out.indiv.data)

}


addWeightCol <- function(each.year.dat, site.col){
    site.ids <- unlist(tapply(each.year.dat[, site.col],
                              each.year.dat[, site.col],
                              function(x) 1:length(x)))


    names(site.ids) <- NULL
    each.year.dat$SiteIDs <- site.ids
    each.year.dat$Weights <- each.year.dat$SiteIDs
    each.year.dat$Weights[each.year.dat$Weights > 1] <- 0
    return(each.year.dat)
}

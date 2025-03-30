## models
library(brms)
library(bayestestR)

## general
library(tidyverse)

## plotting
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(ggthemes)
library(ggtext)
library(stringr)
library(viridis)
library(gridExtra)
library(gridExtra)
library(cowplot)
library(modelr)

## trees
library(phytools)
library(ggtree)
library(ape)

save.dir0 <- "saved/"
if(!dir.exists(save.dir0)) {
  dir.create(save.dir0, showWarnings = FALSE)
}

save.dir <- "saved/tables"
if(!dir.exists(save.dir)) {
  dir.create(save.dir, showWarnings = FALSE)
}

fig.dir <- "figures"
if(!dir.exists(fig.dir)) {
  dir.create(fig.dir, showWarnings = FALSE)
}


fig.dir2 <- "figures/diagnostics"
if(!dir.exists(fig.dir2)) {
  dir.create(fig.dir2, showWarnings = FALSE)
}

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
load.packages <- function(){
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("here",
                "readr",
                "readxl",
                "writexl", 
                "openxlsx",
                "janitor",
                "filesstrings",
                "textshaping",
                "extrafont",
                "qrcode", 
                "lubridate", 
                "naniar", 
                "Hmisc",
                "sf",
                "maps",
                "broom",
                #"rgdal",
                "raster",
                "dismo",
                "terra",
                "spData",
                "tmap",
                "leaflet",
                #"rgeos",
                #"maptools",
                "mapview",
                "stars",
                "landscapemetrics",
                "ggrepel",
                "vctrs",
                "measurements", 
                "arsenal", 
                "vegan",
                "taxize",
                "iNEXT", 
                "tidyverse",
                "stringr",
                "ggplot2", 
                "purrr",
                "magrittr",
                "dplyr",
                "tidylog",
                "data.table",
                "diffr",
                "RColorBrewer")
  
  ipak(packages)
}

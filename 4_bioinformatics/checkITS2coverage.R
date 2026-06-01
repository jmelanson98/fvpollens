### Check ITS2 database coverage
# June 1 2026
# J Melanson


# Load packages
library(taxize)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(rgbif)

bombus_path = "/Users/jenna1/Documents/UBC/bombus_project/"


# Load in plant list data
observedplants = read.csv("3_data/cleandata/plantdata/observedplants.csv")
databaseplants = read.csv("3_data/ITS2_taxonomy", header = FALSE, sep = ";")
colnames(databaseplants) = c("ID", "phylum", "class", "family", "genus", "species")


# Get classifications
observedplants = observedplants[c("scientific_name", "common_name", "plant_code")]

# GNA verifier function (safe version)
gna_verifier_safe = function(names, batch_size = 50) {
  do.call(rbind, lapply(seq(1, length(names), by = batch_size), function(i) {
    gna_verifier(names[i:min(i + batch_size - 1, length(names))])
  }))
}


# Get canonical names
resolved = gna_verifier_safe(observedplants$scientific_name)

resolved_df = tibble(
  input_plants     = observedplants$scientific_name,
  canonicalname    = resolved$matchedCanonicalFull)
resolved_df = resolved_df[!is.na(resolved_df$canonicalname), ]

# Get full classification from GBIF
get_classification = function(input_name, canonical_name) {
  tryCatch({
    
    match = name_backbone(name = canonical_name)
    
    # handle NONE match by falling back to name_suggest
    if (match$matchType == "NONE") {
      print(paste0(canonical_name, " not found. Trying name lookup."))
      suggestion = name_lookup(q = canonical_name, 
                               higherTaxonKey = 6, limit = 1)$data
      if (is.null(suggestion) || nrow(suggestion) == 0) stop("No match")
      key = suggestion$key[1]
      rank_out = "GENUS"
      canonical_out = canonical_name  # just use what we passed in
      rank_out     = suggestion$rank[1]
      status_out   = NA_character_
      authority_out = NA_character_
    } else {
      key           = match$usageKey
      canonical_out = match$canonicalName
      rank_out      = match$rank
      status_out    = match$status
      authority_out = match$authorship
    }
    
    if (is.null(key) || is.na(key)) stop("No usable key")
    
    cls = classification(key, db = "gbif")[[1]]
    if (is.null(cls) || nrow(cls) == 0) stop("No classification returned")
    
    ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    tax = setNames(rep(NA_character_, length(ranks)), ranks)
    for (r in ranks) {
      val = cls$name[cls$rank == r]
      if (length(val) > 0) tax[r] = val
    }
    
    tibble(
      input_plants  = input_name,
      canonicalname = canonical_out,
      taxa_rank     = rank_out,
      taxa_status   = status_out,
      synonym_name  = ifelse(!is.na(status_out) & status_out == "SYNONYM", 
                             canonical_name, NA_character_),
      authority     = authority_out,
      kingdom       = tax["kingdom"],
      phylum        = tax["phylum"],
      class         = tax["class"],
      order         = tax["order"],
      family        = tax["family"],
      genus         = tax["genus"],
      species       = tax["species"]
    )
    
  }, error = function(e) {
    message("Failed for: ", canonical_name, " — ", e$message)
    tibble(
      input_plants  = input_name,
      canonicalname = NA_character_,
      taxa_rank     = NA_character_,
      taxa_status   = NA_character_,
      synonym_name  = NA_character_,
      authority     = NA_character_,
      kingdom       = NA_character_,
      phylum        = NA_character_,
      class         = NA_character_,
      order         = NA_character_,
      family        = NA_character_,
      genus         = NA_character_,
      species       = NA_character_
    )
  })
}


# map to input plants
taxonomy_df = bind_rows(
  Map(get_classification, resolved_df$input_plants, resolved_df$canonicalname))

# add tayberry because it's a weirdo
tayberry_row = tibble(
  input_plants  = "Rubus fruticosus x Rubus idaeus",
  canonicalname = "Rubus fruticosus x Rubus idaeus",  # hybrid cross, not really recognized
  taxa_rank     = "UNRANKED",
  taxa_status   = "UNRESOLVED",
  synonym_name  = NA,
  kingdom       = "Plantae",
  phylum        = "Tracheophyta",
  class         = "Magnoliopsida",
  order         = "Rosales",
  family        = "Rosaceae",
  genus         = "Rubus",
  species       = NA,
  authority = NA)   # no species epithet --- hybrid
taxonomy_df = rbind(taxonomy_df, tayberry_row)
taxonomy_df = distinct(taxonomy_df)

# join back to original data
observedplants = observedplants %>%
  left_join(taxonomy_df, by = c("scientific_name" = "input_plants"))

observedplants$family[observedplants$scientific_name == "Hosta spp."] = "Asparagaceae"


# Check against PLANtITS database

getPhyloMatrix <- function(studyspecies){
  ##Prep phylogenetic tree to include species in 'studyspecies'
  
  #data from:
  ## Henriquez Piskulich, Patricia Andrea; Hugall,
  ## Andrew F.; Stuart-Fox, Devi (2023).  A supermatrix phylogeny of the
  ## world’s bees (Hymenoptera: Anthophila) [Dataset].
  ## Dryad. https://doi.org/10.5061/dryad.80gb5mkw1
  
  #download file BEE_mat7_fulltree_tplo35_sf20lp.nwk
  fulltree = ape::read.tree('data/BEE_mat7_fulltree.nwk')
  
  #make subset tree for study species
  studytips = fulltree$tip.label[grepl(paste(studyspecies, collapse = "|"), 
                                       fulltree$tip.label)]
  studytree = keep.tip(fulltree, studytips)
  
  #make a variance covariance matrix
  studycov = vcv.phylo(studytree)
  
  #change row and column names to "genus_species" format
  colnames(studycov) = sapply(strsplit(colnames(studycov), "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
  rownames(studycov) = sapply(strsplit(rownames(studycov), "_"), function(x){paste(x[[1]], x[[2]], sep="_")})
  
  return(studycov)
  
  #lots of repeat values--is that just because this phylogeny doesn't have a huge amount of resolution wrt bombus?
  
  
  
}



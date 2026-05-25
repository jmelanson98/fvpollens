library(biomformat)

biom = read_biom("3_data/merged_exported/feature-table.biom")
asv_table = as.data.frame(as.matrix(biom_data(biom)))

taxonomy = read.delim("3_data/merged_exported_taxonomy/taxonomy.tsv")
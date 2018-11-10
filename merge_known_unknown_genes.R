known_genes <- read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\training_all_genesinfo_corrected_classes_missing_values.csv", header=TRUE)
unknown_genes <- read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\merged_gene_protein_features_with_ids_missing_values.csv", header=TRUE)

data.merged <- rbind(known_genes, unknown_genes)

write.csv(data.merged,file="C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known_unknown_genes\\merged_known_unknown_genes.csv")

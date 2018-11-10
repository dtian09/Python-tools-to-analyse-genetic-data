data <- read.csv('/home/david/Dropbox/datasets/essential genes prediction/human essential genes/genes_of_different_essentiality_ensembleids.csv')
#check duplicate gene names
table(duplicated(data[,1]))
#check duplicate ensemble ids
table(duplicated(data[,2]))
#get the duplicate gene names
subset(data,duplicated(data[,1]),select = c(Associated.Gene.Name, Ensembl.Gene.ID))
#check unique gene names
genenames <- data[,1]#genenames is a vector, use genenames[1] to access 1st element of genename
genenames.unique <- unique(genenames)#unique works on a vector, dataframe or array (https://stat.ethz.ch/R-manual/R-devel/library/base/html/unique.html)





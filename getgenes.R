allgenes <- read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\allmousegenesfromMouseMine.csv", header=TRUE)
data <- read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\proteincodingmousegenesfromMouseMine.csv", header=TRUE)

#count no. of genes
#result: 22944 genes
nrow(data)

#check how many missing values
table(is.na(data))

#check duplicate records
table(duplicated(data))
#check duplicate MGI IDs
#result: 22944 unique MGI IDs, no duplicate MGI IDs
table(duplicated(data[,1]))
#check duplicate Gene Symbol
#result: 22942 unique gene symbols, 2 duplicate gene symbols
table(duplicated(data[,2]))
#get the 2 duplicate gene symbols
#result: Mar-01, Mar-02
subset(data,duplicated(data[,2]),select = c(Gene.ID,Gene.Symbol))

#get the original gene symbols and their duplicates
#result:
#      Gene.ID          Gene.Symbol    Gene.Name                                   Gene.Type            Gene.Organism
#5356  MGI:1913362      Mar-01       mitochondrial amidoxime reducing component 1  protein_coding_gene  Mus musculus
#6129  MGI:1914497      Mar-02       mitochondrial amidoxime reducing component 2  protein_coding_gene  Mus musculus
#8651  MGI:1920175      Mar-01       membrane-associated ring finger (C3HC4) 1     protein_coding_gene  Mus musculus
#10395 MGI:1925915      Mar-02       membrane-associated ring finger (C3HC4) 2     protein_coding_gene  Mus musculus
#
#
subset(data,data[,2]=="Mar-01" | data[,2] == "Mar-02")
 
#check duplicate Gene Name
#result: 22935 unique gene names, 9 duplicate gene names
table(duplicated(data[,3]))

#get the 9 duplicate gene names
#result: 
# Gene.ID                    Gene.Name
#MGI:1914819                     MGAT4 family, member C
#MGI:1915749        mitochondrial ribosomal protein L53
#MGI:1916719    transmembrane and coiled-coil domains 2
#MGI:1927343 ribosomal protein S6 kinase, polypeptide 2
#MGI:2135756                       heat shock protein 8
#MGI:2443419  ribosomal protein S6 kinase polypeptide 1
#MGI:3802945                sin3 associated polypeptide
#MGI:96243                       heat shock protein 2
#MGI:97606  prolactin family 3, subfamily d, member 1

subset(data,duplicated(data[,3]),select = c(Gene.ID,Gene.Name))
#get the original gene names and their duplicates
#result:
#          Gene.ID Gene.Symbol                                  Gene.Name           Gene.Type Gene.Organism
#    MGI:104558     Rps6ka1  ribosomal protein S6 kinase polypeptide 1 protein_coding_gene  Mus musculus
#    MGI:105384       Hspa8                       heat shock protein 8 protein_coding_gene  Mus musculus
#    MGI:1342290     Rps6ka2 ribosomal protein S6 kinase, polypeptide 2 protein_coding_gene  Mus musculus
#    MGI:1914805      Mgat4d                     MGAT4 family, member C protein_coding_gene  Mus musculus
#    MGI:1914819      Mgat4c                     MGAT4 family, member C protein_coding_gene  Mus musculus
#  MGI:1915090      Mrpl57        mitochondrial ribosomal protein L53 protein_coding_gene  Mus musculus
#  MGI:1915749      Mrpl53        mitochondrial ribosomal protein L53 protein_coding_gene  Mus musculus
#  MGI:1916125       Tmcc2    transmembrane and coiled-coil domains 2 protein_coding_gene  Mus musculus
#  MGI:1916503       Hspb2                       heat shock protein 2 protein_coding_gene  Mus musculus
#  MGI:1916719       Tmco2    transmembrane and coiled-coil domains 2 protein_coding_gene  Mus musculus
#  MGI:1927343     Rps6kb2 ribosomal protein S6 kinase, polypeptide 2 protein_coding_gene  Mus musculus
#  MGI:1929129       Sap30                sin3 associated polypeptide protein_coding_gene  Mus musculus
#  MGI:2135756       Hspb8                       heat shock protein 8 protein_coding_gene  Mus musculus
#  MGI:2443419     Rps6kc1  ribosomal protein S6 kinase polypeptide 1 protein_coding_gene  Mus musculus
#  MGI:2660935      Prl3d2  prolactin family 3, subfamily d, member 1 protein_coding_gene  Mus musculus
#  MGI:3802945       Sap25                sin3 associated polypeptide protein_coding_gene  Mus musculus
#  MGI:96243       Hspa2                       heat shock protein 2 protein_coding_gene  Mus musculus
#  MGI:97606      Prl3d1  prolactin family 3, subfamily d, member 1 protein_coding_gene  Mus musculus

subset(data,data[,3]=="MGAT4 family, member C" | data[,3]=="mitochondrial ribosomal protein L53" | data[,3]=="transmembrane and coiled-coil domains 2"| data[,3]=="ribosomal protein S6 kinase, polypeptide 2"| data[,3]=="heat shock protein 8" | data[,3]=="ribosomal protein S6 kinase polypeptide 1" | data[,3]=="sin3 associated polypeptide" | data[,3]=="heat shock protein 2" | data[,3] == "prolactin family 3, subfamily d, member 1")


#check duplicate Gene type (should be the same type)
table(duplicated(data[,4]))
#check duplicate gene organism (should be the same organism)
table(duplicated(data[,5]))

#get unique gene symbols
#data.unique <- unique(data[,2])
#length(data.unique) #22942 unique gene symbols, 2 duplicate gene symbols

#get unique gene names
#data.unique2 <- unique(data[,3])
#length(data.unique2) #22935 unique gene names, 9 duplicate gene names

#####get the unknown essentiality protein coding genes ####

lethalfinal <-read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\lethal_final.csv", header=TRUE)
viablefinal <-read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\viable_final.csv", header=TRUE)
#check duplicate MGI IDs in the lethal final genes
#result: 1308 MGI IDs, no duplicates
table(duplicated(lethalfinal[,2]))

#check duplicate MGI IDs in the viable final genes
#result: 3459 MGI IDs, no duplicates
table(duplicated(viablefinal[,2]))

#find positions of the lethal final genes in the data frame of all the protein coding genes
#result: all 1308 genes are found in the data frame
lethalfinalpos <- match(lethalfinal[,2],data[,1])
#check number of matched MGI IDs
length(lethalfinalpos)

#find positions of the viable final genes in the data frame of all the protein coding genes
#result: all 3459 genes are found in the data frame
viablefinalpos <- match(viablefinal[,2],data[,1])
#check number of matched MGI IDs
length(viablefinalpos)

#total number of known essentiality genes
#result: 4767
length(c(lethalfinalpos,viablefinalpos))

#get the positions of the genes with unknown essentiality in the data frame
allgenespos <- c(1:22944)
unknowngenespos <- allgenespos[-c(lethalfinalpos,viablefinalpos)]

#get the genes with unknown essentiality in the data frame

unknowngenes <- data[unknowngenespos,]

#total number of unknown essentiality genes
#result: 18179

nrow(unknowngenes)

#write the unknown essentiality genes to a csv file

write.csv(unknowngenes,file="C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\genesunknownessentiality.csv")

#write the unknown essentiality genes to a space delimited file
#write.table(unknowngenes,file="C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\genesunknownessentiality2.csv")

unknowngenes <- read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\genesunknownessentiality.csv", header=TRUE)

#check duplicate Gene Symbol
#result: 18178 unique gene symbols, 1 duplicate gene symbol
table(duplicated(unknowngenes[,2]))

#check duplicate Gene Name
#result: 18172 unique gene names, 7 duplicate gene names
table(duplicated(unknowngenes[,3]))



#create a training set by adding 153 new lethal genes and 153 new viable genes to balanced_missing_vals.csv (balanced known genes data set)
new_lethal <- read.csv('/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/93 features data/new_lethal_genes.csv', header=TRUE)
new_viable <- read.csv('/home/david/Dropbox/datasets/essential genes prediction/new viable genes/93 features data/new_viable_genes.csv', header=TRUE)
balanced <- read.csv('/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/balanced2_missing_vals.csv', header=TRUE)
lethal_indices = sample(1:nrow(new_lethal),153,replace=F)
viable_indices = sample(1:nrow(new_viable),153,replace=F)
new_lethal2 = new_lethal[lethal_indices,]
new_viable2 = new_viable[viable_indices,]
balanced2 <- rbind(balanced,new_lethal2,new_viable2)
#check duplicate genes
table(duplicated(balanced2))
#remove duplicate genes
balanced3 = unique(balanced2)
write.csv(balanced3,file="/home/david/Dropbox/datasets/essential genes prediction/train set/new_lethal_new_viable_balanced2_missing_vals.csv")
#create a test set
data1 <- read.csv("/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/103 features data/new_lethal_genes_not_in_train_set_missing_vals.csv", header=TRUE)
data2 <- read.csv("/home/david/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals.csv", header=TRUE)
data.merged <- rbind(data1,data2)
write.csv(data.merged,file="/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/new_lethal_new_viable_genes_not_in_train_set_missing_vals.csv")


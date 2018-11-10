#label genes or change labels of genes
#input: a csv data file containing genes
#output: a csv data file containing genes with class labels

data<-read.csv("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\GO\\data_GOannots.csv",header=TRUE)
#rows from 1 to 1308 are lethal genes: class 1
#rows from 1309 to end of data frame are viable genes: class 0

#get the column with name 'class'
#class <- data["class"]

for(i in 1:nrow(data)){
	data[i,"class"]<-'?'
}

for(i in 1:1308){
	data[i,"class"]<-1
}

for(i in 1309:nrow(data)){
	data[i,"class"]<-0
}

#write data to file skipping the MGI id field
write.csv(data[,-1], file = "C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\GO\\data_GOannots_labelled.csv", row.names=FALSE))

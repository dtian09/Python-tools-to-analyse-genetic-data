#input: a csv data file with class attribute the last column
#output: features of the data ranked in ascending order of p-value by Mann-Whitney U test 
#
#usage: R --slave --args '/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/new_lethal_new_viable_balanced_missing_vals.csv' < ./man_whitney_feature_rank.R > ranked_features.mann_whitney

args <- commandArgs(TRUE)

print(args[1])

datafile <- read.csv(args[1],header=TRUE)

class_index <- ncol(datafile)
last_feature_index <- ncol(datafile)-1

#get feature names
features <- names(datafile)

pvals <- c()
map <- new.env(hash=TRUE)#key: feature, value: p-value

#Mann-Whitney U test
for (i in 1:last_feature_index)
{
 result <- wilcox.test(features[i] ~ features[class_index], data=datafile)
 p <- result$p.value
 pvals[i] <- p
 map[[features[i]]] <- p
}

pvals <- sort(pvals)

#rank features by p-value
for(i in 1:last_feature_index)
{
  p <- pvals[i]
  for(j in 1:last_feature_index)
  {
   feature <- features[j]
   p2 <- map[[feature]]
   if(p == p2)
     print(paste0(paste0(feature,','),p))
  }
}

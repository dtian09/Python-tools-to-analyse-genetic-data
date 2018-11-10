#merge two files by rows into a new file
lethal_ids<-read.table("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\lethal_final_mgi",header=FALSE)
viable_ids<-read.table("C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\viable_final_mgi",header=FALSE)

lethal_viable_ids <- rbind(lethal_ids,viable_ids)

write.csv(lethal_viable_ids, file = "C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\lethal_viable_mgi_ids.csv")

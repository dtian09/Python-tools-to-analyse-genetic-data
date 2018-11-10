#program to collect Age feature of unknown genes
#
#input: Mouse_Genes_All_Chen_AB.csv
#	genenames_mgiids_uniprotids_ensembleids.csv
#output: age feature file
#
#format of Mouse_Genes_All_Chen_AB.csv
#
#Mouse_Gene,Age_Duplicate,Age_Singleton,Duplicate/Singleton,SCA,DCA,MRD
#ENSMUSG00000046709,104,NA, duplicate,,Opisthokonta,Eutheria
#ENSMUSG00000021936,400,NA, duplicate,,Opisthokonta,Euteleostomi
#ENSMUSG00000053137,535,NA, duplicate,,Opisthokonta,Vertebrata
#
def main():
	agefile = '/home/david/Dropbox/datasets/essential genes prediction/age/Mouse_Genes_All_Chen_AB.csv'
	idsfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/age/age.csv'
	ageL = [line.strip() for line in open(agefile)]
	idsL = [line.strip() for line in open(idsfile)]
	age_hash = {}
	ensembleids_of_idsfile = set()
	for line in ageL[1:len(ageL)]:
		vals = line.split(',')
		ensembleid = vals[0]
		age_hash[ensembleid]=line
	for line in idsL[1:len(idsL)]:
		vals = line.split(',')
		ensembleid = vals[3]
		ensembleids_of_idsfile.add(ensembleid)
	get_Age(ageL,ensembleids_of_idsfile,age_hash,outfile)

def get_Age(ageL,ensembleids_of_idsfile,age_hash,outfile):
	fw = open(outfile,'w')
	fw.write('Ensemble_ID,Age\n')
	for ensembleid in list(ensembleids_of_idsfile):
		if age_hash.get(ensembleid) != None :#get the age of the ensemble ids in the idsfile
			line = age_hash[ensembleid]
			vals = line.split(',')
			age_duplicate = vals[1]
			age_singleton = vals[2]
			if age_duplicate != 'NA': #get age duplicate
				age = age_duplicate
			elif age_singleton != 'NA':#get age singleton
				age = vals[2] #age singleton
			else:
				age = '?'
			fw.write(ensembleid+','+age+'\n')
		else:#age files does not have age of this ensemble id
			age = '?'
			fw.write(ensembleid+','+age+'\n')
	fw.close()

if __name__ == "__main__":
	main()


							
			
			
			




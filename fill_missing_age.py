#fill in missing age values of genes using ensemble ids
#
#input: newly_known_lethal_genes.csv
#	Mouse_Genes_All_Chen_AB.csv
#output: newly_known_lethal_genes.csv with missing age filled in

import sys

def main():
	#infile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known lethal genes\\newly_known_lethal_genes.csv'
	#infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\age\\Mouse_Genes_All_Chen_AB.csv'
	#outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known lethal genes\\newly_known_lethal_genes.csv'

	#infile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known viable genes\\newly_known_viable_genes.csv'
	#infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\age\\Mouse_Genes_All_Chen_AB.csv'
	#outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known viable genes\\newly_known_viable_genes.csv'
	infile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\merged_gene_protein_features_with_ids_missing_values.csv'
	infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\age\\Mouse_Genes_All_Chen_AB.csv'
	outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\merged_gene_protein_features_with_ids_missing_values.csv'
	
	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)

	hash_ages = {}
	for line in infile2L[1:len(infile2L)]:
		vals = line.split(',')
		ensemble_id = vals[0]
		if vals[1] != 'NA':
			age = vals[1]
		elif vals[2] != 'NA':
			age = vals[2]
		else:
			#age = '-1'
			age = '?'
		hash_ages[ensemble_id] = age

	#print(hash_ages)
	
	fw = open(outfile,'w')
	fw.write(infileL[0]+'\n')
	
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		age = vals[len(vals)-2]
		ensemble_id = vals[3]
		if ensemble_id == '?' or ensemble_id == '-1':
			fw.write(line+'\n')
		elif age == '-1' or age == '?':
			for val in vals[0:len(vals)-2]:
				fw.write(val+',')
			if hash_ages.get(ensemble_id) != None:
				fw.write(hash_ages[ensemble_id]+',')
			else:
				#fw.write('-1,')
				fw.write('?,')
			fw.write(vals[len(vals)-1]+'\n')
		else:
			fw.write(line+'\n')
	fw.close()

if __name__ == "__main__":
	main()				
			
		

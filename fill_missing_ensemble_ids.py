#fill in missing ensemble ids of newly known lethal genes dataset
#
#input: newly_known_lethal_genes.csv 
#	missing_ensemble_ids_of_new_lethal_genes.txt
#output: newly_known_lethal_genes.csv with ensemble ids filled in
import sys
import re

def main():
	#infile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known lethal genes\\newly_known_lethal_genes.csv'
	#infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known lethal genes\\missing_ensemble_ids_of_new_lethal_genes.txt'
	#outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known lethal genes\\newly_known_lethal_genes.csv'

	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/unknown genes16999/genenames16999_ensembleids.txt'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/unknown genes16999/12genenames_ensembleids_in_ensemble.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals2.csv'
	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	
	hash_ensemble_ids = {}
	for line in infile2L[1:len(infile2L)]:
		#m = re.match('^(\w+).+(ENSMUSG\d+)$',line)
		#Ensembl Gene ID,MGI symbol,MGI ID
		#ENSMUSG00000105096,Gbp10,MGI:4359647
		m = re.match('^(ENSMUSG\d+),([\w-]+\.*[\w-]+),MGI:\d+$',line)
		if m:
			#gene_name = m.group(1)
			gene_name = m.group(2)
			ensemble_id = m.group(1)			
			#ensemble_id = m.group(2)
			hash_ensemble_ids[gene_name] = ensemble_id
		else:
			print(line+" does not match pattern\n")
	#print(hash_ensemble_ids)
	
	#Fill in missing ensemble ids in data file
	fw=open(outfile,'w')
	fw.write(infileL[0]+"\n")
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		ensemble_id = vals[3]
		#if ensemble_id == '-1':
		if ensemble_id == '?':
			gene_name = vals[0]
			if hash_ensemble_ids.get(gene_name) != None:
				#GeneName,MGI_ID,UniProt_ID,Ensemble_ID
				fw.write(vals[0]+','+vals[1]+','+vals[2]+','+hash_ensemble_ids[gene_name])
				#,GeneLength,...,Class
				for val in vals[4:len(vals)]:
					fw.write(','+val)
				fw.write('\n')
			else:
				fw.write(line+'\n')
		else:
			fw.write(line+'\n')
	fw.close()

if __name__ == "__main__":
	main()
	 	


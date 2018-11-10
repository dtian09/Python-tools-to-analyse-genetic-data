#get new lethal genes from all the new lethal genes file from merged_gene_protein_features_with_ids.csv and 
#	training_all_genesinfo_corrected_classes_missing_values.csv
#input: gene_variants_with_phen_MP_0011100.csv (newly-known lethal genes containing duplicate gene names)
#	unknown genes dataset
#	known genes dataset
#output: a csv file containing the features of the newly known lethal genes	
#
import sys
import re

def main():
	infile = "C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\newly known lethal genes\\gene_variants_with_phen_MP_0011100.csv"
	infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\merged_gene_protein_features_with_ids_no_missing_vals.csv'
	infile3 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\training_all_genesinfo_corrected_classes_missing_values.csv'
	outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\newly known lethal genes\\newly_known_lethal_genes.csv'
	fw = open(outfile,"w")

	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)

	infile3L = [line.strip() for line in open(infile3)]

	if len(infile3L)==0:
		print(infile3+" is empty.")
		sys.exit(-1)

	lethal_gene_names = set()
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		lethal_gene_names.add(valsL[0])
			
	hash_unknown_genes_features = {} #key=gene name, value = features of gene
	for line in infile2L[1:len(infile2L)]:
		vals = line.split(',')
		gene_name = vals[0]
		line = ''
		for val in vals[0:len(vals)-1]:
			line += val+','
		line += 'Lethal'	
		hash_unknown_genes_features[gene_name] = line
		
	hash_known_genes_features = {} #key=gene name, value = features of gene
	for line in infile3L[1:len(infile3L)]:
		vals = line.split(',')
		gene_name = vals[0]
		line = ''
		for val in vals[0:len(vals)-1]:
			line += val+','
		line += 'Lethal'
		hash_known_genes_features[gene_name] = line
		
	#write lethal genes to file	
	fw.write(infile2L[0]+"\n")
	genes_in_known_genes_set = set()
	for gene_name in lethal_gene_names:
		if hash_known_genes_features.get(gene_name)!= None:
			genes_in_known_genes_set.add(gene_name)
		if hash_unknown_genes_features.get(gene_name)!= None:
			fw.write(hash_unknown_genes_features[gene_name]+"\n")
		elif hash_known_genes_features.get(gene_name)!= None:
			fw.write(hash_known_genes_features[gene_name]+"\n")
		else:
			print("gene name: "+gene_name+" is not in "+infile2+" and "+infile3+"\n")
	fw.close()
	print(str(len(genes_in_known_genes_set))+" lethal genes are in known genes data set : "+infile2)
		
if __name__ == "__main__":
	main()

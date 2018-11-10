#get the lethal genes which are not in the known gene set from all the newly-known lethal genes file.
#1. The duplicate gene names are removed from all the newly-known lethal genes file
#2. The new lethal genes which are not in the training set training_all_genesinfo_no_missing_vals.csv are retrieved from merged_gene_protein_features_with_ids.csv
#
#input: gene_variants_with_phen_MP_0011100.csv (newly-known lethal genes containing duplicate gene names)
#	unknown genes dataset
#	known genes dataset (training set)
#output: a csv file containing the features of the newly known lethal genes which are not in the training set	
#
import sys
import re

def main():

	infile = './new lethal genes/gene_variants_with_phen_MP_0011100.csv'
	infile2 = './unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	infile3 = './known essentiality genes/training_all_genesinfo.csv'#known gene set
	outfile = './new lethal genes/new_lethal_genes_not_in_train_set_ids.csv'

	#infile = "C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\newly known lethal genes\\gene_variants_with_phen_MP_0011100.csv"
	#infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\merged_gene_protein_features_with_ids_missing_vals.csv'
	#infile3 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\training_all_genesinfo.csv'
	#outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\newly known lethal genes\\new_lethal_genes_not_in_train_set_missing_vals.csv'
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
		gene_name = valsL[0]
		gene_name = gene_name.lower()
		lethal_gene_names.add(gene_name)
			
	hash_unknown_genes_features = {} #key=gene name, value = features of gene
	original_gene_names_of_unknown_genes = {}
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		original_gene_name = valsL[0]
		gene_name = original_gene_name.lower()
		original_gene_names_of_unknown_genes[gene_name] = original_gene_name
		line = ''
		for val in valsL[0:len(valsL)-1]:
			line += val+','
		hash_unknown_genes_features[gene_name] = line
		
	gene_names_of_known_genes = set() 
	for line in infile3L[1:len(infile3L)]:
		valsL = line.split(',')
		gene_name = valsL[0]
		gene_name = gene_name.lower()
		gene_names_of_known_genes.add(gene_name)
		
	#get the new lethal genes which are not in train set	
	#fw.write(infile2L[0]+"\n")
	fw.write('GeneName\n')
	genes_not_in_train_set = set()
	for gene_name in lethal_gene_names:
		if gene_name not in gene_names_of_known_genes:
				genes_not_in_train_set.add(gene_name)
				if hash_unknown_genes_features.get(gene_name) != None:
					#fw.write(hash_unknown_genes_features[gene_name]+"Lethal\n")				
					fw.write(original_gene_names_of_unknown_genes[gene_name]+'\n')
				else:
					print("gene name: "+gene_name+" is not in "+infile2+" and "+infile3+"\n")
	fw.close()
	print(str(len(lethal_gene_names)-len(genes_not_in_train_set))+" new lethal genes are in train set: "+infile3+"\n")	
	print(str(len(genes_not_in_train_set))+" new lethal genes are not in train set: "+infile3+"\n")
	
if __name__ == "__main__":
	main()

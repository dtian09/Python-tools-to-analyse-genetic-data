#Get the new viable genes which are not in the known gene set from all the newly-known unknown genes
#input: IMPC_all_genes_may2015.csv (all newly-known genes)
#	gene_variants_with_phen_MP_0011100.csv (newly-known lethal genes)
#	merged_gene_protein_features_with_ids.csv
#	training_all_genesinfo.csv
#output: a csv file containing the features of the newly known viable genes which are not in the training set	
#
#note: frrs1 and frrs1l are different genes. frrs1 and frrs1l are not in the known gene set and is in the unknown gene set. 

import sys
import re

def main():
	#infile = sys.argv[1]#./unknown essentiality genes/all newly known genes/IMPC_all_genes_may2015.csv
	#infile2 = sys.argv[2]#./unknown essentiality genes/newly known lethal genes/gene_variants_with_phen_MP_0011100.csv
	#outfile = sys.argv[3]]#./unknown essentiality genes/newly known viable genes/newly_known_viable_genes.csv
	infile = './all new known genes/IMPC_all_genes_may2015.csv'
	infile2 = './new lethal genes (lethal genes known since May 2015)/gene_variants_with_phen_MP_0011100.csv'
	infile3 = './unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	infile4 = './known essentiality genes/training_all_genesinfo.csv'#known gene set
	outfile = './new viable genes (viable genes known since May 2015)/new_viable_genes_not_in_train_set_genenames.csv'

	#infile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\all newly known genes\\IMPC_all_genes_may2015.csv'
	#infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\newly known lethal genes\\gene_variants_with_phen_MP_0011100.csv'
	#infile3 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\merged_gene_protein_features_with_ids_missing_vals.csv'
	#infile4 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\training_all_genesinfo.csv'
	#outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\newly known viable genes\\new_viable_genes_not_in_train_set_missing_vals.csv'
	
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
	
	infile4L = [line.strip() for line in open(infile4)]

	if len(infile4L)==0:
		print(infile4+" is empty.")
		sys.exit(-1)
	
	all_original_gene_names={}
	all_gene_names = set()
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_gene_name = valsL[0]
		gene_name = original_gene_name.lower()
		all_original_gene_names[gene_name] = original_gene_name
		all_gene_names.add(gene_name)
	
	lethal_gene_names = set()
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		original_gene_name = valsL[0]
		gene_name = original_gene_name.lower()
		lethal_gene_names.add(gene_name)
	
	viable_gene_names = all_gene_names.difference(lethal_gene_names)
	
	hash_unknown_genes_features = {} #key=gene name, value = features of gene
	for line in infile3L[1:len(infile3L)]:
		vals = line.split(',')
		gene_name = vals[0]
		gene_name = gene_name.lower()
		line = ''
		for val in vals[0:len(vals)-1]:
			line += val+','		
		hash_unknown_genes_features[gene_name] = line
		
	gene_names_of_known_genes = set() 
	for line in infile4L[1:len(infile4L)]:
		vals = line.split(',')
		gene_name = vals[0]
		gene_name = gene_name.lower()
		gene_names_of_known_genes.add(gene_name)
		
	#get the new viable genes which are not in the train set	
	#fw.write(infile3L[0]+"\n")
	fw.write('GeneName\n')
	genes_not_in_train_set = set()
	genes_not_in_known_and_unknown_genes_set = set()
	genes_not_in_known_genes_set_in_unknown_genes_set = set()
	for gene_name in viable_gene_names:
		if gene_name not in gene_names_of_known_genes:
			genes_not_in_train_set.add(gene_name)	
			if hash_unknown_genes_features.get(gene_name) != None:#if a gene is in the unknown gene set, take it
				#fw.write(hash_unknown_genes_features[gene_name]+"Viable\n")
				fw.write(all_original_gene_names[gene_name]+'\n')
				print(all_original_gene_names[gene_name])
				genes_not_in_known_genes_set_in_unknown_genes_set.add(all_original_gene_names[gene_name])			
			else:
				genes_not_in_known_and_unknown_genes_set.add(all_original_gene_names[gene_name])
				print("gene name: "+gene_name+" is not in "+infile3+" and "+infile4+"\n")
			
	fw.close()
	print(str(len(genes_not_in_known_genes_set_in_unknown_genes_set))+" new viable genes not in known genes set and in unknown genes set\n")
	print(str(len(genes_not_in_known_and_unknown_genes_set))+" new viable genes are not in known genes set and unknown genes set: "+infile4+" and "+infile3+"\n")
	print(str(len(genes_not_in_train_set))+" new viable genes are not in known genes set: "+infile4+"\n")
	print(str(len(viable_gene_names)-len(genes_not_in_train_set))+" new viable genes are in known genes set: "+infile4+"\n")
		
if __name__ == "__main__":
	main()


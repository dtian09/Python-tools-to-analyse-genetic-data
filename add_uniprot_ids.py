#Add Uniprot ids of genes to a csv file using their gene names
#
#input: new_viable_genes_not_in_train_set.csv 
#	new_viable_genes_uniprot_ids.csv
#output: new_viable_genes_not_in_train_set_missing_vals2.csv with a new column uniprot id 
import sys
import re

def main():
	
	#infile = 'c:/Users/David/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals.csv'
	#infile2 = 'c:/Users/David/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_uniprot_ids.csv'
	#outfile = 'c:/Users/David/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals2.csv'
	
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/103 features data/new_lethal_genes_not_in_train_set_missing_vals.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_uniprot_ids.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/103 features data/new_lethal_genes_not_in_train_set_missing_vals2.csv'
	
	infile = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_uniprot_ids.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals2.csv'
	

	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	
	hash_uniprot_ids = {}#key: gene name, value: uniprot id
	
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		gene_name = valsL[0]
		uniprot_id = valsL[1]
		if uniprot_id != '':
			hash_uniprot_ids[gene_name.lower()]=uniprot_id
		
	#Add uniprot ids to data file
	fw=open(outfile,'w')
	fieldsL = infileL[0].split(',')
	#write feature names
	#GeneName,MGI_ID,Uniprot_ID,Ensemble_ID,...
	fw.write(fieldsL[0]+','+fieldsL[1]+',UniProt_ID')
	for field in fieldsL[2:len(fieldsL)]:
		fw.write(','+field)
	fw.write('\n')
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_gene_name = valsL[0]
		gene_name = original_gene_name.lower()
		if hash_uniprot_ids.get(gene_name)!='None':
			uniprot_id = hash_uniprot_ids[gene_name]
		else:
			uniprot_id = '?'
			print(gene_name+' is not in '+infile2)
		line_new = valsL[0]#gene name feature
		line_new += ','+valsL[1]#mgi id feature
		line_new += ','+uniprot_id#Uniprot id feature
		for val in valsL[2:len(valsL)]:
			line_new += ','+val
		fw.write(line_new+'\n')
	fw.close()

if __name__ == "__main__":
	main()
	 	


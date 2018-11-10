#get the new viable genes from all the new known genes by subtracting the new lethal genes from all the new known genes
#input: IMPC_all_genes_may2015.csv (all newly-known unknown genes)
#	gene_variants_with_phen_MP_0011100.csv (newly-known lethal genes)
#	merged_gene_protein_features_with_ids.csv
#	training_all_genesinfo.csv
#output: a csv file containing gene names of the newly known viable genes
#
import sys
import re

def main():
	#infile = sys.argv[1]#./unknown essentiality genes/all newly known genes/IMPC_all_genes_may2015.csv
	#infile2 = sys.argv[2]#./unknown essentiality genes/newly known lethal genes/gene_variants_with_phen_MP_0011100.csv
	#outfile = sys.argv[3]]#./unknown essentiality genes/newly known viable genes/newly_known_viable_genes.csv
	#infile = './unknown essentiality genes/all newly known genes/IMPC_all_genes_may2015.csv'
	#infile2 = './unknown essentiality genes/newly known lethal genes/gene_variants_with_phen_MP_0011100.csv'
	#outfile = './unknown essentiality genes/newly known viable genes/newly_known_viable_genes.csv'

	infile = './all new known genes/IMPC_all_genes_may2015.csv'
	infile2 = './new lethal genes/gene_variants_with_phen_MP_0011100.csv'
	infile3 = './unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	infile4 = './known essentiality genes/training_all_genesinfo.csv'
	outfile = './new viable genes (viable genes known since May 2015)/new_known_viable_genes.csv'
	
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
			
	all_gene_names = set()
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		all_gene_names.add(valsL[0])
	
	lethal_gene_names = set()
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		lethal_gene_names.add(valsL[0])
	
	viable_gene_names = all_gene_names.difference(lethal_gene_names)
	
	hash_unknown_genes_features = {} #key=gene name, value = features of gene
	for line in infile3L[1:len(infile3L)]:
		vals = line.split(',')
		gene_name = vals[0]
		line = ''
		for val in vals[0:len(vals)-1]:
			line += val+','
		line += 'Viable'	
		hash_unknown_genes_features[gene_name] = line
	known_viable_genes = set()
	known_lethal_genes = set()	
	hash_known_genes_features = {} #key=gene name, value = features of gene
	for line in infile4L[1:len(infile4L)]:
		vals = line.split(',')
		gene_name = vals[0]
		line = ''
		for val in vals[0:len(vals)-1]:
			line += val+','
		line += 'Viable'
		if vals[len(vals)-1]=='Viable':
			known_viable_genes.add(gene_name)
		else:
			known_lethal_genes.add(gene_name)
		hash_known_genes_features[gene_name] = line
		
	#write viable genes to file	
	fw.write(infile3L[0]+"\n")
	genes_in_known_genes_set = set()
	genes_not_in_known_and_unknown_genes_set = set()
	k=0
	for gene_name in viable_gene_names:
		if gene_name in known_lethal_genes:
			print(gene_name+' is viable in new known gene set, but is lethal in know gene set')
			k += 1		
		if hash_known_genes_features.get(gene_name)!= None:
			genes_in_known_genes_set.add(gene_name)	
		if hash_unknown_genes_features.get(gene_name) != None:
			fw.write(hash_unknown_genes_features[gene_name]+"\n")
		elif hash_known_genes_features.get(gene_name)!= None:
			fw.write(hash_known_genes_features[gene_name]+"\n")
		else:
			genes_not_in_known_and_unknown_genes_set.add(gene_name)
			#print("gene name: "+gene_name+" is not in "+infile3+" and "+infile4+"\n")				
	fw.close()
	print('no. of viable genes in new known gene set which are labelled as lethal in known gene set: '+str(k))
	print("total no. of new viable genes: "+str(len(viable_gene_names))+"\n")
	print(str(len(genes_not_in_known_and_unknown_genes_set))+" viable genes are not in known and unknown genes data sets: "+infile3+" and "+infile4+"\n")
	print(str(len(genes_in_known_genes_set))+" viable genes are in known genes data set : "+infile4)

if __name__ == "__main__":
	main()


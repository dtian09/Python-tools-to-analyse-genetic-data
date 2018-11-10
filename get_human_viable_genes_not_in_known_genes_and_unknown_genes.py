#Get the human viable genes which have different gene names to the genes of the known essentiality mouse gene set and the unknown essentiality mouse gene set

import sys

def main():	
	infile='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'
	infile2 ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	infile3='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes.csv'
	infile4='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes2.csv'
		
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
	
	known_mouse_genes_and_unknown_mouse_genes = set()
	human_genes = set()
	original_gene_names = {}
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		known_mouse_genes_and_unknown_mouse_genes.add(genename)
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		known_mouse_genes_and_unknown_mouse_genes.add(genename)
	for line in infile3L:
		original_genename = line
		genename = original_genename.lower()
		human_genes.add(genename)
		original_gene_names[genename] = original_genename
	for line in infile4L:
		original_genename = line
		genename = original_genename.lower()
		human_genes.add(genename)
		original_gene_names[genename] = original_genename
	#Get the human viable genes which have different gene names to the genes of the known essentiality mouse gene set and the unknown essentiality mouse gene set
	human_genes_not_in_known_mouse_genes_and_unknown_mouse_genes = human_genes.difference(known_mouse_genes_and_unknown_mouse_genes)
	k = 0
	for gene in list(human_genes_not_in_known_mouse_genes_and_unknown_mouse_genes):
		print(original_gene_names[gene])
		k += 1
	print(str(k)+' human viable genes are not in known mouse gene set and unknown mouse gene set.')

if __name__ == "__main__":
        main()


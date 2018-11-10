#Compare the gene name or another gene id of 2 files

import sys

def main():
	#infile='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes2.csv'
	#infile='c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	#infile2 = 'c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_unigeneids.csv'
	#infile2 = 'c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/mgiids_unigeneids.csv'
	#infile2 = 'c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ensembleids_unigeneids.csv'
	#infile2 = 'c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/uniprotids_unigeneids.txt'
	infile = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_all_features_1031.csv'
	infile2='C:/Users/David/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15.csv'
	id_indx=0
	id_indx2=0
	infileL = [line.strip() for line in open(infile)]
	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	
	infile2L = [line.strip() for line in open(infile2)]
	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)	
	#mouse_genes_classes = {}
	genes_infile = set()
	duplicate_genes_infile = set()
	original_gene_names_infile = {}
	genes_infile2 = set()
	duplicate_genes_infile2 = set()
	original_gene_names_infile2 = {}
	(genes_infile,duplicate_genes_infile,original_gene_names_infile) = get_gene_names_from_csv_file_with_headers(infileL,id_indx,genes_infile,duplicate_genes_infile,original_gene_names_infile)
	(genes_infile2,duplicate_genes_infile2,original_gene_names_infile2) = get_gene_names_from_csv_file_with_headers(infile2L,id_indx2,genes_infile2,duplicate_genes_infile2,original_gene_names_infile2)
	print('infile contains '+str(len(genes_infile)))
	print('infile2 contains '+str(len(genes_infile2)))		
	if genes_infile == genes_infile2:
		print('infile and infile2 have the same genes')
	else:
		genes_of_infile_not_in_infile2 = genes_infile.difference(genes_infile2)
		genes_of_infile2_not_in_infile = genes_infile2.difference(genes_infile)
		common_genes_of_infile_and_infile2 = genes_infile.intersection(genes_infile2)
		print('infile contains '+str(len(genes_of_infile_not_in_infile2))+' genes which are not in infile2: '+str(genes_of_infile_not_in_infile2))
		print('infile2 contans '+str(len(genes_of_infile2_not_in_infile))+' genes which are not in infile: '+str(genes_of_infile2_not_in_infile))
		print('common genes of infile and infile2: '+str(common_genes_of_infile_and_infile2))
	print('duplicate genes of infile: '+str(duplicate_genes_infile))
	print('duplicate genes of infile2: '+str(duplicate_genes_infile2))
	'''
	print('unique gene names in infile: '+str(len(genes_infile)))
	print('unique gene names in infile2: '+str(len(genes_infile2)))
	print('duplicate gene names in infile: '+str(len(duplicate_genes_infile)))
	if len(duplicate_genes_infile) > 0:
		for gene in list(duplicate_genes_infile):
			print(original_gene_names_infile[gene])
	print('duplicate gene names in infile2: '+str(len(duplicate_genes_infile2)))
	#if len(duplicate_genes_infile2) > 0:
	#	for gene in list(duplicate_genes_infile2):
	#		print(original_gene_names_infile[gene])
	#print(str(len(genes_in_infile_not_in_infile2))+' genes of infile are not in infile2.')
	#for gene in list(genes_in_infile_not_in_infile2):
	#	print(original_gene_names_infile[gene])
	#print(str(len(genes_of_infile2_not_in_infile))+' genes of infile2 are not in infile.')
	#for gene in list(genes_of_infile2_not_in_infile):
	#	print(original_gene_names_infile[gene])
	if genes_infile2.issubset(genes_infile):
		print('infile contains all the genes of infile2')
	k = 0
	lethal_mouse_genes = set()
	print('common gene names of infile and infile2 '+str(len(common_genes_of_infile_and_infile2)))
	print('GeneName\tEssentiality of Human Gene\tEssentiality of Mouse Gene')
	for gene in list(common_genes_of_infile_and_infile2):
		if mouse_genes_classes[gene] == 'Lethal':
			k += 1
			lethal_mouse_genes.add(original_gene_names_infile[gene])
		print(original_gene_names_infile[gene]+'\tViable\t'+mouse_genes_classes[gene])
	print('no. of lethal mouse genes with the same gene names as the human viable genes: '+str(k))
	for gene in list(lethal_mouse_genes):
		print(gene)
	'''
	
def get_gene_names_from_csv_file_with_headers(infileL,id_indx,genes_infile,duplicate_genes_infile,original_gene_names_infile):
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_genename = valsL[id_indx]
		genename = original_genename.lower()
		if genename in genes_infile:
			duplicate_genes_infile.add(genename)
		genes_infile.add(genename)
		original_gene_names_infile[genename] = original_genename
	return (genes_infile,duplicate_genes_infile,original_gene_names_infile)

def each_line_is_a_gene_name_file(infileL,genes_infile,duplicate_genes_infile,original_gene_names_infile):
	for line in infileL:
		original_genename = line
		genename = original_genename.lower()
		if genename in genes_infile:
			duplicate_genes_infile.add(genename)
		else:
			genes_infile.add(genename)
		original_gene_names_infile[genename] = original_genename
	return (genes_infile,duplicate_genes_infile,original_gene_names_infile)
		
if __name__ == "__main__":
        main()


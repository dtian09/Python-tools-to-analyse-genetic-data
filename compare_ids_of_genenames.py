#compare the gene ids of two files

import sys

def main():
	id_indx = 0 #infile id index
	id_indx2 = 0 #infile2 id index

	#infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals_non-mouse_genes_removed_and_duplicates_removed.csv'
	
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/Match_LethalViableTest_MyVsDavid_Lethal.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/Lethal_blindTest_IDs.csv'
	#infile='/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/103 features data/new_lethal_genes_not_in_train_set_missing_vals2.csv'
	
	#infile2='C:/Users/David/Dropbox/datasets/essential genes prediction/test set/NewTestSet_Lethal_UniProtID&ProteinLen.csv'
	
	#infile='/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_uniprot_ids.csv'	
	#infile='C:/Users/David/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_uniprot_ids.csv'
	#infile='/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_uniprot_ids2.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/new viable genes/Match_LethalViableTest_MyVsDavid_Viable.csv'
	#infile='/home/david/Dropbox/datasets/essential genes prediction/new viable genes/gene features/gene_start_end_bps.csv'
	#infile= '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/gene features/transcript_length_exon_rank_exon_start_end.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/test set/NewTest_MismatchGeneProperties.csv'
	
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids2.csv'

	#infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_ensembleids.csv'
	#infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids.csv'
	#infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_all_uniprotids.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_all_features_1031.csv'
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_all_features_1031_2.csv'

	#infile='/home/david/Dropbox/datasets/essential genes prediction/test_ids'#example ids file 1
	#infile2='/home/david/Dropbox/datasets/essential genes prediction/test_ids2'#example ids file 2
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes2.csv'
	infileL = [line.strip() for line in open(infile)]
	infile2L = [line.strip() for line in open(infile2)]
	idSet=set()
	idSet2=set()
	gene_names_infile = {}
	gene_names_infile2 = {}
	for line in infileL[1:len(infileL)]:#first contains column names
		vals = line.split(',')
		gene_name = vals[0]
		original_gene_id = vals[id_indx]
		gene_id = original_gene_id.lower()
		idSet.add(gene_id)
		gene_name = gene_name.lower()
		gene_names_infile[gene_name] = gene_id
	for line in infile2L[1:len(infile2L)]:
	#for line in infile2L:#first line is beginning of data
		vals = line.split(',')
		gene_name = vals[0]
		original_gene_id = vals[id_indx2]
		gene_id = original_gene_id.lower()
		idSet2.add(gene_id)
		gene_name = gene_name.lower()
		gene_names_infile2[gene_name] = gene_id
	genenames = set(list(gene_names_infile.keys()))
	genenames2 = set(list(gene_names_infile2.keys()))
	print('no. of gene names in infile: '+str(len(genenames)))
	print('no. of gene names in infile2: '+str(len(genenames2)))
	common_genes_of_infile_and_infile2 = genenames.intersection(genenames2)
	print('common gene names of infile and infile2: '+str(len(common_genes_of_infile_and_infile2)))
	genes_of_infile_not_in_infile2 = genenames.difference(genenames2)
	genes_of_infile2_not_in_infile = genenames2.difference(genenames)
	print(str(len(genes_of_infile_not_in_infile2))+' genes of infile are not in infile2.')
	for gene in list(genes_of_infile_not_in_infile2):
		print(gene)
	print(str(len(genes_of_infile2_not_in_infile))+' genes of infile2 are not in infile.')
	for gene in list(genes_of_infile2_not_in_infile):
		print(gene)
	print('Compare ids of infile and infile2.')
	features1 = infileL[0]
	features1L = features1.split(',')
	features2 = infile2L[0]
	features2L = features2.split(',')
	print('id of infile: '+features1L[id_indx]+', id_indx: '+str(id_indx))
	print('id of infile2: '+features2L[id_indx2]+', id_indx2: '+str(id_indx2))
	if idSet == idSet2:
		print(infile+" and "+infile2+" have the same ids")
		sys.exit(0)
	d = 0#no. of gene names in infile with different ids to infile2
	d2 = 0#no. of gene names in infile2 with different ids to infile
	genenames = list(gene_names_infile.keys())
	id_comparisons = ''
	for genename in genenames:
		gene_id_infile = gene_names_infile[genename]
		if gene_names_infile2.get(genename) != None:
			gene_id_infile2 = gene_names_infile2[genename]
		else:
			gene_id_infile2 = 'missing_id_infile2'
		if gene_id_infile != gene_id_infile2:
			id_comparisons += genename+'\t'+gene_id_infile+'\t'+gene_id_infile2+'\n'
			d += 1
	print('Compare ids by the gene names of infile1')
	print('There are '+str(d)+' id mismatches.')	
	print('gene names\tid of infile\tid of infile2')
	print('---------------------------------------')
	print(id_comparisons)
	id_comparisons2=''
	genenames = list(gene_names_infile2.keys())
	for genename in genenames:
		gene_id_infile2 = gene_names_infile2[genename]
		if gene_names_infile.get(genename) != None:
			gene_id_infile = gene_names_infile[genename]
		else:
			gene_id_infile = 'missing_id_infile'
		if gene_id_infile != gene_id_infile2:
			id_comparisons2 += genename+'\t'+gene_id_infile+'\t'+gene_id_infile2+'\n'
			d2 += 1
	print('Compare ids by the gene names of infile2')	
	print('There are '+str(d2)+' id mismatches.')
	print('gene names\tid of infile\tid of infile2')
	print('---------------------------------------')
	print(id_comparisons2)

if __name__ == "__main__":
	main()

	

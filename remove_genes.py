#remove genes from a dataset
#input: gene data set in csv format
#output: gene data set with genes removed

import sys

def main():
	'''
	#remove ambiguous genes from the known genes set
	infile = '/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo_no_ambiguous_genes.csv'
	fw = open(outfile,'w')
	infileL = [line.strip() for line in open(infile)]
	ambiguous_viable_genes = ['Pcgf2', 'Ift88', 'Eya4', 'Abca1', 'Pthlh', 'Clcf1', 'Bag3', 'Satb1']
	ambiguous_lethal_genes = ['Syt1', 'Lama4', 'Pdpn', 'Hipk2','Afmid']
	ambiguous_viable_genes2 = set()
	ambiguous_lethal_genes2 = set()
	for gene in ambiguous_viable_genes:
		ambiguous_viable_genes2.add(gene.lower())
	for gene in ambiguous_lethal_genes:
		ambiguous_lethal_genes2.add(gene.lower())
	lethal_genes_removed = 0
	viable_genes_removed = 0
	fw.write(infileL[0]+'\n')
	for line in infileL[1:len(infileL)]:
		line2 = line.split(',')
		gene_name = line2[0]
		if gene_name.lower() in ambiguous_viable_genes2:
			viable_genes_removed +=1
		elif gene_name.lower() in ambiguous_lethal_genes2:
			lethal_genes_removed +=1
		else:
			fw.write(line+'\n')										
	fw.close()
	print(str(lethal_genes_removed)+' ambiguous lethal genes removed')
	print(str(viable_genes_removed)+' ambiguous viable genes removed')

	'''
	#remove the 29 non-mouse gene names and the 10 duplicate gene names from the unknown gene set
	#infile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/29_genenames_not_mouse_genes'
	infile3 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/10_duplicate_genenames'
	outfile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	fw=open(outfile,'w')
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
	non_mouse_genes = set()
	duplicate_genes = set()
	non_mouse_genes_removed = set()

	k = 0#no. of duplicate genes removed
	k2 = 0 #no. of non-mouse gene removed
	for gene in infile2L:
		non_mouse_genes.add(gene.lower())
	for gene in infile3L:
		duplicate_genes.add(gene.lower())
	genes_to_remove = non_mouse_genes.union(duplicate_genes)
	common_genes = non_mouse_genes.intersection(duplicate_genes)
	print('There are '+str(len(non_mouse_genes))+' non-mouse genes.')
	print('There are '+str(len(duplicate_genes))+' duplicate genes.')
	print('The non-mouse genes and duplicate genes have '+str(len(common_genes))+' common genes.')		
	fw.write(infileL[0]+'\n')
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		if genename in duplicate_genes:
			#print('duplicate gene: '+genename+' is removed')
			duplicate_genes.remove(genename)#remove the first occurrance of a duplicate gene
			k += 1		
		elif genename not in non_mouse_genes:#gene name is a mouse gene
			fw.write(line+'\n')
		else:
			#print('non-mouse gene: '+genename+' is removed')
			non_mouse_genes_removed.add(genename)
			k2 += 1							
	fw.close()
	print(str(len(genes_to_remove))+' genes to be removed in total.')
	print(str(k)+' duplicate genes have been removed.')
	print(str(k2)+' non-mouse genes have been removed.')
	
if __name__ == "__main__":
	main()
		

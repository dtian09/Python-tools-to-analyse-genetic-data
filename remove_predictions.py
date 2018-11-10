#remove the predictions of the non-mouse genes and the duplicate gene names from a prediction file output by a weka classifier
#input: a prediction file
#output: the prediction file with non-mouse genes and the duplicate gene names removed
#file format:
#
#
#=== Predictions on test data ===
#
#gene names  inst#     actual  predicted error prediction
#Plekhg2  1        1:?   2:Viable       0.522
#1700006E09Rik  2        1:?   2:Viable       0.835

#The non-mouse_gene name 'C2cd4cC2CD4 family' is not removed by this program because there is a space in the gene name. Remove this gene name manually.

import sys

def main():
	#remove the 29 non-mouse gene names and the 10 duplicate gene names from the unknown gene set
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_data_with_gene_names.random_forest_230_trees'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/29_genenames_not_mouse_genes'
	infile3 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/10_duplicate_genenames'
	outfile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_data_with_non-mouse_genes_removed_and_duplicates_removed.random_forest_230_trees'
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
		valsL = line.split('\t')
		original_genename = valsL[0]
		genename = original_genename.lower()
		if genename in duplicate_genes:
			print('duplicate gene: '+genename+' is removed')
			duplicate_genes.remove(genename)#remove the first occurrance of a duplicate gene
			k += 1
		elif genename not in non_mouse_genes:#gene name is a mouse gene
			fw.write(line+'\n')
		else:
			print('non-mouse gene: '+genename+' is removed')
			non_mouse_genes_removed.add(genename)
			k2 += 1							
	fw.close()
	print(str(len(genes_to_remove))+' genes to be removed in total.')
	print(str(k)+' duplicate genes have been removed.')
	print(str(k2)+' non-mouse genes have been removed.')
	if len(non_mouse_genes) > k2:
		print('non mouse genes not removed: '+str(non_mouse_genes.difference(non_mouse_genes_removed)))
	
if __name__ == "__main__":
	main()


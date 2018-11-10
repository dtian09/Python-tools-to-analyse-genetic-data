#check whether human viable genes are in known genes data set
#output: a table as follows
#
#Human Gene Name | in known genes set | Essentiality | Known Essentiality(Y/N)
#-----------------------------------------------------------------------------
#gene1		 | Yes		      | Essential    | Y
#gene2           | No                 | unknown      | N
#gene3		 | Yes		      | non-essential| Y

import sys
import re

def main():
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo_no_ambiguous_classes.csv'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/new_lethal_new_viable_genes_not_in_train_set_missing_vals.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes.csv'
	
	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	
	lethal=set()
	viable=set()
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')		
		gene = vals[0]
		essentiality = vals[len(vals)-1]
		if essentiality == 'Lethal':
			lethal.add(gene)
		else:
			viable.add(gene)
	human_viable_genes=set()
	k = 0 #count no. of human genes in known gene set
	print('Human Gene Names  | in known genes set | Essentiality | Known Essentiality(Yes/No)|')
	print('-----------------------------------------------------------------------------------')
	for gene in infile2L:
		if gene in human_viable_genes:
			print('duplicate gene: '+gene+' in human viable genes set')
		human_viable_genes.add(gene)
		if gene in lethal:
			l = gene+' | Yes      | lethal  |  Yes |'
			k += 1
		elif gene in viable:
			l = gene+' | Yes      | viabal  |  Yes |'
			k += 1
		else:
			l = gene+' | No       | unknown |  No  |'
	
	percentage = float(k)/float(len(human_viable_genes))*100
	print(str(percentage)+'% of human viable genes are in known genes set')
	print('no. of lethal genes in known genes set: '+str(len(lethal)))
	print('no. of viable genes in known genes set: '+str(len(viable)))
	print('no. of human viable genes: '+str(len(human_viable_genes)))

if __name__ == "__main__":
	main()


	

#get human viable genes from supplementary table 3 in MacArthur paper, Supplementary Tables 5 and 6 in the Sulem et al paper 
#Supplemental table 3 in the MacArthur et al paper contains viable genes.
#With Supplementary table 5, only consider genes with 0 in the het and fail categories as true viable genes. 
#For Supplementary table 6, consider only genes with 0 in the mismatch and fail categories as true viable genes.

import sys
import re

def main():
	table3='c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/MacArthur_supplementary table 3.txt'
	table5='c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/Sulem_supplementary table 5.txt'
	table6='c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/Sulem_supplementary table 6.txt'
	outfile='c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes.csv'
	fw=open(outfile,'w')

	table3L = [line.strip() for line in open(table3)]

	if len(table3L)==0:
		print(table3+" is empty.")
		sys.exit(-1)

	table5L = [line.strip() for line in open(table5)]

	if len(table5L)==0:
		print(table5+" is empty.")
		sys.exit(-1)
	
	table6L = [line.strip() for line in open(table6)]

	if len(table6L)==0:
		print(table6+" is empty.")
		sys.exit(-1)
	#get genes from table 3
	genes = set()
	for line in table3L:
		vals = line.split(' ')
		gene = vals[0]
		genes.add(gene)
	print('no. of human viable genes in table 3: '+str(len(genes)))
	genesL = list(genes)
	for gene in genesL:
		fw.write(gene+'\n')
	#get genes from table 5
	#fw.write('\n===table 5===\n')
	genes2 = set()
	for line in table5L:
		vals = line.split(' ')
		het = vals[len(vals)-2]
		fail = vals[len(vals)-1]
		if het == '0' and fail == '0':
			gene = vals[0]
			if gene not in genes:
				genes2.add(gene)
	print('no. of human viable genes in table 5 not in table 3: '+str(len(genes2)))
	genes2L = list(genes2)
	for gene in genes2L:		
		fw.write(gene+'\n')
	#fw.write('\n===table 6===\n')
	genes3 = set()
	for line in table6L:
		vals = line.split(' ')
		mismatch = vals[len(vals)-2]
		fail = vals[len(vals)-1]
		if mismatch == '0' and fail == '0':
			gene = vals[0]
			if gene not in genes2 and gene not in genes:
				genes3.add(gene)
	print('no. of human viable genes in table 6 not in tables 3 and 5: '+str(len(genes3)))
	genes3L = list(genes3)
	for gene in genes3L:
		fw.write(gene+'\n')
	fw.close()
	t = len(genes) + len(genes2) + len(genes3)
	print('total no. of human viable genes: '+str(t))

if __name__ == "__main__":
	main()
	 	


	



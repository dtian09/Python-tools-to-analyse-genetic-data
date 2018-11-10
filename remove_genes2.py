#remove genes from a file

import sys

def main():
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals5.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals6.csv'
	infileL = [line.strip() for line in open(infile)]
	fw = open(outfile,'w')
	
	genes_to_removeL = ['Gm10427','Gm10047','Gm10811','4930423O20Rik','Gm9778','F530104D19Rik','Spi15','Spi16','Neurod5','Gm1971','Gm10484']
	genes_to_removeSet = set()
	removed_genes = set()
	for gene in genes_to_removeL:
		genes_to_removeSet.add(gene.lower())

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	fw.write(infileL[0]+'\n')
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		gene_name = valsL[0]
		gene_name = gene_name.lower()
		if gene_name not in genes_to_removeSet:
			fw.write(line+'\n')
		else:
			removed_genes.add(gene_name)
	fw.close()
	if genes_to_removeSet != removed_genes:
		print("genes_to_removeSet is not equal to removed_genes")
	else:
		print("genes_to_removeSet is equal to removed_genes")
	print(str(len(genes_to_removeL))+' genes should be removed')
	print(str(len(removed_genes))+' genes were removed')

if __name__ == "__main__":
	main()

	


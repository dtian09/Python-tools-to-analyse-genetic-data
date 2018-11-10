#collect GO terms from MGI database
#
#input: a file containing the GO annotation of genes from MGI database
#	essentiality of the genes (1, 0 or ?)
#output: a data file containing the GO annotations of genes in csv format
#
#GO terms file format:
#
#MGI Gene/Marker ID,GO ID,Term
#MGI:101758,GO:0008150,biological_process
#MGI:101758,GO:0005575,cellular_component
#MGI:101758,GO:0003674,molecular_function
#MGI:101760,GO:0006397,mRNA processing
#MGI:101760,GO:0048025,"negative regulation of mRNA splicing, via spliceosome"
#MGI:101760,GO:0006355,"regulation of transcription, DNA-templated"
#MGI:101760,GO:0006396,RNA processing
#MGI:101760,GO:0008380,RNA splicing
#
#This program uses collectfeatures_GO.py
#
#author: David Tian	tiand@cs.man.ac.uk

import sys
import re
import collectfeatures_GO

def main():
	unknown_GOannot_file = sys.argv[1]
	essentiality = sys.argv[2]
	outfile = sys.argv[3]
	GOterms = set()
	
	unknowngenesL = [line.strip() for line in open(unknown_GOannot_file)]

	if len(unknowngenesL)==0:
		print(unknown_GOannot_file+" is empty.")
		sys.exit(-1)
		
	unknown_hash = {}#unknown_hash: key=ensemble gene id, value = a set of GO terms
	unknown_genesAll = collectfeatures_GO.get_genes(unknowngenesL[1:len(unknowngenesL)])#return set of all genes
	unknowngenesL = collectfeatures_GO.remove_lines_with_no_GOannotations(unknowngenesL)#return list of genes with GO annotations
	GOterms = collectfeatures_GO.get_GOterms2(unknowngenesL[1:len(unknowngenesL)],GOterms)
	GOtermsL = list(GOterms)
	unknown_hash = collectfeatures_GO.get_GOterms_Only(unknowngenesL[1:len(unknowngenesL)],unknown_hash)
	#####write GO terms to an output file
	fw = open(outfile,"w")
	fw.write("MGI_Gene_ID,")
	for GOterm in GOtermsL:
		fw.write(GOterm+",")
	fw.write("class\n")
	collectfeatures_GO.write_GOannotations2(unknown_hash,GOtermsL,fw,essentiality)
	fw.close()
	ks = list(unknown_hash.keys())
	genes_with_no_GOannots = len(unknown_genesAll)-len(ks)
	print("no. of genes in the file "+unknown_GOannot_file+": "+str(len(unknown_genesAll)))
	print("no. of genes with GO annotations in the file "+unknown_GOannot_file+": "+str(len(ks)))
	print("no. of genes with no GO annotations in the file "+unknown_GOannot_file+": "+str(genes_with_no_GOannots))
	print("no. of GO annotations of genes in the file "+unknown_GOannot_file+": "+str(len(GOtermsL)))

	outfileL = [line.strip() for line in open(outfile)]
	print("no. of genes in file "+outfile+": "+str(len(outfileL)-1))
	
if __name__ == "__main__":
	main()

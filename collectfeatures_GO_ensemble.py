#Collect Go annotation, GO domains of genes from Ensembles database
# 
#input: a file containing the GO annotation of genes from ensemble database
#output: a data file containing the GO annotations of genes in csv format
#
#GO terms file format:
#
#MGI ID,GO Term Name,GO domain
#MGI:101787,cellular_component,cellular_component
#MGI:101787,molecular_function,molecular_function
#MGI:101787,oxidoreductase activity,molecular_function
#MGI:101787,intracellular,cellular_component
#MGI:101787,cell,cellular_component
#
#This program uses collectfeatures_GO.py
#
#author: David Tian	tiand@cs.man.ac.uk

import sys
import re
import collectfeatures_GO

def main():
	unknown_GOannot_file = sys.argv[1]
	outfile = sys.argv[2]
	GOterms = set()
	
	#create unlablled data file
	unknowngenesL = [line.strip() for line in open(unknown_GOannot_file)]

	if len(unknowngenesL)==0:
		print(unknown_GOannot_file+" is empty.")
		sys.exit(-1)
		
	unknown_hash = {}#unknown_hash: key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	unknown_genesAll = collectfeatures_GO.get_genes(unknowngenesL[1:len(unknowngenesL)])#return set of all genes
	unknowngenesL = collectfeatures_GO.remove_lines_with_no_GOannotations(unknowngenesL)#return list of genes with GO annotations
	GOterms = collectfeatures_GO.get_GOterms(unknowngenesL[1:len(unknowngenesL)],GOterms)
	GOtermsL = list(GOterms)
	unknown_hash = collectfeatures_GO.get_GOterms_GOdomains(unknowngenesL[1:len(unknowngenesL)],unknown_hash)
	(bio,cell,mole) = GOtermsDomains(unknowngenesL[1:len(unknowngenesL)])
	#printGOtermsOfGODomains(bio,cell,mole)
	#####write GO terms to an output file
	fw = open(outfile,"w")
	fw.write("MGI_Gene_ID,")
	for GOterm in GOtermsL:
		fw.write(GOterm+",")
	fw.write("class\n")
	collectfeatures_GO.write_GOannotations(unknown_hash,GOtermsL,fw,"?")
	fw.close()
	ks = list(unknown_hash.keys())
	genes_with_no_GOannots = len(unknown_genesAll)-len(ks)
	print("no. of genes in the file "+unknown_GOannot_file+": "+str(len(unknown_genesAll)))
	print("no. of genes with GO annotations in the file "+unknown_GOannot_file+": "+str(len(ks)))
	print("no. of genes with no GO annotations in the file "+unknown_GOannot_file+": "+str(genes_with_no_GOannots))
	print("no. of GO annotations of genes in the file "+unknown_GOannot_file+": "+str(len(GOtermsL)))
	print("#######################################################################################")
	print("no. of GO annotations of biological process in the file "+unknown_GOannot_file+": "+str(len(bio)))
	print("no. of GO annotations of cellular component in the file "+unknown_GOannot_file+": "+str(len(cell)))
	print("no. of GO annotations of molecular function in the file "+unknown_GOannot_file+": "+str(len(mole)))
	
	outfileL = [line.strip() for line in open(outfile)]
	print("no. of genes in file "+outfile+": "+str(len(outfileL)-1))
	
def printGOtermsOfGODomains(bio,cell,mole):
	print("GO terms of biological process: ")
	print(bio)
	print("#######################################################################################")
	print("GO terms of cellular component: ")
	print(cell)
	print("#######################################################################################")
	print("GO terms of molecular function: ")
	print(mole)
	
def GOtermsDomains(GOannot_fileL):
	#input: a list with each element a line in a GO annotation file
	#output: the set of GO annotations of biological process, 
	#	 the set of GO annotations of cellular component, 
	#	 the set of GO annotations of molecular function
	
	biological = set()
	cellular = set()
	molecular = set()
	for line in GOannot_fileL:
		m = re.match("[\w:]+,(.+),(cellular_component|biological_process|molecular_function)",line)
		if m != None:
			GOterm = m.group(1)
			GOterm2 = re.sub(',','_',GOterm)
			GOterm2 = re.sub(' ','_',GOterm2)
			GOterm2 = re.sub('"','',GOterm2)
			GOdomain = m.group(2)
			if GOterm2 != '':
				if GOdomain == 'cellular_component':
					cellular.add(GOterm2)
				elif GOdomain == 'biological_process':
					biological.add(GOterm2)
				else:
					molecular.add(GOterm2)
	return (biological,cellular,molecular)
	
if __name__ == "__main__":
	main()
		

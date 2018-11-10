#Collect Go annotations, GO domains of genes
#input: a file containing the GO annotations of lethal genes
#	a file containing the GO annotations of viable genes
#	a file containing the GO annotations of unknown essentiality genes (unknown_essentiality_GO_terms2.csv)
#output: a labelled data file containing the common GO annotations of lethal, viable and unknown essentiality genes
#	 a unlabelled data file containing the common GO annotations of lethal, viable and unknown essentiality genes
#
#ensemble GO terms file format (using ensemble ids):
#
#Ensembl Gene ID,GO Term Name,GO domain
#ENSMUSG00000064341,mitochondrion,cellular_component
#ENSMUSG00000064341,mitochondrial respiratory chain complex I,cellular_component
#ENSMUSG00000064341,NADH dehydrogenase (ubiquinone) activity,molecular_function
#
#ensemble GO terms file format (using MGI ids):
#
#MGI ID,GO Term Name,GO domain
#MGI:101787,cellular_component,cellular_component
#MGI:101787,molecular_function,molecular_function
#MGI:101787,oxidoreductase activity,molecular_function
#MGI:101787,intracellular,cellular_component
#MGI:101787,cell,cellular_component

#author: David Tian	tiand@cs.man.ac.uk

import sys
import re

def main():
	lethal_GOannot_file = sys.argv[1]
	viable_GOannot_file = sys.argv[2]
	unknown_GOannot_file = sys.argv[1]
	labelled_data = sys.argv[4]
	unlabelled_data = sys.argv[2]
	GOterms = set()
	
	lethal_hash = {}#lethal_hash: key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	viable_hash = {}#viable_hash: key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	
	lethal = [line.strip() for line in open(lethal_GOannot_file)]
	viable = [line.strip() for line in open(viable_GOannot_file)]
	
	if len(lethal)==0:
		print(lethal_GOannot_file+" is empty.")
		sys.exit(-1)
	if len(viable)==0:
		print(viable_GOannot_file+" is empty.")
		sys.exit(-1)
		
	lethal_genes = get_genes(lethal[1:len(lethal)])#return set of all lethal genes
	viable_genes = get_genes(viable[1:len(viable)])#return set of all viable genes
	
	lethal = remove_lines_with_no_GOannotations(lethal)
	viable = remove_lines_with_no_GOannotations(viable)
	
	GOterms = get_GOterms(lethal[1:len(lethal)],GOterms)
	GOterms = get_GOterms(viable[1:len(viable)],GOterms)
	GOtermsL = list(GOterms)
	print("total no. of GO terms of all genes: "+str(len(GOterms)))

	lethal_hash = get_GOterms_GOdomains(lethal[1:len(lethal)],lethal_hash)#skip the 0th element, the header line in list lethal
	viable_hash = get_GOterms_GOdomains(viable[1:len(viable)],viable_hash)
	
	#create unlablled data file
	unknowngenes = [line.strip() for line in open(unknown_GOannot_file)]
	unknown_hash = {}#viable_hash: key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	
	unknowngenes = remove_lines_with_no_GOannotations(unknowngenes)
	unknown_hash = get_GOterms_GOdomains(unknowngenes[1:len(unknowngenes)],unknown_hash)
	frequentGOtermsL = get_frequent_GOterms(unknown_hash,GOtermsL)#get GO terms of lethal and viable genes with frequencies > 0 
	print("total no. of frequent GO annotations: "+str(len(frequentGOtermsL)))
	fw2 = open(unlabelled_data,"w")
	fw2.write("Ensemble_Gene_ID,")
	for GOterm in frequentGOtermsL:
		fw2.write(GOterm+",")
	fw2.write("class\n")
	write_GOannotations(unknown_hash,frequentGOtermsL,fw2,"?")
	fw2.close()

	#create labelled data file
	fw = open(labelled_data,"w")
	fw.write("Ensemble_Gene_ID,")
	for GOterm in frequentGOtermsL:
		fw.write(GOterm+",")
	fw.write("class\n")
	write_GOannotations(lethal_hash,frequentGOtermsL,fw,"1")
	write_GOannotations(viable_hash,frequentGOtermsL,fw,"0")
	fw.close()
	ks = list(lethal_hash.keys())	
	ks2 = list(viable_hash.keys())
	lethal_no_GOannot = len(lethal_genes) - len(ks)
	viable_no_GOannot = len(viable_genes) - len(ks2)
	print("total no. of lethal genes with no GO annotations: "+str(lethal_no_GOannot))
	print("total no. of viable genes with no GO annotations: "+str(viable_no_GOannot))
	print("total no. of genes with no GO annotations: "+str(lethal_no_GOannot+viable_no_GOannot))

def get_frequent_GOterms(GOannot_hash,GOtermsL):
	#get GO terms of the genes with frequencies > 0 	
	frequentGOterms = set()
	ensemble_ids = list(GOannot_hash.keys())
	for ensemble_id in ensemble_ids:
		(GOterms,GOdomains) = GOannot_hash.get(ensemble_id)
		for GOterm in GOterms:
			if GOterm in GOtermsL:
				frequentGOterms.add(GOterm)
	return list(frequentGOterms)
		
def get_genes(GOannot_fileL):
    #get ensemble ids of all genes
	genes = set()
	for line in GOannot_fileL:
		m = re.match("([\w:]+),.+",line)
		if m != None:
			genes.add(m.group(1))
	return genes
	
def remove_lines_with_no_GOannotations(GOannot_fileL):
	#remove genes with no GO annotations and lines with no GO annotations
	i=0
	for line in GOannot_fileL:
		m = re.match("[\w:]+,\s*,.*",line)
		if m != None:	
			del GOannot_fileL[i]
		i+=1
	return GOannot_fileL
	
def get_GOterms(GOannot_fileL,GOterms):
	#input: a list with each element a line in a GO annotation file
	#	a empty set of GO annotations
	#output: a set of GO annotations in the GO annotation file
	#GO terms file format:
	#Ensembl Gene ID,GO Term Name,GO domain
	#ENSMUSG00000064341,mitochondrion,cellular_component
	#ENSMUSG00000064341,mitochondrial respiratory chain complex I,cellular_component
	#ENSMUSG00000064341,NADH dehydrogenase (ubiquinone) activity,molecular_function
	
	for line in GOannot_fileL:
		m = re.match("([\w:]+),(.+),.*",line)
		if m != None:
			GOterm = m.group(2)
			GOterm2 = re.sub(',','_',GOterm)
			GOterm2 = re.sub(' ','_',GOterm2)
			GOterm2 = re.sub('"','',GOterm2)
			if GOterm2 != '':
				GOterms.add(GOterm2)
	return GOterms

def get_GOterms_GOdomains(GOannot_fileL,GOannot_hash):
	#input: a list with each element a line in a GO annotation file
	#	a empty hashtable
	#output: a hashtable: key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	#GO terms file format:
	#Ensembl Gene ID,GO Term Name,GO domain
	#ENSMUSG00000064341,mitochondrion,cellular_component
	#ENSMUSG00000064341,mitochondrial respiratory chain complex I,cellular_component
	#ENSMUSG00000064341,NADH dehydrogenase (ubiquinone) activity,molecular_function
	for line in GOannot_fileL:
		#m = re.match("([\w]+),(.+),(cellular_component|biological_process|molecular_function)",line)
		m = re.match("([\w:]+),(.+),(.*)",line)
		if m != None:#add the gene to hash table if it has GO annotations; if gene has no GO annotation, then do not add the gene to hash table
			ensemble_id = m.group(1)
			GOterm = m.group(2)
			GOdomain = m.group(3)
			GOterm2 = re.sub(',','_',GOterm)
			GOterm2 = re.sub(' ','_',GOterm2)
			GOterm2 = re.sub('"','',GOterm2)
			if GOannot_hash.get(ensemble_id)==None:
				GOterms = set()
				GOdomains = set()
				if GOterm2 != '':
					GOterms.add(GOterm2)
					GOdomains.add(GOdomain)
					GOannot_hash[ensemble_id]=(GOterms,GOdomains)
			else:
				v = GOannot_hash[ensemble_id]
				GOterms = v[0]
				GOdomains = v[1]
				if GOterm2 != '':
					GOterms.add(GOterm2)
					GOdomains.add(GOdomain)
					GOannot_hash[ensemble_id]=(GOterms,GOdomains)
	return GOannot_hash
			
def write_GOannotations(GOannot_hash,GOtermsL,fw,essentiality):
	ensemble_ids = list(GOannot_hash.keys())
	for ensemble_id in ensemble_ids:
		fw.write(ensemble_id+',')
		v = GOannot_hash[ensemble_id]
		GOtermsGene = v[0]
		for GOterm in GOtermsL:
			if GOterm in GOtermsGene:
				fw.write('1,')
			else:
				fw.write('0,')
		fw.write(essentiality+"\n")

if __name__ == "__main__":
	main()
		

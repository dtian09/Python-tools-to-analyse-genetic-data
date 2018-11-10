#Get the frequency of GO annotations of genes
#
#input: a csv file containing the GO annotation of genes (output by Ensemble database)
#
#output: GO annotations ranked by their frequencies
#
#GO annotation file format:
#
#MGI ID,GO Term Name,GO domain
#MGI:64341,mitochondrion,cellular_component
#MGI:64341,mitochondrial respiratory chain complex I,cellular_component
#MGI:64341,NADH dehydrogenase (ubiquinone) activity,molecular_function
#
#This program uses collectfeatures_GO.py
#
import sys
import collectfeatures_GO
import re

def main():
	
	#GOannot_lethal = sys.argv[1]
	#GOannot_viable = sys.argv[2]
	#GOannot_unknown = sys.argv[3]
	
	GOannot_unknown = sys.argv[1]

	GOannot_unknownL = [line.strip() for line in open(GOannot_unknown)]
	
	if len(GOannot_unknown)==0:
		print(GOannot_unknown+" is empty.")
		sys.exit(-1)

	#1. find most frequent GO terms of unknown genes
	#freq_hash: hashtable(key=a GO annotation, value = frequency of the GO annotation)
	#freq_hash2: hashtable (key=frequency of GO annotation, value = set of GO annotations with the frequency)
	#####
	print("#### Find most frequent GO annots of genes ####")
	genesL = collectfeatures_GO.remove_lines_with_no_GOannotations(GOannot_unknownL)#return list of genes with GO annotations
	freq_hash = create_GOannot_hash(genesL[1:len(genesL)])
	freq_hash2 = create_GOannot_hash2(freq_hash)
	freq_L = list(freq_hash2.keys())		
	freq_L.sort()#sort no. of GO annotations in ascending order
	k=len(freq_L)
	mostFreqGOannots(freq_L,freq_hash2,k)	
	#2. count frequency of GO annots of genes
	print("\n#### count frequency of GO annots of genes ####")
	count_frequency_of_GOterms_of_genes(GOannot_unknownL)
	
def count_frequency_of_GOterms_of_genes(genesGOannotL):

	GOannot_hash = {} #key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	GOannot_hash2 = {} #key=no. of GO annotation of a gene, value = no. of genes having that no. of GO terms
	
	genesL = collectfeatures_GO.remove_lines_with_no_GOannotations(genesGOannotL)#return list of genes with GO annotations
	GOannot_hash = collectfeatures_GO.get_GOterms_GOdomains(genesL[1:len(genesL)],GOannot_hash)
	
	GOannot_hash2 = get_freq(GOannot_hash)
	
	no_of_GOannot_L = list(GOannot_hash2.keys())		
	no_of_GOannot_L.sort()#sort no. of GO annotations in ascending order
	printGOannots(no_of_GOannot_L,GOannot_hash2)
	ks = list(GOannot_hash.keys())
	printGOannots2(no_of_GOannot_L,GOannot_hash2,len(ks))
	
def printGOannots(no_of_GOannot_L,GOannot_hash):
	print("no. of GO terms in a gene,no. of genes")
	for no_of_GOannot in no_of_GOannot_L:
		no_of_genes = GOannot_hash[no_of_GOannot]
		print(str(no_of_GOannot)+','+str(no_of_genes))			

def printGOannots2(no_of_GOannot_L,GOannot_hash,total_no_of_genes):
	print("no. of GO terms in a gene,percentage of genes")
	for no_of_GOannot in no_of_GOannot_L:
		no_of_genes = GOannot_hash[no_of_GOannot]
		print(str(no_of_GOannot)+','+str(no_of_genes/total_no_of_genes*100))			

def get_freq(GOannot_hash):
	#input: hashtable (key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains)
	#output: hashtable (key= no. of GO terms, value = no. of genes having that no. of GO terms)
	freq_hash ={}#key= no. of GO terms of a gene, value = no. of genes having that no. of GO terms 
	ensemble_ids = list(GOannot_hash.keys())
	for ensemble_id in ensemble_ids:
		v=GOannot_hash[ensemble_id]
		no_of_GOterms = len(v[0])
		if freq_hash.get(no_of_GOterms) == None:
			freq_hash[no_of_GOterms]=1
		else:
			freq_hash[no_of_GOterms]+=1	
	return freq_hash			

def create_GOannot_hash(GOannot_fileL):
	#input: a list with each element a line in a GO annotation file
	#output: hashtable(key=a GO annotation, value = frequency of the GO annotation)
	freq_hash={}
	for line in GOannot_fileL:
		m = re.match("^[\w:]+,(.+),.+$",line)
		if m != None:
			GOterm = m.group(1)
			GOterm2 = re.sub(',','_',GOterm)
			GOterm2 = re.sub(' ','_',GOterm2)
			GOterm2 = re.sub('"','',GOterm2)
			if GOterm2 != '' and freq_hash.get(GOterm2) == None:
				freq_hash[GOterm2]=1
			elif GOterm2 != '' and freq_hash.get(GOterm2) != None:
				freq_hash[GOterm2] += 1			
	return freq_hash			
				
def create_GOannot_hash2(GOannot_hash):
	#input: hashtable (key=a GO annotation, value = frequency of the GO annotation)
	#output: hashtable (key=frequency of GO annotation, value = set of GO annotations with the frequency)
	GOannots = list(GOannot_hash.keys())
	freq_hash = {}
	for GOannot in GOannots:
		freq = GOannot_hash[GOannot]
		if freq_hash.get(freq)==None:
			s = set()
			s.add(GOannot)
			freq_hash[freq] = s
		else:
			s = freq_hash[freq]
			s.add(GOannot)
			freq_hash[freq]=s
	return freq_hash

def mostFreqGOannots(sorted_L,GOannot_hash,k):
	#input: a list of frequencies in ascending order,
	#	hashtable (key=frequency, value = the set of the GO terms with the frequency)
	#	k (k most frequent GO terms)
	#
	#Print top k most frequent GO terms
	#frequency: no. of GO terms with the frequency
	i = 0
	
	sorted_L.reverse()
	for freq in sorted_L:
		if i<k:
			GOannots = GOannot_hash[freq]
			if len(GOannots) <= 10:#print the GO annots if there are <= 10 GO annots
				print(str(freq)+':'+str(len(GOannots))+" "+str(GOannots))
			else:
				print(str(freq)+':'+str(len(GOannots)))
				
		i += 1
	
if __name__ == "__main__":
	main()

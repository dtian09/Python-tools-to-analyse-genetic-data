#Provide functions to collect Go terms, GO domains of genes from Ensemble and MGI databases
#
#Ensemble GO terms file format:
#
#MGI ID,GO Term Name,GO domain
#MGI:101787,cellular_component,cellular_component
#MGI:101787,molecular_function,molecular_function
#MGI:101787,oxidoreductase activity,molecular_function
#MGI:101787,intracellular,cellular_component
#MGI:101787,cell,cellular_component
#
#MGI GO terms file format:
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
#author: David Tian	tiand@cs.man.ac.uk

import sys
import re

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
	#get gene ids
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
	#get GO terms from ensemble GO annot file
	#input: a list with each element a line in a GO annot file
	#	a empty set of GO annotations
	#output: a set of GO annotations in the GO annot file
		
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
	
def get_GOterms2(GOannot_fileL,GOterms):
	#get GO terms from MGI GO annot file
	#input: a list with each element a line in a GO annot file
	#	a empty set of GO annotations
	#output: a set of GO annotations in the GO annot file
		
	for line in GOannot_fileL:
		m = re.match("([\w:]+),.+,(.+)",line)
		if m:
			GOterm = m.group(2)
			if GOterm != 'biological_process' and GOterm != 'molecular_function' and GOterm != 'cellular_component':
				GOterm2 = re.sub(',','_',GOterm)
				GOterm2 = re.sub(' ','_',GOterm2)
				GOterm2 = re.sub('"','',GOterm2)
				if GOterm2 != '':
					GOterms.add(GOterm2)
	return GOterms

def get_GOterms_Only(GOannot_fileL,GOannot_hash):
	#input: a list with each element a line in a GO annotation file
	#	a empty hashtable
	#output: a hashtable: key=MGI gene id, value = a set of GO terms of the gene
	
	for line in GOannot_fileL:
		m = re.match("([\w:]+),.+,(.+)",line)
		if m != None:#add the gene to hash table if it has GO annotations; if gene has no GO annotation, then do not add the gene to hash table
			mgi_id = m.group(1)
			GOterm = m.group(2)
			GOterm2 = re.sub(',','_',GOterm)
			GOterm2 = re.sub(' ','_',GOterm2)
			GOterm2 = re.sub('"','',GOterm2)
			if GOannot_hash.get(mgi_id) == None:
				GOterms = set()
				if GOterm2 != '':
					GOterms.add(GOterm2)
					GOannot_hash[mgi_id] = GOterms
			else:
				GOterms = GOannot_hash[mgi_id]
				if GOterm2 != '':
					GOterms.add(GOterm2)
					GOannot_hash[mgi_id] = GOterms
	return GOannot_hash
	
def get_GOterms_GOdomains(GOannot_fileL,GOannot_hash):
	#input: a list with each element a line in a GO annotation file
	#	a empty hashtable
	#output: a hashtable: key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	
	for line in GOannot_fileL:
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
	#input: a hashtable: key=ensemble gene id, value = a tuple of set of GO terms and set of GO domains
	#output: a data file
	gene_ids = list(GOannot_hash.keys())
	for gene_id in gene_ids:
		fw.write(gene_id+',')
		v = GOannot_hash[gene_id]
		GOtermsGene = v[0]
		for GOterm in GOtermsL:
			if GOterm in GOtermsGene:
				fw.write('1,')
			else:
				fw.write('0,')
		fw.write(essentiality+"\n")

def write_GOannotations2(GOannot_hash,GOtermsL,fw,essentiality):
	#input: a hashtable: key=mgi gene id, value = a set of GO terms
	#output: a data file
	gene_ids = list(GOannot_hash.keys())
	for gene_id in gene_ids:
		fw.write(gene_id+',')
		GOtermsGene = GOannot_hash[gene_id]
		for GOterm in GOtermsL:
			if GOterm in GOtermsGene:
				fw.write('1,')
			else:
				fw.write('0,')
		fw.write(essentiality+"\n")

def write_GOannotations3(hash1,hash2,GOtermsAllL,fw,essentiality):
	#input: a hashtable: key=mgi id, value = a set of GO terms
	#	a hashtable: key=mgi id, value = a set of GO terms
	#	list of all GO terms
	#	output file
	#	essentiality '1', '0' or '?'
	#output: a merge data file
	#
	#Merge the GO annots of the genes contained in hash1 and hash2 into a data file so that the data file contains the genes which have GO annots stored in hash1 and GO annots stored in hash2
	#
	#example merged data file:
	#MGI id,T1,T2,T3,T4,...,TK,class
	#MGI:1,1,0,1,0,...,1,?
	#MGI:2,0,0,1,1,...,1,?
	#MGI:3,1,0,0,1,...,0,?
	#MGI:4,1,0,1,0,...,1,?
	#MGI:5,1,1,1,1,...,0,?
	#MGI:6,0,0,0,0,...,0,?
	#MGI:7,1,0,1,0,...,1,?
	#MGI:8,0,1,1,1,...,1,?
	#MGI:9,1,0,0,1,...,0,?
	#
	#where genes MGI:1,MGI:2,MGI:3 are contained in hash1 only;
	#      genes MGI:4,MGI:5,MGI:6,MGI:7 are contained in hash2 only;
	#      genes MGI:8 and MGI:9 are contained in both hash1 and hash2.
	ks1 = hash1.keys()
	gene_ids = list(ks1)
	gene_idsset = set(gene_ids)
	ks2 = hash2.keys()
	gene_ids2set = set(list(ks2))
	gene_ids2diff = gene_ids2set.difference(gene_idsset)#mgi ids of hash2 which are not contained in hash1
	gene_ids2diffL = list(gene_ids2diff)
	#write GO annots of all the genes contained in hash1. Some of these genes are contained in hash1 only. Some are contained in both hash1 and hash2
	for gene_id in gene_ids:
		fw.write(gene_id+',')
		GOterms1 = hash1[gene_id]
		if hash2.get(gene_id) == None:
			GOterms3 = GOterms1
		else:
			GOterms2 = hash2[gene_id]
			GOterms3 = GOterms1.union(GOterms2)
			
		for GOterm in GOtermsAllL:
			if GOterm in GOterms3:
				fw.write('1,')
			else:
				fw.write('0,')
		fw.write(essentiality+"\n")#unknown class
	#write the annots of those genes contained in hash2 only but are not contained in hash1  
	for gene_id in gene_ids2diffL:
		fw.write(gene_id+',')
		GOterms = hash2[gene_id]		
		for GOterm in GOtermsAllL:
			if GOterm in GOterms:
				fw.write('1,')
			else:
				fw.write('0,')
		fw.write(essentiality+"\n")#unknown class
	

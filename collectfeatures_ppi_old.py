#program to get known ppi(protein-protein interactions) and known and predicted ppi of longest proteins encoded by unknown mouse genes
#
#input: a PPI file downloaded from I2D database (e.g. i2d.2_3.Public.MOUSE.tab containing known and predicted PPIs of all mouse genes from the I2D database)
#	uniprot-mapped_unknown_essentiality_genes.pepstats (a csv file containing the Uniprot ids of the longest proteins encoded by unknown mouse genes)
#output: a known PPI file (file containing the known PPIs retrieved from known PPI databases)
#	 a predicted and known PPI file (file containing the known PPIs and predicted PPIs from known and predicted PPI databases)
#
#databases containing known PPIs of mouse genes
#
#BioGrid_Mouse 	  	
#BIND_Mouse 	  	
#Chen_PiwiScreen 	  
#DIP_Mouse 	  	
#I2D-c_Fiona_MOUSE 	
#I2D-c_Mouse 	  	
#IntAct_Mouse 	  	
#INNATEDB_Mouse 	
#KIM_MYC 	  	
#MGI 	  	
#MINT_Mouse 	
#WangEScmplx 	
#WangEScmplxlow 
#WangEScoIP 	
#
#other databases contain the predicted PPIs of mouse genes
#
#format of the PPI file downloaded from l2D database
#Dataset	SwissProt1	SwissProt2
#BIND_Mouse	Q5U413	A0A5E3
#BioGrid	Q6PDY0	A0AUP1
#MINT	Q6PDY0	A0AUP1
#IntAct	Q9D6I9	A0AUP1
#IntAct	Q9D845	A0AUP1
#IntAct_Mouse	Q9QY53	A0AUV1
#IntAct_Mouse	Q9QY53	A0JLV3
#IntAct_Mouse	Q8BP00	A0JNT0
#
#format of output PPI file
#Uniprot1	Uniprot2
#Q9QY53	A0AUV1
#Q9QY53	A0JLV3
#Q8BP00	A0JNT0

import sys
import re

def main():
	infile = sys.argv[1]#a ppi file downloaded from l2D database
	infile2 = sys.argv[2]#uniprot-mapped_unknown_essentiality_genes.pepstats
	outfile = sys.argv[3]#a known ppi file
	outfile2 = sys.argv[4]#a known and predicted ppi file

	ppiL = [line.strip() for line in open(infile)]

	if len(ppiL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	pepstatsL = [line.strip() for line in open(infile2)]

	if len(pepstatsL)==0:
		print(infile2+" is empty.")
		sys.exit(-1)

	known_databases = set(['BioGrid_Mouse','BIND_Mouse','Chen_PiwiScreen','DIP_Mouse','I2D-c_Fiona_Mouse,I2D-c_Mouse','IntAct_Mouse','INNATEDB_Mouse','KIM_MYC','MGI','MINT_Mouse','WangEScmplx','WangEScmplxlow','WangEScoIP'])

	(all_known_ppiL,all_known_predicted_ppiL) = get_all_known_ppi_and_known_and_predicted_ppi(ppiL,known_databases)
	#remove any duplicate ppis and self-loops
	all_known_ppiL = remove_duplicate_selfloops(all_known_ppiL)
	all_known_predicted_ppiL = remove_duplicate_selfloops(all_known_predicted_ppiL)

	#get the known PPIs and known and predicted PPIs between longest proteins of unknown genes only
	(known_ppiL,known_predicted_ppiL) = get_longest_protein_ppi(all_known_ppiL,all_known_predicted_ppiL,pepstatsL[1:len(pepstatsL)])
	
	print("no. of known ppi: "+str(len(known_ppiL)-1))
	print("no. of known and predicted ppi: "+str(len(known_predicted_ppiL)-1))

	fw = open(outfile,'w')
	fw.write('uniprot1\tuniprot2\n')
	for ppi in known_ppiL:
		fw.write(ppi+'\n')
	fw.close()

	fw2 = open(outfile2,'w')
	fw2.write('uniprot1\tuniprot2\n')
	for ppi in known_predicted_ppiL:
		fw2.write(ppi+'\n')
	fw2.close()

def get_all_known_ppi_and_known_and_predicted_ppi(ppiL,known_databases):
	all_known_ppiL = []
	all_known_predicted_ppiL = []
	lowercase_names = set()
	#change database names to lowercases
	for database in known_databases:
		lowercase_name = database.lower()
		lowercase_names.add(lowercase_name)
	for ppi in ppiL:
		l = ppi.split('\t')
		if l[0].lower() in lowercase_names:
			all_known_ppiL.append(ppi)
		else:
			all_known_predicted_ppiL.append(ppi)
	return (all_known_ppiL,all_known_predicted_ppiL)

def remove_duplicate_selfloops(ppiL):
	collected_ppi = set()
	ppi_no_duplicate_selfloop = []

	for ppi in ppiL:
		ppi_L = ppi.split('\t')
		a_ppi = ppi_L[1]+'\t'+ppi_L[2]
		if a_ppi not in collected_ppi and ppi_L[1] != ppi_L[2]:#check whether the ppi has been collected and is a self interaction e.g. Q6PDY0 Q6PDY0
			ppi_no_duplicate_selfloop.append(ppi_L)
			collected_ppi.add(a_ppi)
	return ppi_no_duplicate_selfloop

def get_longest_protein_ppi(all_known_ppiL,all_known_predicted_ppiL,pepstatsL):
	#get the ppi of the longest protein encoded by unknown genes
	#input: all_known_ppiL (list of all known ppi)
	#	all_known_and_predicted_ppiL (list of all known ppi and known and predicted ppi)
	#output: a list of known ppi of longest proteins of unknown genes
	#	 a list of known ppi and known and predicted ppi of longest proteins of unknown genes

	uniprot_ids = set()
	known_longest_protein_ppiSet = set()
	known_and_predicted_longest_protein_ppiSet = set()

	#get longest proteins encoded by unknown genes
	for line in pepstatsL:
		l = line.split(',')
		uniprot_id = l[1]#uniprot id of a longest protein encoded by an unknown gene
		uniprot_ids.add(uniprot_id)
		
	#get known and predicted ppi of longest proteins encoded by unknown genes from known interaction databases
	for l in all_known_ppiL:
		if l[1] in uniprot_ids:
			if l[2] in uniprot_ids:
				known_longest_protein_ppiSet.add(l[1]+'\t'+l[2])
	for l in all_known_predicted_ppiL:
		if l[1] in uniprot_ids:
			if l[2] in uniprot_ids:
				known_and_predicted_longest_protein_ppiSet.add(l[1]+'\t'+l[2])
	return (list(known_longest_protein_ppiSet),list(known_and_predicted_longest_protein_ppiSet))
	
if __name__ == "__main__":
	main()
	

#program to get known ppi(protein-protein interactions) of longest reviewed proteins encoded by mouse genes
#
#input: i2d.2_3.Public.MOUSE.tab (PPI file containing known and predicted PPIs of all mouse genes from the I2D database)
#	a csv file containing the Uniprot ids of the longest reviewed proteins of the mouse genes retrieved from Uniprot (e.g. new_lethal_new_viable_ids.csv or uniprot-mapped_unknown_essentiality_genes.pepstats)
#	a fasta file containing all the proteins encoded by each mouse gene (the fasta file is obtained by mapping MGI ids or gene names of the genes to Uniprot ids using Uniprot and downloading the fasta file from Uniprot)
#	a pepstats file containing the properties of all the proteins encoded by each mouse gene (the pepstats file is output by Pepstats tool taking the above fasta file as input)  
#	
#output: a known PPI file (file containing the known PPIs retrieved from known PPI databases)
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
#format of output PPI file
#Uniprot1	Uniprot2
#Q9QY53	A0AUV1
#Q9QY53	A0JLV3
#Q8BP00	A0JNT0

import sys
import re
import genes

def main():
	ppi_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/PPI network features/i2d.2_3.Public.MOUSE.tab'
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_ids.csv'
	fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_all_proteins_of_each_gene.fasta'
	pepstats_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_all_proteins_of_each_gene.pepstats'
	known_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_known_ppi.tab'
	predicted_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_predicted_ppi.tab'
	ppiL = [line.strip() for line in open(ppi_file)]
	genesidL = [line.strip() for line in open(genesid_file)]
	fastaL = [line.strip() for line in open(fasta_file)]
	pepstatsL = [line.strip() for line in open(pepstats_file)]
	known_ppi_databases = set(['BioGrid_Mouse','BIND_Mouse','Chen_PiwiScreen','DIP_Mouse','I2D-c_Fiona_Mouse,I2D-c_Mouse','IntAct_Mouse','INNATEDB_Mouse','KIM_MYC','MGI','MINT_Mouse','WangEScmplx','WangEScmplxlow','WangEScoIP'])
	longest_proteins_of_all_genes = get_longest_proteins_of_all_genes(genesidL,fastaL,pepstatsL)	
	(known_ppis,predicted_ppis) = retrieve_known_and_predicted_ppi(ppiL,genesidL,longest_proteins_of_all_genes,known_ppi_databases)
	known_ppiL = remove_duplicate_self_ppi(list(known_ppis))#remove any duplicate ppis and self ppi of known ppis
	predicted_ppiL = remove_duplicate_self_ppi(list(predicted_ppis))#remove any duplicate ppis and self ppi of predicted ppis
	fw = open(known_ppi_outfile,'w')
	fw.write('uniprot1\tuniprot2\n')
	for ppi in known_ppiL:
		fw.write(ppi[0]+'\t'+ppi[1]+'\n')
	fw.close()
	fw2 = open(predicted_ppi_outfile,'w')
	fw2.write('uniprot1\tuniprot2\n')
	for ppi in predicted_ppiL:
		fw2.write(ppi[0]+'\t'+ppi[1]+'\n')
	fw2.close()

def get_longest_proteins_of_all_genes(genesidL,fastaL,pepstatsL):
	genenames = set()
	longest_proteins_of_all_genes = set()
	for line in genesidL[1:len(genesidL)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		genenames.add(genename)	
	(proteins,proteins_of_genes_not_to_collect,proteins_with_no_gene_names) = genes.get_proteins(list(genenames),fastaL)#proteins is a hashtable: key=gene_name, value=set of tuples (protein name, UniProt id, UniProt entry name)
	genesL = list(proteins.keys())
	for gene in genesL:
		longest_proteinsSet = get_longest_proteins(gene,proteins,pepstatsL) #longest_proteinsSet = a set of uniprot ids of the longest proteins of a gene
		longest_proteins_of_all_genes = longest_proteins_of_all_genes.union(longest_proteinsSet)
	return 	longest_proteins_of_all_genes
	
def retrieve_known_and_predicted_ppi(ppiL,genesidL,longest_proteins_of_all_genes,known_ppi_databases):
	known_ppis = set()#set of known ppis between longest reviewed proteins and longest proteins of other genes
			  #{(p1,p2),(p2,p3),...,(pl,pk)}
	predicted_ppis = set()#predicted ppis between longest reviewed proteins and longest proteins of other genes
			  #{(p1,p2),(p2,p3),...,(pl,pk)}
	(all_known_ppi,all_predicted_ppi) = get_all_known_ppi_and_known_and_predicted_ppi(ppiL,list(known_ppi_databases))
	retrieved_longest_reviewed_proteins = get_retrieved_longest_reviewed_proteins(genesidL)
	#loop:
	#	Take a PPI (pi,pj) from ppiL (file containing all the PPIs) 
	#	If pi is a retrieved longest reviewed protein and pj is a longest protein of another gene
	#          or
	#   	   pj is a retrieved longest reviewed protein and pi is a longest protein of another gene
	#	Then If (pi,pj) is a known PPI
	#     	     Then keep (pi,pj) as a known PPI
	#     	     Else keep (pi,pj) as a predicted PPI	
	#until ppiL is an empty list
	for ppi in ppiL[1:length(ppiL)]:
		l = ppi.split('\t')
		pi = l[0]
		pj = l[1]
		if pi in retrieved_longest_reviewed_proteins and pj in longest_proteins_of_all_genes:
			if (pi,pj) in all_known_ppi or (pj,pi) in all_known_ppi:
				known_ppis.add((pi,pj))
			else:
				predicted_ppis.add((pi,pj))
		elif pj in retrieved_longest_reviewed_proteins and pi in longest_proteins_of_all_genes:
			if (pj,pi) in all_known_ppi or (pi,pj) in all_known_ppi:
				known_ppis.add((pj,pi))
			else:
				predicted_ppis.add((pj,pi))
	return (known_ppis,predicted_ppis)			
						
def get_retrieved_longest_reviewed_proteins(genesidL):
	#get the uniprot ids of the retrieved longest reviewed proteins of the genes in the genesidL
	retrieved_longest_reviewed_proteins = set()
	for genesidL[1:len(genesidL)]:
		vals = line.split(',')	
		genename = vals[0]
		uniprotid = vals[2]
		retrieved_longest_reviewed_proteins.add(uniprotid)
	return retrieved_longest_reviewed_proteins

def get_all_known_ppi_and_predicted_ppi(ppiL,known_databases):
	all_known_ppis = set()
	all_predicted_ppis = set()
	lowercase_names = set()
	#change database names to lowercases
	for database in known_databases:
		lowercase_name = database.lower()
		lowercase_names.add(lowercase_name)
	for ppi in ppiL:
		l = ppi.split('\t')
		database = l[0]
		p1 = l[1]
		p2 = l[2]
		if database.lower() in lowercase_names:
			all_known_ppis.add((p1,p2))
		else:
			all_predicted_ppis.add((p1,p2))
	return (all_known_ppis,all_predicted_ppis)

def remove_duplicate_self_ppi(ppiL):
	#remove duplicate ppis and self ppis from ppiL
	collected_ppi = set()
	no_duplicate_self_ppi = set()
	for ppi in ppiL:
		p1 = ppi[0]
		p2 = ppi[1]
		if (p1,p2) not in collected_ppi and (p2,p1) not in collected_ppi and p1 != p2:#if a ppi has not been collected and is not a self interaction e.g. Q6PDY0 Q6PDY0
			no_duplicate_self_ppi.add((p1,p2))
			collected_ppi.add((p1,p2))
	return list(no_duplicate_self_ppi)

def get_longest_proteins(gene,proteins,pepstatsL):
	#return: a list of uniprot ids of the longest proteins of a gene
	max_l=0
	proteins_length = {} #key = uniprot_entry_name of a protein, value = length of the protein
	longest_proteinsSet = set()
	proteins_of_geneSet = proteins[gene]
	proteins_of_geneL = list(proteins_of_geneSet)
	#compute the length of each protein of the gene; store the lengths of proteins in the hashtable proteins_length and find the longest length of these proteins
	for protein_of_gene in proteins_of_geneL:
		uniprotid = protein_of_gene[1]
		entry_name = protein_of_gene[2]	
		#print('max_l: '+str(max_l))			
		(proteins_length,l) = get_length_of_protein(proteins_length,entry_name,pepstatsL)
		#print('l: '+str(l))
		if l > max_l:
			max_l = l
	#find the uniprot ids of the longest protein(s) of the gene
	for protein_of_gene in proteins_of_geneL:
		uniprotid = protein_of_gene[1]		
		entry_name = protein_of_gene[2]
		length = proteins_length[entry_name]
		if length == max_l:
			longest_proteinsSet.add(uniprotid)
	return longest_proteinsSet

def get_length_of_protein(proteins_length,uniprot_entry_name,pepstatsL):
	#input a dictionary: key = uniprot id of a protein, value = protein length
	#      uniprot_entry_name	    
	#      a list with each element a line in a pepstats output file
	#output: length of protein sequence
	if proteins_length.get(uniprot_entry_name)==None:
		length = get_protein_length(uniprot_entry_name,pepstatsL)
		proteins_length[uniprot_entry_name] = length
		return (proteins_length,length)
	else:
		length = proteins_length[uniprot_entry_name]
		return (proteins_length,length)
	
def get_protein_length(entry_name,pepstatsL):
	#input: uniprot entry name of a protein sequence e.g. A2AHY8_MOUSE
	#       a list with each element a line in a pepstats output file
	#output: length (residues) of the protein	
	i=0
	residues=''
	protein_found=False
	while(i<len(pepstatsL)):
		line = pepstatsL[i]	
		matchObj = re.match("PEPSTATS of "+entry_name+".+\\s+to\\s+(\d+)",line)	#match line: "PEPSTATS of A2AHY8_MOUSE from 1 to 747"
		if matchObj:
			residues = matchObj.group(1)				
			protein_found=True
			i+=1
			break		
		else:
			i+=1
	if protein_found==False:
		print ('Uniprot entry name '+entry_name+' is not found in pepstats output file')
	return int(residues)	

if __name__ == "__main__":
	main()
	

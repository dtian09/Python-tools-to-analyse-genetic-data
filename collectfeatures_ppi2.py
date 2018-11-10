#program to get the known ppi and the predicted ppi between each longest reviewed protein and the longest protein among all the interactions of the protein
#e.g. p1 is a longest reviewed protein retrieved and has interactions with p2, p3, p4 and p5. P5 is the longest among p2, p3, p4 and p5. The interaction between p1 and p5 is chosen
#input: i2d.2_3.Public.MOUSE.tab
#	gene ids file e.g. new_lethal_new_viable_ids.csv
#	a review status file obtained by mapping MGI ids of all protein-coding mouse genes (proteincodingmousegenesfromMouseMine.csv) to uniprot ids on uniprot website
#	a review status file obtained by mapping gene names of all protein-coding mouse genes (proteincodingmousegenesfromMouseMine.csv) to uniprot ids on uniprot website
#	a review status file containing the Uniprot ids which are not contained in the other 2 review status files
#output: a known ppi file
#	 a predicted ppi file
#notes: The review status files of MGI ids and gene names give different mapping of gene names to Uniprot ids. For some genes, their gene names are not mapped to Uniprot ids, but their MGI ids are mapped to Uniprot ids.
#	A review status file contains more protein information than the fasta file of the same genenames or MGI ids. 
import re
import collectfeatures_all_ppi
import collectfeatures_ppi
#import genes

def main():
	ppi_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/i2d.2_9.Public.MOUSE.tab'
	#genesid_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_ids.csv'
	genesid_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/new_lethal_ids.csv'
	#genesid_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/10_ids.csv'
	review_status_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/review_status_file_of_all_protein_coding_genes_by_mapping_mgiid_to_uniprotid'
	review_status_file2 = 'C:/Users/David/Dropbox/datasets/essential genes prediction/review_status_file_of_all_protein_coding_genes_by_mapping_genename_to_uniprotid'
	review_status_file3 = 'C:/Users/David/Dropbox/datasets/essential genes prediction/review_status_file_of_proteins_which_are_not_in_the_other_2_review_status_files'
	#notes: review_status_file_of_all_protein_coding_genes_by_mapping_mgiid_to_uniprotid and review_status_file_of_all_protein_coding_genes_by_mapping_genename_to_uniprotid do not contain 23 uniprot ids. These uniprot ids are mapped to uniprot ids on Uniprot website separately and the reviewed status file is saved as 'review_status_file_of_proteins_which_are_not_in_the_other_2_review_status_files'
	'''
	fasta_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/all_proteins_of_all_genes.fasta'
	pepstats_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/all_proteins_of_all_genes.pepstats'
	'''
	known_ppi_outfile = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/new_lethal_known_ppi.tab'
	predicted_ppi_outfile = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/new_lethal_predicted_ppi.tab'
	ppiL = [line.strip() for line in open(ppi_file)]
	genesidL = [line.strip() for line in open(genesid_file)]
	review_statusL = [line.strip() for line in open(review_status_file)]
	review_statusL2 = [line.strip() for line in open(review_status_file2)]
	review_statusL3 = [line.strip() for line in open(review_status_file3)]
	#fastaL = [line.strip() for line in open(fasta_file)]
	#pepstatsL = [line.strip() for line in open(pepstats_file)]
	known_ppi_databases = set(['BioGrid_Mouse','BIND_Mouse','Chen_PiwiScreen','DIP_Mouse','I2D-c_Fiona_Mouse,I2D-c_Mouse','IntAct_Mouse','INNATEDB_Mouse','KIM_MYC','MGI','MINT_Mouse','WangEScmplx','WangEScmplxlow','WangEScoIP'])
	(known_ppi,predicted_ppi) = collectfeatures_all_ppi.print_statistics(ppiL,genesidL,known_ppi_databases)#known_ppi: key = uniprot id of a longest reviewed protein, value = set of all known ppis (tupes) of the protein e.g. {(p1,p2),(p2,p3),...}
	'''
	#get proteins length from pepstats file	
	genenames = set()
	for line in genesidL[1:len(genesidL)]:
		valsL = line.split(',')
		genename = valsL[0]
		genenames.add(genename)
	genesL = list(genenames)
	entry_name_hash = get_entry_name(genesL,fastaL)#entry_name_hash: key=uniprot id, value=entry name
	proteins_length = {}#key=uniprot id, value=length of protein
	proteins_length = get_proteins_length(proteins_length,known_ppi,entry_name_hash,pepstatsL)
	proteins_length = get_proteins_length(proteins_length,predicted_ppi,entry_name_hash,pepstatsL)
	longest_proteins_known_ppi = get_longest_proteins(known_ppi,proteins_length)
	longest_proteins_predicted_ppi = get_longest_proteins(predicted_ppi,proteins_length)
	'''
	proteins_length = {}#key=uniprot id, value=length of protein
	proteins_length = get_proteins_length(proteins_length,review_statusL)
	proteins_length = get_proteins_length(proteins_length,review_statusL2)
	proteins_length = get_proteins_length(proteins_length,review_statusL3)
	print('--getting known ppis between longest reviewed proteins and longest proteins--')
	longest_proteins_known_ppi = get_longest_proteins(known_ppi,proteins_length)
	print('--getting predicted ppis between longest reviewed proteins and longest proteins--')
	longest_proteins_predicted_ppi = get_longest_proteins(predicted_ppi,proteins_length)
	longest_proteins_known_ppiL = remove_duplicate_and_self_ppi(longest_proteins_known_ppi)
	print('remove duplicate ppi and self ppi done')
	longest_proteins_predicted_ppiL = remove_duplicate_and_self_ppi(longest_proteins_predicted_ppi)
	print('removal duplicate ppi and self ppi done')
	write_to_file(longest_proteins_known_ppiL,known_ppi_outfile)
	write_to_file(longest_proteins_predicted_ppiL,predicted_ppi_outfile)
	
def get_proteins_length(proteins_length,review_statusL):
	#return: proteins_length: key=uniprot id, value=length of protein
	i=0
	while i < len(review_statusL):
		#print(str(i))
		line = review_statusL[i]
		m = re.match('^ID\s+(.+_.+)\s+[ReviwdUnr]+;*\s*(\d+)\s*AA.*$',line)	
		if m:
			uniprot_entry_name = m.group(1)
			uniprot_entry_name = uniprot_entry_name.strip()
			length = m.group(2)
			i += 1
			uniprotids_of_protein_obtained = False
			#uniprotid = 'none'
			protein_ids = set()
			m2_matched = False #If there are more than one lines which look like 'AC   A6MDD3;' or 'AC   Q8R422; Q8BLT6;'. The uniprot id of the protein is the first uniprot id on the first of these lines. So m2 is matched to the first of these lines only.
			while uniprotids_of_protein_obtained == False and i < len(review_statusL):
				line = review_statusL[i]
				#m2 = re.match('^AC\s+([^;]+);.*$',line)#match 'AC   A6MDD3;' or 'AC   Q8R422; Q8BLT6;'
				m2 = re.match('^AC\s+([\w\d\s;]+)',line)
				#m3 = re.match('^//$',line)#end of a protein
				if m2 and m2_matched == False:
					m2_matched = True
					uniprotids = m2.group(1)
					uniprotidsL = uniprotids.split(';')
					for protein_id in uniprotidsL:
						if protein_id != '':
							protein_id = protein_id.strip()
							protein_ids.add(protein_id)
					uniprotids_of_protein_obtained = True					
					i += 1
				else:
					i += 1
			#if uniprotid == 'none':
			if len(protein_ids) == 0:
				print('Uniprot id is not obtained for '+uniprot_entry_name)
			else:
				for proteinid in list(protein_ids):
					if proteins_length.get(proteinid)==None:
						proteins_length[proteinid] = int(length)
		else:
			i += 1
	return proteins_length
	
def remove_duplicate_and_self_ppi(longest_proteins_ppi):
	ppis = set()
	ks = list(longest_proteins_ppi.keys())
	for k in ks:
		ppi = longest_proteins_ppi[k]
		ppis.add(ppi)
	ppisL = collectfeatures_ppi.remove_duplicate_self_ppi(list(ppis))
	return ppisL
	
def write_to_file(ppisL,ppi_outfile):
	#ppisL:  a list of ppis e.g. [(p1,p2),(p3,p4),...,(pk,pl)]
	fw=open(ppi_outfile,'w')
	fw.write('uniprot1\tuniprot2\n')
	for ppi in ppisL:
		fw.write(ppi[0]+'\t'+ppi[1]+'\n')
	fw.close()

def get_longest_proteins(ppi_hash,proteins_length):
	#ppi_hash: key = uniprot id of a longest reviewed protein, value = set of all known ppis (tupes) of the protein e.g. {(p1,p2),(p2,p3),...}
	#return: longest_proteins: key=uniprot id of a longest reviewed protein, value=longest protein among all the proteins in all the interactions with this longest reviewed protein
	longest_proteins={} #key=uniprot id of a longest reviewed protein, value = longest protein among all the proteins in all the interactions with this longest reviewed protein
	#e.g. p1 is a longest reviewed protein retrieved and has interactions with p2, p3, p4 and p5. P5 is the longest among p2, p3, p4 and p5. The interaction between p1 and p5 is stored in longest_proteins
	proteins_with_unknown_length = set()
	longest_reviewed_proteins_with_some_unknown_length_proteins_in_their_ppis = set()
	uniprotids = list(ppi_hash.keys())
	for uniprotid in uniprotids:
		max_length = 0
		max_length_protein_ppi = ''
		ppis = ppi_hash[uniprotid]
		ppisL = list(ppis)
		for ppi in ppisL:
			p1 = ppi[0]
			p2 = ppi[1]
			if p1 == uniprotid:#p1 is a longest reviewed protein, get the length of p2
				if proteins_length.get(p2)!=None:
					l = proteins_length[p2]
				else:
					l = 'unknown'
					#print(p2+' is not in the review status file(s) and has unknown length')
					other_protein = p2
					proteins_with_unknown_length.add(p2)
					longest_reviewed_proteins_with_some_unknown_length_proteins_in_their_ppis.add(uniprotid)
			else:#p2 is a longest reviewed protein, get the length of p1
				if proteins_length.get(p1)!=None:
					l = proteins_length[p1]
				else:
					l = 'unknown'
					#print(p1+' is not in the review status file(s) and has unknown length')
					other_protein = p1
					proteins_with_unknown_length.add(p1)
					longest_reviewed_proteins_with_some_unknown_length_proteins_in_their_ppis.add(uniprotid)
			if l != 'unknown':
				if l > max_length:
					max_length = l
					max_length_protein_ppi = ppi#this ppi contains the longest protein
			else:
				print(other_protein+' is not considered when comparing protein length as it is not in the review status file(s) and has unknown length.')
		if max_length_protein_ppi != '':				
			longest_proteins[uniprotid] = max_length_protein_ppi
		else:#All proteins in all the ppis of this longest reviewed protein have unknown length. Choose the first ppi as the ppi of the longest protein of this longest reviewed protein.
			longest_proteins[uniprotid] = ppisL[0]
			print('All proteins in all the ppis of '+uniprotid+' have unknown length. Choose the first ppi as the ppi of the longest protein of this longest reviewed protein')
	print(str(len(proteins_with_unknown_length))+' proteins are not in the review status file(s) and have unknown length.')
	print(str(len(longest_reviewed_proteins_with_some_unknown_length_proteins_in_their_ppis))+' longest reviewed proteins have at least 1 ppi in which at least 1 protein has unknown length')
	return longest_proteins
'''	
def get_proteins_length(proteins_length,ppis_hash,entry_name_hash,pepstatsL):
	all_ppis = set()
	uniprotids = list(ppis_hash.keys())
	for uniprotid in uniprotids:
		ppis = ppis_hash[uniprotid]
		all_ppi = all_ppis.union(ppis)
	for ppi in list(all_ppis):
		p1 = ppi[0]
		p2 = ppi[1]
		if proteins_length.get(p1) == None:
			entry_name = entry_name_hash[p1]
			length = collectfeatures_ppi.get_protein_length(entry_name,pepstatsL)
			proteins_length[p1] = length
		if proteins_length.get(p2) == None:
			entry_name = entry_name_hash[p2]
			length = collectfeatures_ppi.get_protein_length(entry_name,pepstatsL)
			proteins_length[p2] = length
	return proteins_length
	
def get_entry_name(genesL,fastaL):
	#input: genesL is a list of gene names
	#	fastaL is a list representing fasta file
	entry_name_hash = {}#key=uniprot id, value=entry name
	(proteins,proteins_of_genes_not_to_collect,proteins_with_no_gene_names) = genes.get_proteins(genesL,fastaL)#proteins is a hashtable: key=gene_name, value=set of tuples (protein name, UniProt id, UniProt entry name)
	genenames = list(proteins.keys())
	for genename in genenames:
		proteins_of_gene = proteins[genename]
		for protein in list(proteins_of_gene):
			uniprotid = protein[1]
			entry_name = protein[2]
			entry_name_hash[uniprotid] = entry_name
	return entry_name_hash
'''

if __name__ == "__main__":
	main()

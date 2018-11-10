#program to get all known ppi(protein-protein interactions) and all predicted ppi between longest reviewed proteins encoded by mouse genes and longest proteins of other genes
#
#input: i2d.2_3.Public.MOUSE.tab (PPI file containing known and predicted PPIs of all mouse genes from the I2D database)
#	a csv file containing the Uniprot ids of the longest reviewed proteins of the mouse genes retrieved from Uniprot (e.g. new_lethal_new_viable_ids.csv or uniprot-mapped_unknown_essentiality_genes.pepstats)
#	a review status file obtained by mapping MGI ids of all protein-coding mouse genes (proteincodingmousegenesfromMouseMine.csv) to uniprot ids on uniprot website
#	a review status file obtained by mapping gene names of all protein-coding mouse genes (proteincodingmousegenesfromMouseMine.csv) to uniprot ids on uniprot website
#	a review status file containing the Uniprot ids which are not contained in the other 2 review status files
'''
#	a fasta file containing all the proteins encoded by each mouse gene (the fasta file is obtained by mapping MGI ids or gene names of all the mouse genes to Uniprot ids using Uniprot and downloading the fasta file from Uniprot)
#	a pepstats file containing the properties of all the proteins encoded by each mouse gene (the pepstats file is output by Pepstats tool taking the above fasta file as input)  
#	
'''
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
#import genes

def main():
	ppi_file = '/home/david/Dropbox/datasets/essential genes prediction/i2d.2_9.Public.MOUSE.tab'
	#genesid_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_ids.csv'
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/mitra_new_lethal_ids.csv'
	#genesid_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/10_ids.csv'
	review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/review_status_file_of_all_protein_coding_genes_by_mapping_mgiid_to_uniprotid'
	review_status_file2 = '/home/david/Dropbox/datasets/essential genes prediction/review_status_file_of_all_protein_coding_genes_by_mapping_genename_to_uniprotid'
	review_status_file3 = '/home/david/Dropbox/datasets/essential genes prediction/review_status_file_of_proteins_which_are_not_in_the_other_2_review_status_files'
	#notes: review_status_file_of_all_protein_coding_genes_by_mapping_mgiid_to_uniprotid and review_status_file_of_all_protein_coding_genes_by_mapping_genename_to_uniprotid do not contain 23 uniprot ids. These uniprot ids are mapped to uniprot ids on Uniprot website separately and the reviewed status file is saved as 'review_status_file_of_proteins_which_are_not_in_the_other_2_review_status_files'
	'''
	fasta_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/all_proteins_of_all_genes.fasta'
	pepstats_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/all_proteins_of_all_genes.pepstats'
	'''
	known_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_known_longest_ppi.tab'
	predicted_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_predicted_longest_ppi.tab'
	ppiL = [line.strip() for line in open(ppi_file)]
	genesidL = [line.strip() for line in open(genesid_file)]
	review_statusL = [line.strip() for line in open(review_status_file)]
	review_statusL2 = [line.strip() for line in open(review_status_file2)]
	review_statusL3 = [line.strip() for line in open(review_status_file3)]
	#fastaL = [line.strip() for line in open(fasta_file)]
	#pepstatsL = [line.strip() for line in open(pepstats_file)]
	known_ppi_databases = set(['BioGrid_Mouse','BIND_Mouse','Chen_PiwiScreen','DIP_Mouse','I2D-c_Fiona_Mouse,I2D-c_Mouse','IntAct_Mouse','INNATEDB_Mouse','KIM_MYC','MGI','MINT_Mouse','WangEScmplx','WangEScmplxlow','WangEScoIP'])
	#longest_proteins_of_all_genes = get_longest_proteins_of_all_genes(genesidL,fastaL,pepstatsL)
	longest_proteins_of_all_genes = get_longest_proteins_of_all_genes(genesidL,review_statusL,review_statusL2,review_statusL3)
	#print('longest_proteins_of_all_genes: '+str(longest_proteins_of_all_genes))	
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
	print('after removing duplicate PPIs and self interactions, there are '+str(len(known_ppiL))+' known longest PPIs')
	print('after removing duplicate PPIs and self interactions, there are '+str(len(predicted_ppiL))+' predicted longest PPIs')

def get_longest_proteins_of_all_genes(genesidL,review_statusL,review_statusL2,review_statusL3):
	proteins_hash = {}#proteins_hash: key = a gene name in infile, value = set of tuples (protein name, UniProt id, UniProt entry name)
	proteins_length_hash={}#proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
	longest_proteins_of_all_genes = set()
	proteins_considered = set() #a set of proteins whose lengths have been considered. This happens if more than one gene encodes the same protein
	genenames = set()
	for line in genesidL[1:len(genesidL)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		genenames.add(genename)
	genes_to_collect = list(genenames)
	(proteins_hash,proteins_length_hash) = get_proteins(genes_to_collect,proteins_hash,proteins_length_hash,review_statusL)
	(proteins_hash,proteins_length_hash) = get_proteins(genes_to_collect,proteins_hash,proteins_length_hash,review_statusL2)
	(proteins_hash,proteins_length_hash) = get_proteins(genes_to_collect,proteins_hash,proteins_length_hash,review_statusL3)
	k=1
	for gene in genes_to_collect:
		(longest_proteinsSet,proteins_considered) = get_longest_proteins(gene,proteins_hash,proteins_length_hash,proteins_considered) #longest_proteinsSet = a set of uniprot ids of the longest proteins of a gene
		print('get_longest_proteins of gene '+str(k)+' done')
		longest_proteins_of_all_genes = longest_proteins_of_all_genes.union(longest_proteinsSet)
		k += 1
	return 	longest_proteins_of_all_genes
	
def get_proteins(genes_to_collect,proteins_hash,proteins_length_hash,review_statusL):
	#return: 
	#	proteins_hash: key = a gene name in infile, value = set of tuples (protein name, UniProt id, UniProt entry name)
	#	proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
	genenames_hash = {} #stores all the gene names of a gene encoding each protein. 
			    #key = (protein name, uniprot id, uniprot entry name), value = the set of all the gene names (recommended gene name and other gene names) of the gene which encodes this protein	
	i=0
	while i < len(review_statusL):
		line = review_statusL[i]
		m = re.match('^ID\s+(.+_.+)\s+[ReviwdUnr]+;*\s*(\d+)\s*AA.*$',line)	
		if m:
			uniprot_entry_name = m.group(1)
			uniprot_entry_name = uniprot_entry_name.strip()
			length = m.group(2)
			genenames = set()
			i += 1
			all_details_of_protein_obtained = False
			protein_name = 'none'
			uniprot_id = 'none'
			gene_has_a_recommended_name = False
			m2_matched = False #If there are more than one lines which look like 'AC   A6MDD3;' or 'AC   Q8R422; Q8BLT6;'. The uniprot id of the protein is the first uniprot id on the first of these lines. So m2 is matched to the first of these lines only.
			#some genes have more than one MGI ids. e.g.
			#
			#DR   MGI; MGI:1196450; Tmem254a.
			#DR   MGI; MGI:3710397; Tmem254b.
			#DR   MGI; MGI:3711260; Tmem254c.
			while all_details_of_protein_obtained == False and  i < len(review_statusL):
				line = review_statusL[i]
				m2 = re.match('^AC\s+([^;]+);.*$',line)#AC   A6MDD3; or AC   Q8R422; Q8BLT6;
				m3 = re.match('^DE\s+[RecNamSub]+\:\s+Full=(.+)$',line)#DE   RecName: Full=Band 4.1-like protein 4A;
				m4 = re.match('^GN\s+Name=(.+);\s+Synonyms=(.+)$',line)#GN   Name=Epb41l4a; Synonyms=Epb4.1l4, Epb4.1l4a, Epb41l4; or GN   Name=Eef2kmt; Synonyms=Fam86, Fam86a;
				m5 = re.match('^GN\s+Name=([^\{\};]+)\{*.*\}*;*$',line)#GN   Name=Epb4.1l4a {ECO:0000313|MGI:MGI:103007};
				m6 = re.match('^GN\s+Synonyms=([^\{\};,]+).*$',line)#GN   Synonyms=Tmem254b {ECO:0000313|Ensembl:ENSMUSP00000098369},
				m7 = re.match('GN\s+([^\{\};,]+).*$',line)#GN   Tmem254c {ECO:0000313|Ensembl:ENSMUSP00000098381};
				m8 = re.match('^GN\s+ORFNames=([^\{\};]+)\{*.*\}*;*$',line)#GN   ORFNames=MNCb-2990; or GN   ORFNames=mCG_116815 {ECO:0000313|EMBL:EDL02759.1};
				m9 = re.match('^DR\s+MGI;\s+MGI:\d+;\s+(.+)$',line)#DR   MGI; MGI:103007; Epb4.1l4a.
				m10 = re.match('^//$',line)#end of a protein
				if m2 and m2_matched == False:
					m2_matched = True
					uniprot_id = m2.group(1)
					proteins_length_hash[(uniprot_entry_name,uniprot_id)] = int(length)
					i += 1
				elif m3:
					protein_name = m3.group(1)
					protein_name = protein_name.strip(';')
					i += 1
				elif m4:
					gene_has_a_recommended_name = True
					genename = m4.group(1)
					genename = genename.strip(';')
					genenames.add(genename.lower())
					other_genenames = m4.group(2)
					other_genenames = other_genenames.strip(';')
					other_genenamesL = other_genenames.split(',')
					for other_genename in other_genenamesL:
						genenames.add(other_genename.lower())
					i += 1
				elif m5:
					gene_has_a_recommended_name = True
					genename = m5.group(1)
					genename = genename.strip()
					genenames.add(genename.lower())
					i += 1
				elif m6:
					gene_has_a_recommended_name = True
					genename = m6.group(1)
					genename = genename.strip()
					genenames.add(genename.lower())
					i += 1	
				elif m7:
					gene_has_a_recommended_name = True
					genename = m7.group(1)
					genename = genename.strip()
					genenames.add(genename.lower())
					i += 1				
				elif m8:
					genename = m8.group(1)
					genename = genename.strip()
					genenames.add(genename.lower())
					i += 1	      					
				elif m9:
					genename = m9.group(1)
					genename = genename.strip('.')
					genenames.add(genename.lower())							
					i += 1
				elif m10:
					all_details_of_protein_obtained = True
					if genenames_hash.get((protein_name,uniprot_id,uniprot_entry_name)) != None:
						genenames2 = genenames_hash[(protein_name,uniprot_id,uniprot_entry_name)]
						genenames2 = genenames2.union(genenames)
						genenames_hash[(protein_name,uniprot_id,uniprot_entry_name)] = genenames2
						print(uniprot_id+' is encoded by numerous genes in the review_status file. The gene names of all the genes are stored in a set.')
					else:
						genenames_hash[(protein_name,uniprot_id,uniprot_entry_name)] = genenames					
				else:
					i += 1
			if gene_has_a_recommended_name == False:
				print('The gene encoding '+uniprot_entry_name+' does not have a recommeded gene name in Uniprot. It has '+str(len(genenames))+' other gene name(s).')
			if uniprot_id == 'none':
				print('Uniprot id is not obtained for '+uniprot_entry_name)
			if protein_name == 'none':
				print('Protein name is not obtained for '+uniprot_entry_name)
			if genenames == set():
				print('Gene names are not obtained for '+uniprot_entry_name)
		else:
			i += 1
	ks = list(genenames_hash.keys())
	for gene_to_collect in genes_to_collect:
		for k in ks:	
			genenames = genenames_hash[k]
			if gene_to_collect in genenames:
				if proteins_hash.get(gene_to_collect)!= None:
					proteins = proteins_hash[gene_to_collect]
					proteins.add(k)	
					proteins_hash[gene_to_collect] = proteins
				else:
					proteins = set()
					proteins.add(k)
					proteins_hash[gene_to_collect] = proteins
	return (proteins_hash,proteins_length_hash)
'''
def get_longest_proteins_of_all_genes(genesidL,fastaL,pepstatsL):
	genenames = set()
	longest_proteins_of_all_genes = set()
	proteins_considered = set() #a set of proteins whose lengths have been considered. This happens if more than one gene encodes the same protein
	for line in genesidL[1:len(genesidL)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		genenames.add(genename)
	genesL = list(genenames)
	(proteins,proteins_of_genes_not_to_collect,proteins_with_no_gene_names) = genes.get_proteins(genesL,fastaL)#proteins is a hashtable: key=gene_name, value=set of tuples (protein name, UniProt id, UniProt entry name)
	k=1
	for gene in genesL:
		(longest_proteinsSet,proteins_considered) = get_longest_proteins(gene,proteins,pepstatsL,proteins_considered) #longest_proteinsSet = a set of uniprot ids of the longest proteins of a gene
		print('get_longest_proteins of gene '+str(k)+' done')
		longest_proteins_of_all_genes = longest_proteins_of_all_genes.union(longest_proteinsSet)
		k += 1
	return 	longest_proteins_of_all_genes
'''	
def retrieve_known_and_predicted_ppi(ppiL,genesidL,longest_proteins_of_all_genes,known_ppi_databases):
	#Retrieve all the known ppis and all the predicted ppis between each retrieved longest reviewed protein and all longest proteins of other genes
	known_ppis = set()#set of known ppis between longest reviewed proteins and longest proteins of other genes
			  #{(p1,p2),(p2,p3),...,(pl,pk)}
	predicted_ppis = set()#predicted ppis between longest reviewed proteins and longest proteins of other genes
			  #{(p1,p2),(p2,p3),...,(pl,pk)}
	known_ppi_dbs = set()
	retrieved_longest_reviewed_proteins = get_uniprotids(genesidL)
	#loop:
	#	Take a PPI (pi,pj) from ppiL (file containing all the PPIs) 
	#	If pi is a retrieved longest reviewed protein and pj is a longest protein of another gene
	#          or
	#   	   pj is a retrieved longest reviewed protein and pi is a longest protein of another gene
	#	Then If (pi,pj) is a known PPI
	#     	     Then keep (pi,pj) as a known PPI
	#     	     Else keep (pi,pj) as a predicted PPI	
	#until ppiL is an empty list
	for known_ppi_db in list(known_ppi_databases):
		known_ppi_db = known_ppi_db.lower()
		known_ppi_dbs.add(known_ppi_db)
	for ppi in ppiL[1:len(ppiL)]:
		l = ppi.split('\t')
		db = l[0]
		p1 = l[1]
		p2 = l[2]
		db = db.lower()
		if p1 in retrieved_longest_reviewed_proteins and p2 in longest_proteins_of_all_genes:
			if db in known_ppi_dbs:
				known_ppis.add((p1,p2))
			else:
				predicted_ppis.add((p1,p2))
		else:
			if p2 in retrieved_longest_reviewed_proteins and p1 in longest_proteins_of_all_genes:
				if db in known_ppi_dbs:
					known_ppis.add((p2,p1))
				else:
					predicted_ppis.add((p2,p1))
	print('all known PPIs between longest reviewed proteins and all longest proteins of genes: ')
	h1 = {}#longest reviewed proteins which have known PPI with a longest protein of another gene
	for ppi in list(known_ppis):
		p1 = ppi[0]
		p2 = ppi[1]
		if p1 in retrieved_longest_reviewed_proteins:
			if h1.get(p1)==None:
				h1[p1] = set()
				h1[p1].add(ppi)
			else:
				h1[p1].add(ppi)
		if p2 in retrieved_longest_reviewed_proteins:
			if h1.get(p2)==None:
				h1[p2] = set()
				h1[p2].add(ppi)
			else:
				h1[p2].add(ppi)
	print(str(len(h1))+' longest reviewed proteins have known PPIs with longest proteins of genes')
	n=0 #total no. of known PPIs including duplicate PPIs and self interactions
	ps = list(h1.keys())
	for p in ps:
		ppis = h1[p]
		n += len(ppis)
		print('longest reviewed protein: '+p+' has '+str(len(ppis))+' known longest PPIs: '+str(ppis))
	print('total no. of known longest PPIs between longest reviewed proteins and all longest proteins of genes (including duplicate PPIs and self interactions): '+str(n))
	h2 = {}#longest reviewed proteins which have predicted PPI with a longest protein of another gene
	for ppi in list(predicted_ppis):
		p1 = ppi[0]
		p2 = ppi[1]
		if p1 in retrieved_longest_reviewed_proteins:
			if h2.get(p1)==None:
				h2[p1] = set()
				h2[p1].add(ppi)
			else:
				h2[p1].add(ppi)
		if p2 in retrieved_longest_reviewed_proteins:
			if h2.get(p2)==None:
				h2[p2] = set()
				h2[p2].add(ppi)
			else:
				h2[p2].add(ppi)
	print(str(len(h2))+' longest reviewed proteins have predicted PPIs with longest proteins of genes')
	n2=0 #total no. of predicted PPIs including duplicate PPIs and self interactions
	ps = list(h2.keys())
	for p in ps:
		ppis = h2[p]
		n2 += len(ppis)
		#print('longest reviewed protein: '+p+' has '+str(len(ppis))+' predicted longest PPIs: '+str(ppis))
	print('total no. of predicted longest PPIs between longest reviewed proteins and all longest proteins of genes (including duplicate PPIs and self interactions): '+str(n2))
	return (known_ppis,predicted_ppis)			
						
def get_uniprotids(genesidL):
	#get the uniprot ids of the retrieved longest reviewed proteins of the genes in the genesidL
	#return: a set of uniprot ids of the retrieved longest reviewed proteins of the genes in the genesidL
	retrieved_longest_reviewed_proteins = set()
	for line in genesidL[1:len(genesidL)]:
		vals = line.split(',')	
		genename = vals[0]
		uniprotid = vals[2]
		retrieved_longest_reviewed_proteins.add(uniprotid)
	return retrieved_longest_reviewed_proteins

def remove_duplicate_self_ppi(ppiL):
	#remove duplicate ppis and self ppis from ppiL
	#input: ppiL (a list of (p1,p2),(p2,p3),(p1,p2),...,(pm,pn))
	#output: ppiL with duplicates and self interactions removed
	duplicates = []
	self_ppi = set()
	no_duplicate_self_ppi = set()
	for ppi in ppiL:
		p1 = ppi[0]
		p2 = ppi[1]
		if (p1,p2) in no_duplicate_self_ppi or (p2,p1) in no_duplicate_self_ppi:
			duplicates.append((p1,p2))
		elif p1 == p2:
			self_ppi.add((p1,p2))
		else:#The ppi has not been collected and is not a self interaction e.g. Q6PDY0 Q6PDY0
			no_duplicate_self_ppi.add((p1,p2))
	'''
	print(str(len(duplicates))+' duplicates: '+str(duplicates)+'\n')
	print(str(len(self_ppi))+' self interactions: '+str(self_ppi)+'\n')
	print(str(len(ppiL))+' PPIs before removing duplicate PPIs and self PPIs\n')
	print(str(len(no_duplicate_self_ppi))+' PPIs after removing duplicate PPIs and self PPIs\n')
	'''
	return list(no_duplicate_self_ppi)

def get_longest_proteins(gene,proteins_hash,proteins_length_hash,proteins_considered):
	#Get all the longest proteins of a gene
	#return: a set of uniprot ids of the longest proteins of a gene
	#longest_proteinsSet = a set of uniprot ids of the longest proteins of a gene
	#proteins_hash: key = a gene name in infile, value = set of tuples (protein name, UniProt id, UniProt entry name)
	#proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
	max_l=0
	longest_proteinsSet = set()
	if proteins_hash.get(gene) == None: #proteins hastable does not contain proteins of the gene
		return (set(),proteins_considered)
	else:
		proteins_length = {} #hashtable: key = protein length, value = set of uniprot ids of proteins with that length
		proteins_of_geneSet = proteins_hash[gene]
		proteins_of_geneL = list(proteins_of_geneSet)
		#compute the length of each protein of the gene; store the lengths of proteins in the hashtable proteins_length and find the length of the longest proteins
		for protein_of_gene in proteins_of_geneL:
			uniprotid = protein_of_gene[1]
			entry_name = protein_of_gene[2]
			if uniprotid not in proteins_considered:
				proteins_considered.add(uniprotid)
				length = proteins_length_hash[(entry_name,uniprotid)]
				if proteins_length.get(length)==None:
					proteins_length[length] = set([uniprotid])
				else:
					proteins = proteins_length[length]
					proteins.add(uniprotid)
					proteins_length[length] = proteins
		lengthL = list(proteins_length.keys())
		lengthL.sort(reverse=True)
		max_length = lengthL[0]
		longest_proteinsSet = proteins_length[max_length]
		print('gene '+gene+' has longest proteins: '+str(longest_proteinsSet)+' length: '+str(max_length))
		return (longest_proteinsSet,proteins_considered)
'''		
def get_longest_proteins(gene,proteins,pepstatsL,proteins_considered):
	#Get all the longest proteins of a gene
	#return: a set of uniprot ids of the longest proteins of a gene
	max_l=0
	longest_proteinsSet = set()
	if proteins.get(gene) == None: #proteins hastable does not contain proteins of the gene
		return (set(),proteins_considered)
	else:
		proteins_length = {} #hashtable: key = protein length, value = set of uniprot ids of proteins with that length
		proteins_of_geneSet = proteins[gene]
		proteins_of_geneL = list(proteins_of_geneSet)
		#compute the length of each protein of the gene; store the lengths of proteins in the hashtable proteins_length and find the length of the longest proteins
		for protein_of_gene in proteins_of_geneL:
			uniprotid = protein_of_gene[1]
			entry_name = protein_of_gene[2]
			if uniprotid not in proteins_considered:
				proteins_considered.add(uniprotid)
				length = get_protein_length(entry_name,pepstatsL)
				if proteins_length.get(length)==None:
					proteins_length[length] = set([uniprotid])
				else:
					proteins = proteins_length[length]
					proteins.add(uniprotid)
					proteins_length[length] = proteins
		lengthL = list(proteins_length.keys())
		lengthL.sort(reverse=True)
		max_length = lengthL[0]
		longest_proteinsSet = proteins_length[max_length]
		return (longest_proteinsSet,proteins_considered)
'''
'''		
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
'''
if __name__ == "__main__":
	main()
	

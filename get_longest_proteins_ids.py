#Find Uniprot ids of longest reviewed proteins of genes from a review status file
#
#input: a csv file containing the gene names of the genes whose longest reviewed proteins are to be collected 
#	a review status file containing the proteins of the genes
#output: Uniprot Ids and Unprot entry names of longest reviewed proteins of the genes
#
#review status file format:
#
#ID   CERS4_MOUSE             Reviewed;         393 AA.
#AC   Q9D6J1; Q8BZA6; Q8C151; Q9CX09;
#
#ID   Q8CDK1_MOUSE            Unreviewed;       855 AA.
#AC   Q8CDK1;
#
#note: A protein can have the same protein name as another protein.
#      Each protein has an unique Uniprot entry name.
#      Pepstats features of proteins are collected using Uniprot entry names
#import genes
import re
import sys

def main():
	'''
	infile='/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/103 features data/new_lethal_genes_not_in_train_set_missing_vals.csv'
	review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_review_status'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_uniprot_ids2.csv'
	'''
	#infile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals_non-mouse_genes_removed_and_duplicates_removed.csv'
	#review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/review_status_obtained_by_mapping_mgiids_to_uniprotids.txt'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids.csv'
	#review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/review_status_obtained_by_mapping_genenames_to_uniprotids.txt'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids2.csv'	
	'''
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_with_missing_uniprotids'
	review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/review_status_obtained_by_mapping_mgiids_of_genenames_with_missing_uniprotids_to_uniprotids.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids3.csv'
	'''
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_with_missing_uniprotids'
	#review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/review_status_obtained_by_mapping_genenames_with_missing_uniprotids_to_uniprotids.txt'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids4.csv'
	'''
	infile='/home/david/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals.csv'
	review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_review_status'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_uniprot_ids2.csv'
	'''
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/58proteins_not_in_fasta_file.ids'
	review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/review_status_obtained_by_mapping_58mgiids_to_uniprotids.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/58proteins_not_in_fasta_file.ids2'
	infileL = [line.strip() for line in open(infile)]
	review_statusL = [line.strip() for line in open(review_status_file)]
	hash_genes = {}#stores the genes in infile i.e. the genes to collect: key = gene name in lower case, value = original gene name		
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_gene_name = valsL[0]
	#for gene in infileL:
		#original_gene_name = gene
		gene_name_lower_case = original_gene_name.lower()
		hash_genes[gene_name_lower_case] = original_gene_name		
	genes_to_collect = list(hash_genes.keys())
	t = len(genes_to_collect) # total no. of genes to collect
	#proteins_hash: key=gene_name, value=set of tuples (protein name, UniProt id, UniProt entry name)
	#review_status_hash: key=(uniprot entry name,uniprot id), value=review status of protein
	#proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
	(proteins_hash,review_status_hash,proteins_length_hash) = get_proteins(genes_to_collect,review_statusL)
	fw = open(outfile,"w")
	fw.write("GeneName,UniProtID,UniprotEntryName,ReviewStatus,ProteinLength,NumberofLongestProteins\n");
	genenamesL = list(proteins_hash.keys())
	for genename in genenamesL:
		genename_lower_case = genename.lower()
		if genename_lower_case in genes_to_collect:
			genes_to_collect.remove(genename_lower_case)
			(longest_proteins_uniprot_entry_names_uniprot_idsL,max_protein_length) = get_longest_proteins(genename,proteins_hash,proteins_length_hash)
			(longest_reviewed_protein_uniprot_entry_name,longest_reviewed_protein_uniprot_id,review_status) = get_a_reviewed_protein(longest_proteins_uniprot_entry_names_uniprot_idsL,review_status_hash)
			fw.write(hash_genes[genename_lower_case]+','+longest_reviewed_protein_uniprot_id+','+longest_reviewed_protein_uniprot_entry_name+','+review_status+','+str(max_protein_length)+','+str(len(longest_proteins_uniprot_entry_names_uniprot_idsL))+'\n')		
	fw.close()
	print('total no. of genes: '+str(t))
	if len(genes_to_collect) > 0:
		print('review_status_file does not have Uniprot ids of '+str(len(genes_to_collect))+' genes.')
		#print('The fasta file does not have protein information of these genes: ')
		#for genename in genes_to_collect:
		#	print(hash_genes[genename])
	else:
		print('Uniprot ids of all genes are collected.')

def get_proteins(genes_to_collect,review_statusL):
	#return: 
	#	proteins_hash: key = a gene name in infile, value = set of tuples (protein name, UniProt id, UniProt entry name)
	#	review_status_hash: key=(uniprot entry name,uniprot id), value=review status of protein
	#	proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
	genenames_hash = {} #stores all the gene names of a gene encoding each protein. 
			    #key = (protein name, uniprot id, uniprot entry name), value = the set of all the gene names (recommended gene name and other gene names) of the gene which encodes this protein	
	proteins_hash = {}#proteins_hash: key = a gene name in infile, value = set of tuples (protein name, UniProt id, UniProt entry name)	
	review_status_hash={}#review_status_hash: key=(uniprot entry name,uniprot id), value=review status of protein
	proteins_length_hash={}#proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
 	i=0
	while i < len(review_statusL):
		line = review_statusL[i]
		m = re.match('^ID\s+(.+_.+)\s+([ReviwdUnr]+);*\s*(\d+)\s*AA.*$',line)	
		if m:
			uniprot_entry_name = m.group(1)
			uniprot_entry_name = uniprot_entry_name.strip()
			review_status = m.group(2)
			length = m.group(3)
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
			while all_details_of_protein_obtained == False:
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
					review_status_hash[(uniprot_entry_name,uniprot_id)] = review_status
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
	''' === display proteins_hash ===
	ks = list(proteins_hash.keys())
	for k in ks:
		print(str(k)+': '+str(proteins_hash[k]))
	'''
	return (proteins_hash,review_status_hash,proteins_length_hash)	

def get_review_status_and_length_of_proteins(review_statusL):
	#review_status_hash: key=(uniprot entry name,uniprot id), value=review status of protein
	#proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
	review_status_hash={}
	proteins_length_hash={}
	i=0
	while i < len(review_statusL):
		line = review_statusL[i]
		m = re.match('^ID\s+(\w+_\w+)\s+Reviewed;*\s*(\d+)\s*AA.*$',line)	
		if m:
			uniprot_entry_name = m.group(1)
			length = m.group(2)
			i += 1
			line = review_statusL[i]
			m3 = re.match('^AC\s+(\w+);.*',line)
			if m3:
				uniprot_id = m3.group(1)
				review_status_hash[(uniprot_entry_name,uniprot_id)]='Reviewed'
				proteins_length_hash[(uniprot_entry_name,uniprot_id)]=int(length)
				i += 1
			else:
				print(line+' does not match pattern in m3')
				sys.exit(-1)
		else:
			m2 = re.match('^ID\s+(\w+_\w+)\s+Unreviewed;*\s*(\d+)\s*AA.*$',line)
			if m2:
				uniprot_entry_name = m2.group(1)
				length = m2.group(2)
				i += 1
				line = review_statusL[i]
				m4 = re.match('^AC\s+(\w+);.*',line)
				if m4:
					uniprot_id = m4.group(1)
					review_status_hash[(uniprot_entry_name,uniprot_id)]='Unreviewed'
					proteins_length_hash[(uniprot_entry_name,uniprot_id)]=int(length)
					i += 1
				else:
					print(line+' does not match pattern in m4')
					sys.exit(-1)
			else:
				i += 1
	return (review_status_hash,proteins_length_hash)

def get_a_reviewed_protein(uniprot_entry_names_uniprot_idsL,review_status_hash):
	#return a reviewed protein or an unreviewed protein if there is no reviewed protein
	a_reviewed_protein = None
	for uniprot_entry_name_uniprot_id in uniprot_entry_names_uniprot_idsL:
		review_status = review_status_hash[uniprot_entry_name_uniprot_id]
		if review_status == 'Reviewed':
			a_reviewed_protein = (str(uniprot_entry_name_uniprot_id[0]),str(uniprot_entry_name_uniprot_id[1]),'Reviewed')
			break
	if a_reviewed_protein != None:		
		return a_reviewed_protein
	else:#if all proteins are unreviewed, return the first protein
		uniprot_entry_name_uniprot_id = uniprot_entry_names_uniprot_idsL[0]
		return (str(uniprot_entry_name_uniprot_id[0]),str(uniprot_entry_name_uniprot_id[1]),'Unreviewed')

def get_longest_proteins(genename,proteins_hash,proteins_length_hash):
	#return: a list of uniprot ids of the longest proteins of a gene
	max_length_so_far=0
	longest_proteinsSet = set()
	proteins_of_geneSet = proteins_hash[genename]
	proteins_of_geneL = list(proteins_of_geneSet)		
	for protein_of_gene in proteins_of_geneL:
		uniprot_id = protein_of_gene[1]			
		uniprot_entry_name = protein_of_gene[2]	
		l = proteins_length_hash[(uniprot_entry_name,uniprot_id)]
		#print('l: '+str(l))
		if l > max_length_so_far:
			max_length_so_far = l
	for protein_of_gene in proteins_of_geneL:		
		uniprot_id = protein_of_gene[1]			
		uniprot_entry_name = protein_of_gene[2]	
		l = proteins_length_hash[(uniprot_entry_name,uniprot_id)]
		if l == max_length_so_far:
			longest_proteinsSet.add((uniprot_entry_name,uniprot_id))
	return (list(longest_proteinsSet),max_length_so_far)
	
if __name__ == "__main__":
	main()
				
	'''=== obsolete ===
	(proteins_hash,proteins_of_genes_not_to_collect,proteins_with_no_gene_names) = genes.get_proteins(genes_to_collect,fasta_fileL)#returns a dictionary: key=gene_name, value=set of tuples (protein name, UniProt id, UniProt entry name)
	print(str(len(proteins_of_genes_not_to_collect))+' proteins in the fasta file are encoded by genes not to be collected.')
	for protein in list(proteins_of_genes_not_to_collect):
		print(protein)
	print(str(len(proteins_with_no_gene_names))+' proteins have no gene names in the fasta file.')
	for protein in list(proteins_with_no_gene_names):
		print(protein)
	#print(str(proteins_hash))
	#review_status_hash: key=(uniprot entry name,uniprot id), value=review status of protein
	#proteins_length_hash: key=(uniprot entry name,uniprot id), value=length of protein
	#(review_status_hash,proteins_length_hash) = get_review_status_and_length_of_proteins(review_statusL)
	'''


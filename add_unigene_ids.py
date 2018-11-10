#Add unigene ids of genes to a csv file after the ensemble ids column by their Uniprot ids, MGI ids, Ensemble ids, gene names in the mapping files and gene names in Mm.data (Unigene database file) in this order.
#
#input: Mm.data
#	new_viable_genes_not_in_train_set.csv (contains fields GeneName,MGI_ID,Uniprot_ID,Ensemble_ID etc)
#	new_viable_genes_genenames_unigeneids.csv
#	new_viable_genes_mgiids_unigeneids.csv
#	new_viable_genes_ensembleids_unigeneids.csv
#	new_viable_genes_uniprot_ids_to_unigene_ids
#output: new_viable_genes_not_in_train_set.csv with a new column for unigene ids
#
#format of new_viable_genes_not_in_train_set.csv:
#
#GeneName,MGI_ID,Uniprot_ID,Ensemble_ID,f1,f2,...,class
#...
#
import sys
import re

def main():
	get_unigene_ids_of_reference_sequence = True
	Mmfile='/home/david/Dropbox/datasets/essential genes prediction/Mm.data/Mm.data'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_unigeneids.csv'
	infile3 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/mgiids_unigeneids.csv'
	infile4 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ensembleids_unigeneids.csv'
	infile5 =  '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/all_uniprotids_unigeneids.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids_unigeneids.csv'
	'''
	Mmfile='/home/david/Dropbox/datasets/essential genes prediction/Mm.data/Mm.data'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_genenames_unigeneids2'
	infile3 = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_mgiids_unigeneids2'
	infile4 = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_ensembleids_unigeneids2'
	infile5 =  '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_uniprotids_unigeneids2'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/103 features data/new_viable_genes_not_in_train_set_missing_vals5.csv'
		
	infile = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/103 features data/new_lethal_genes_not_in_train_set_missing_vals.csv'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/numerous_unigeneids_examples.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_genenames_unigeneids2'
	infile3 = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_mgiids_unigeneids2'
	infile4 = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_ensembleids_unigeneids2'
	infile5 =  '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_uniprotids_unigeneids2'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/103 features data/new_lethal_genes_not_in_train_set_missing_vals2.csv'
	'''
	if get_unigene_ids_of_reference_sequence == True:
		print('Get the unigene ids with reference sequences.\nIf a gene name has more than one unigene id, get the unigene id with a reference sequence and the largest number of ESTs.')
	else:
		print('Get the unigene ids with the largest number of ESTs.')

	MmfileL = [line.strip() for line in open(Mmfile)]
	if len(MmfileL)==0:
		print(Mmfile+" is empty.")
		sys.exit(-1)

	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	
	infile3L = [line.strip() for line in open(infile3)]

	if len(infile3L)==0:
		print(infile3+" is empty.")
		sys.exit(-1)

	infile4L = [line.strip() for line in open(infile4)]

	if len(infile4L)==0:
		print(infile4+" is empty.")
		sys.exit(-1)

	infile5L = [line.strip() for line in open(infile5)]

	if len(infile5L)==0:
		print(infile5+" is empty.")
		sys.exit(-1)
	#A gene name, Uniprot id, MGI id or Ensemble id can map to more than one Unigene ids. One of the mapped Unigene ids is the id of the gene in infile. The other ids are not the id of the gene in infile. In this case, get the Unigene id of the gene in infile.
	#Get the Unigene ids and gene names of all the non-transcribed locus genes in Mm.data
	#hash_non_transcribed_locus: key=unigene id, value=gene name
	hash_non_transcribed_locus = {}#key = unigene id, value = gene name
	hash_non_transcribed_locus_with_max_EST_counts = {}#key=gene name, value = Unigene id of non-transcribed locus with max EST counts
	hash_EST = {}#key = unigene id, value = number of EST
	hash_reference_sequence = {}#key = unigene id, value = reference sequence of gene (e.g. NM_016661.3)
	hash_genenames_unigeneids = {}
	hash_mgiids_unigeneids = {}
	hash_ensembleids_unigeneids = {}
	hash_uniprotids_unigeneids = {}
	genenames_in_infile = set()
	hash_mgiids_genenames_infile = {}
	hash_uniprotids_genenames_infile = {}
	hash_ensembleids_genenames_infile = {}
	hash_genenames_unigeneids_with_max_EST_counts = {}
	hash_mgiids_unigeneids_with_max_EST_counts = {}
	hash_ensembleids_unigeneids_with_max_EST_counts = {}
	hash_uniprotids_unigeneids_with_max_EST_counts = {}
	non_transcribed_locus_mappings_in_Mmfile = set() #set of non-transcribed loci mappings (gene name, unigene id)s in Mm.data
	(hash_non_transcribed_locus,hash_non_transcribed_locus_with_max_EST_counts,hash_reference_sequence,hash_EST,non_transcribed_locus_mappings_in_Mmfile) = get_non_transcribed_locus_and_unigeneids_with_ref_seq_and_EST_count_and_non_transcribed_locus_mappings(MmfileL,hash_non_transcribed_locus,hash_non_transcribed_locus_with_max_EST_counts,hash_reference_sequence,hash_EST,non_transcribed_locus_mappings_in_Mmfile,get_unigene_ids_of_reference_sequence)
	#print('hash_non_transcribed_locus_with_max_EST_counts: '+str(hash_non_transcribed_locus_with_max_EST_counts))
	#print('non_transcribed_locus_mappings_in_Mmfile: '+str(non_transcribed_locus_mappings_in_Mmfile))
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		genename = valsL[0]
		mgiid = valsL[1]
		uniprotid = valsL[2]
		ensembleid = valsL[3]
		genename_lowercase = genename.lower()
		genenames_in_infile.add(genename_lowercase)
		hash_mgiids_genenames_infile[mgiid] = genename_lowercase
		hash_uniprotids_genenames_infile[uniprotid] = genename_lowercase
		hash_ensembleids_genenames_infile[ensembleid] = genename_lowercase
	#get gene name to unigene id mappings
	for line in infile2L[1:len(infile2L)]:
		m = re.match('^(\w+.*\w+)[\t,]+(Mm.\d+)$',line)
		if m:
			gene_name = m.group(1)
			unigene_id = m.group(2)
			if hash_genenames_unigeneids.get(gene_name) != None:
				unigene_idsL = hash_genenames_unigeneids[gene_name]
				unigene_idsL.append(unigene_id)
				hash_genenames_unigeneids[gene_name] = unigene_idsL
			else:
				hash_genenames_unigeneids[gene_name] = [unigene_id]
	gene_namesL = list(hash_genenames_unigeneids.keys())
	for gene_name in gene_namesL:
		unigene_idsL = hash_genenames_unigeneids[gene_name]	
		if len(unigene_idsL) == 1:#there is one unigene id,so there is one mgi id to unigene id mapping.
			hash_genenames_unigeneids_with_max_EST_counts[gene_name] = unigene_idsL[0]
		else:
			unigene_id_found = False
			if get_unigene_ids_of_reference_sequence == True:#this gene name maps to numerous unigene ids, take the unigene id with a reference sequence and the largest EST count
				unigeneids_with_ref_seq = set()
				for unigene_id in unigene_idsL:
					if hash_non_transcribed_locus.get(unigene_id)!= None and hash_reference_sequence.get(unigene_id) != None:
						unigeneids_with_ref_seq.add(unigene_id)
				if unigeneids_with_ref_seq != set():
					(hash_genenames_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,list(unigeneids_with_ref_seq),non_transcribed_locus_mappings_in_Mmfile,gene_name,hash_EST,hash_genenames_unigeneids_with_max_EST_counts)
				else:#all unigene ids have no reference sequence, take the unigene id with the largest EST count
					print('All unigene ids '+str(unigene_idsL)+' have no reference sequence, take the unigene id with the largest EST count.')
					(hash_genenames_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,gene_name,hash_EST,hash_genenames_unigeneids_with_max_EST_counts)
			else:#Take the unigene id with the largest EST count
				(hash_genenames_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,gene_name,hash_EST,hash_genenames_unigeneids_with_max_EST_counts)
			if unigene_id_found == False:
				print('Mm.data does not contain mappings of '+gene_name+' to Unigene ids.')
	#get mgi id to unigene id mappings
	for line in infile3L[1:len(infile3L)]:
		m = re.match('^(MGI:\d+)[\t,]+(Mm.\d+)$',line)
		if m:
			mgiid = m.group(1)
			unigene_id = m.group(2)
			if hash_mgiids_unigeneids.get(mgiid) != None:
				unigene_idsL = hash_mgiids_unigeneids[mgiid]
				unigene_idsL.append(unigene_id)
				hash_mgiids_unigeneids[mgiid] = unigene_idsL
			else:
				hash_mgiids_unigeneids[mgiid] = [unigene_id]
	mgiidsL = list(hash_mgiids_unigeneids.keys())
	for mgiid in mgiidsL:
		unigene_idsL = hash_mgiids_unigeneids[mgiid]
		if len(unigene_idsL) == 1:#there is one unigene id,so there is one mgi id to unigene id mapping.
			hash_mgiids_unigeneids_with_max_EST_counts[mgiid] = unigene_idsL[0]
		else:
			unigene_id_found = False
			if get_unigene_ids_of_reference_sequence == True:#this mgi id maps to numerous unigene ids, take the unigene id with a reference sequence and the largest EST count
				unigeneids_with_ref_seq = set()
				for unigene_id in unigene_idsL:
					if hash_non_transcribed_locus.get(unigene_id)!= None and hash_reference_sequence.get(unigene_id) != None:
						unigeneids_with_ref_seq.add(unigene_id)
				if unigeneids_with_ref_seq != set():
					(hash_mgiids_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,list(unigeneids_with_ref_seq),non_transcribed_locus_mappings_in_Mmfile,hash_mgiids_genenames_infile[mgiid],hash_EST,hash_mgiids_unigeneids_with_max_EST_counts)
				else:#all unigene ids have no reference sequence, take the unigene id with the largest EST count
					print('All unigene ids '+str(unigene_idsL)+' have no reference sequence, take the unigene id with the largest EST count.')
					(hash_mgiids_unigeneids_with_max_EST_counts,unigene_id_found) =get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,hash_mgiids_genenames_infile[mgiid],hash_EST,hash_mgiids_unigeneids_with_max_EST_counts)
			else:#this gene name maps to numerous unigene ids, take the unigene id with the largest EST count
				(hash_mgiids_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,hash_mgiids_genenames_infile[mgiid],hash_EST,hash_mgiids_unigeneids_with_max_EST_counts,)
			if unigene_id_found == False:
				print('Mm.data does not contain mappings of '+mgiid+' to Unigene ids.')
	#get ensemble id to unigene id mappings
	for line in infile4L[1:len(infile4L)]:
		m = re.match('^(ENSMUSG:\d+)[\t,]+(Mm.\d+)$',line)
		if m:
			ensembleid = m.group(1)
			unigene_id = m.group(2)
			if hash_ensembleids_unigeneids.get(ensembleid) != None:
				unigene_idsL = hash_ensembleids_unigeneids[ensembleid]
				unigene_idsL.append(unigene_id)
				hash_ensembleids_unigeneids[ensembleid] = unigene_idsL
			else:
				hash_ensembleids_unigeneids[ensembleid] = [unigene_id]
	ensembleidsL = list(hash_ensembleids_unigeneids.keys())
	for ensembleid in ensembleidsL:
		unigene_idsL = hash_ensembleids_unigeneids[ensembleid]
		if len(unigene_idsL) == 1:#there is one unigene id,so there is one ensemble id to unigene id mapping.
			hash_ensembleids_unigeneids_with_max_EST_counts[ensembleid] = unigene_idsL[0]
		else:#this ensemble id maps to numerous unigene ids, take the unigene id of a gene name in infile
			unigene_id_found = False
			if get_unigene_ids_of_reference_sequence == True:
				unigeneids_with_ref_seq = set()
				for unigene_id in unigene_idsL:
					if hash_non_transcribed_locus.get(unigene_id)!= None and hash_reference_sequence.get(unigene_id) != None:
						unigeneids_with_ref_seq.add(unigene_id)
				if unigeneids_with_ref_seq != set():
					(hash_ensembleids_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,list(unigeneids_with_ref_seq),non_transcribed_locus_mappings_in_Mmfile,hash_ensembleids_genenames_infile[ensembleid],hash_EST,hash_ensembleids_unigeneids_with_max_EST_counts)
				else:#all unigene ids have no reference sequence, take the unigene id with the largest EST count
					print('All unigene ids '+str(unigene_idsL)+' have no reference sequence, take the unigene id with the largest EST count.')
					unigene_id_found = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,hash_ensembleids_genenames_infile[ensembleid],hash_EST,hash_ensembleids_unigeneids_with_max_EST_counts)
			else:#this gene name maps to numerous unigene ids, take the unigene id with the largest EST count
				(hash_ensembleids_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,hash_ensembleids_genenames_infile[ensembleid],hash_EST,hash_ensembleids_unigeneids_with_max_EST_counts)
			if unigene_id_found == False:
				print('Mm.data does not contain mappings of '+ensembleid+' to Unigene ids.')
	#get Uniprot id to unigene id mappings
	for line in infile5L[1:len(infile5L)]:
		m = re.match('^(\w+)[\t,]+(Mm.\d+)$',line)
		if m:
			uniprot_id = m.group(1)
			unigene_id = m.group(2)
			if hash_uniprotids_unigeneids.get(uniprot_id) != None:
				unigene_idsL = hash_uniprotids_unigeneids[uniprot_id]
				unigene_idsL.append(unigene_id)
				hash_uniprotids_unigeneids[uniprot_id] = unigene_idsL
			else:
				hash_uniprotids_unigeneids[uniprot_id] = [unigene_id]
	uniprot_idsL = list(hash_uniprotids_unigeneids.keys()) 
	for uniprot_id in uniprot_idsL:	
		unigene_idsL = hash_uniprotids_unigeneids[uniprot_id]
		if len(unigene_idsL) == 1:#there is one unigene id,so there is one uniprot id to unigene id mapping.
			hash_uniprotids_unigeneids_with_max_EST_counts[uniprot_id] = unigene_idsL[0]
		else:#this Uniprot id maps to numerous unigene ids, take the unigene id of a gene name in infile
			unigene_id_found = False
			if get_unigene_ids_of_reference_sequence == True:
				unigeneids_with_ref_seq = set()
				for unigene_id in unigene_idsL:
					if hash_non_transcribed_locus.get(unigene_id)!= None and hash_reference_sequence.get(unigene_id) != None:
						unigeneids_with_ref_seq.add(unigene_id)
				if unigeneids_with_ref_seq != set():
					(hash_uniprotids_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,list(unigeneids_with_ref_seq),non_transcribed_locus_mappings_in_Mmfile,hash_uniprotids_genenames_infile[uniprot_id],hash_EST,hash_uniprotids_unigeneids_with_max_EST_counts)
				else:#all unigene ids have no reference sequence, take the unigene id with the largest EST count
					(hash_uniprotids_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,hash_uniprotids_genenames_infile[uniprot_id],hash_EST,hash_uniprotids_unigeneids_with_max_EST_counts)
			else:#this gene name maps to numerous unigene ids, take the unigene id with the largest EST count
				(hash_uniprotids_unigeneids_with_max_EST_counts,unigene_id_found) = get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,hash_uniprotids_genenames_infile[uniprot_id],hash_EST,hash_uniprotids_unigeneids_with_max_EST_counts)
			if unigene_id_found == False:
				print('Mm.data does not contain mappings of '+uniprot_id+' to Unigene ids.')	
	fw=open(outfile,'w')
	#fieldsL = infileL[0].split(',')
	#write feature names
	#GeneName,MGI_ID,Uniprot_ID,Ensemble_ID,Unigene_ID,f1,f2,f3,...,class
	'''
	fw.write(fieldsL[0]+','+fieldsL[1]+','+fieldsL[2]+','+fieldsL[3]+',Unigene_ID')
	for field in fieldsL[4:len(fieldsL)]:
		fw.write(','+field)
	fw.write('\n')
	'''
	fw.write('GeneName,MGI_ID,Uniprot_ID,Ensemble_ID,Unigene_ID\n')
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		gene_name = valsL[0]
		mgiid = valsL[1]		
		uniprotid = valsL[2]
		ensembleid = valsL[3]
		if hash_uniprotids_unigeneids_with_max_EST_counts.get(uniprotid) != None:
			unigeneid = hash_uniprotids_unigeneids_with_max_EST_counts[uniprotid]
			hash_where_unigeneid_is_found = 'The unigene id is retrieved from hash_uniprotids_unigeneids_with_max_EST_counts.'
		elif hash_mgiids_unigeneids_with_max_EST_counts.get(mgiid) != None:
			unigeneid = hash_mgiids_unigeneids_with_max_EST_counts[mgiid]
			hash_where_unigeneid_is_found = 'The unigene id is retrieved from hash_mgiids_unigeneids_with_max_EST_counts.'
		elif hash_ensembleids_unigeneids_with_max_EST_counts.get(ensembleid) != None:
			unigeneid = hash_ensembleids_unigeneids_with_max_EST_counts[ensembleid]
			hash_where_unigeneid_is_found = 'The unigene id is retrieved from hash_ensembleids_unigeneids_with_max_EST_counts.'
		elif hash_genenames_unigeneids_with_max_EST_counts.get(gene_name)!= None:
			unigeneid = hash_genenames_unigeneids_with_max_EST_counts[gene_name]
			hash_where_unigeneid_is_found = 'The unigene id is retrieved from hash_genenames_unigeneids_with_max_EST_counts.'
		elif hash_non_transcribed_locus_with_max_EST_counts.get(gene_name)!= None:
			unigeneid = hash_non_transcribed_locus_with_max_EST_counts[gene_name]
			hash_where_unigeneid_is_found = 'The unigene id is retrieved from hash_non_transcribed_locus_with_max_EST_counts.'
		else:
			unigeneid = 'none'#unknown unigene id as there are no mapping from gene name or mgi id or ensemble id or uniprot id to the unigene id from Ensemble database, Uniprot database and Unigene database (Mm.data file)
			hash_where_unigeneid_is_found = 'The unigene id is not in any hash tables.'
			print('There are no mappings from gene name or MGI ID or Ensemble ID or Uniprot ID of '+gene_name+' to a Unigene ID in Ensemble, Uniprot and Unigene(Mm.data file) databases.')
		print(gene_name+','+uniprotid+', '+unigeneid+', '+hash_where_unigeneid_is_found)
		line_new = gene_name#gene name feature
		line_new += ','+mgiid#mgi id feature
		line_new += ','+uniprotid#Uniprot id feature
		line_new += ','+ensembleid#ensemble id feature
		line_new += ','+unigeneid#unigene id feature
		'''
		for val in valsL[4:len(valsL)]:
			line_new += ','+val
		'''
		fw.write(line_new+'\n')
	fw.close()

def get_non_transcribed_locus_and_unigeneids_with_ref_seq_and_EST_count_and_non_transcribed_locus_mappings(MmfileL,hash_non_transcribed_locus,hash_non_transcribed_locus_with_max_EST_counts,hash_reference_sequence,hash_EST,non_transcribed_locus_mappings_in_Mmfile,get_unigene_ids_of_reference_sequence):
	hash_non_transcribed_locus2 = {}
	unigeneid_of_this_gene = ''
	EST_count = 0
	this_gene_is_a_transcribed_locus = True
	for line in MmfileL[1:len(MmfileL)]:
		m = re.match('^ID\s+(Mm.\d+)\s*$',line)
		if m:
			unigeneid = m.group(1)
			unigeneid_of_this_gene = unigeneid
			EST_count = 0
		else:
			m2 = re.match('^GENE\s+([^\s]+)\s*$',line)
			if m2:
				this_gene_is_a_transcribed_locus = False
				gene_name = m2.group(1)
				hash_non_transcribed_locus[unigeneid_of_this_gene] = gene_name
				non_transcribed_locus_mappings_in_Mmfile.add((gene_name,unigeneid_of_this_gene))
				if hash_non_transcribed_locus2.get(gene_name) == None:
					hash_non_transcribed_locus2[gene_name] = [unigeneid_of_this_gene]
				else:
					unigeneidsL = hash_non_transcribed_locus2[gene_name]
					unigeneidsL.append(unigeneid_of_this_gene)
					hash_non_transcribed_locus2[gene_name] = unigeneidsL
			else:
				m3 = re.match('^TITLE\s+Transcribed\s+locus.*$',line)
				if m3:
					this_gene_is_a_transcribed_locus = True
				else:
					m4 = re.match('^.*ACC=(N[MR]+_\d+.\d+).*$',line)#line has the reference sequence of gene e.g. ACC=NM_213727.2, ACC=NR_002866.2
					if m4 and this_gene_is_a_transcribed_locus == False:
						ref_seq = m4.group(1)
						hash_reference_sequence[unigeneid_of_this_gene] = ref_seq
					else:
						m5 = re.match('^.*SEQTYPE=EST.*$',line)	#line is an Expressed Sequence Tag. pattern: 'SEQTYPE=EST'
						if m5 and this_gene_is_a_transcribed_locus == False:
							EST_count += 1
						else:
							m6 = re.match('^//$',line)#line is end of a EST profile
							if m6 and this_gene_is_a_transcribed_locus == False:
								hash_EST[unigeneid_of_this_gene] = EST_count
	#get the unigene ids with reference sequence and largest EST counts 
	gene_namesL = list(hash_non_transcribed_locus2.keys())
	for gene_name in gene_namesL:
		unigene_idsL = hash_non_transcribed_locus2[gene_name]
		if len(unigene_idsL) == 1:#there is one unigene id, so there is one gene name to unigene id mapping.
			hash_non_transcribed_locus_with_max_EST_counts[gene_name] = unigene_idsL[0]
		else:
			unigene_id_found = False
			if get_unigene_ids_of_reference_sequence == True:#this gene name maps to numerous unigene ids, take the unigene id with a reference sequence and the largest EST count
				unigeneids_with_ref_seq = set()
				for unigene_id in unigene_idsL:
					if hash_reference_sequence.get(unigene_id) != None:
						unigeneids_with_ref_seq.add(unigene_id)
				if unigeneids_with_ref_seq != set():
					(hash_non_transcribed_locus_with_max_EST_counts,unigene_id_found) = get_unigeneid_with_largest_EST_count(unigene_id_found,list(unigeneids_with_ref_seq),gene_name,hash_EST,hash_non_transcribed_locus_with_max_EST_counts)
				else:#all unigene ids have no reference sequence, take the unigene id with the largest EST count
					print('All unigene ids '+str(unigene_idsL)+' have no reference sequence, take the unigene id with the largest EST count.')
					(hash_non_transcribed_locus_with_max_EST_counts,unigene_id_found) = get_unigeneid_with_largest_EST_count(unigene_id_found,unigene_idsL,gene_name,hash_EST,hash_non_transcribed_locus_with_max_EST_counts)
					print('hash_non_transcribed_locus_with_max_EST_counts: '+str(hash_non_transcribed_locus_with_max_EST_counts))
			else:#Take the unigene id with the largest EST count
				(hash_non_transcribed_locus_with_max_EST_counts,unigene_id_found) = get_unigeneid_with_largest_EST_count(unigene_id_found,unigene_idsL,gene_name,hash_EST,hash_non_transcribed_locus_with_max_EST_counts)
				print('hash_non_transcribed_locus_with_max_EST_counts: '+str(hash_non_transcribed_locus_with_max_EST_counts))
			if unigene_id_found == False:
				print('In function get_non_transcribed_locus_and_unigeneids_with_ref_seq_and_EST_count_and_non_transcribed_locus_mappings:')
				print('Mm.data does not contain mappings of '+gene_name+' to Unigene ids.')
	return (hash_non_transcribed_locus,hash_non_transcribed_locus_with_max_EST_counts,hash_reference_sequence,hash_EST,non_transcribed_locus_mappings_in_Mmfile)

def get_unigeneid_with_largest_EST_count(unigene_id_found,unigene_idsL,genename_infile,hash_EST,hash_unigeneids_with_max_EST_counts):
	#return hash_unigeneids_with_max_EST_count: key = gene name, value = unigene id of gene with max EST count
	#	unigene_id_found
	max_EST_count = -1
	for unigene_id in unigene_idsL:
		unigene_id_found = True
		EST_count = hash_EST[unigene_id]
		if EST_count > max_EST_count:
			max_EST_count = EST_count
			hash_unigeneids_with_max_EST_counts[genename_infile] = unigene_id
	return (hash_unigeneids_with_max_EST_counts,unigene_id_found)

def get_unigeneid_of_non_transcribed_locus_with_largest_EST_count(unigene_id_found,unigene_idsL,non_transcribed_locus_mappings_in_Mmfile,genename_infile,hash_EST,hash_unigeneids_with_max_EST_counts):
	#return hash_unigeneids_with_max_EST_count: key = mgi id or ensemble id or uniprot id, value = unigene id with max EST count
	#	unigene_id_found
	max_EST_count = -1
	for unigene_id in unigene_idsL:
		if (genename_infile,unigene_id) in non_transcribed_locus_mappings_in_Mmfile:
			unigene_id_found = True
			EST_count = hash_EST[unigene_id]
			if EST_count > max_EST_count:
				max_EST_count = EST_count
				hash_unigeneids_with_max_EST_counts[genename_infile] = unigene_id
	return (hash_unigeneids_with_max_EST_counts,unigene_id_found)

if __name__ == "__main__":
	main()
	 	


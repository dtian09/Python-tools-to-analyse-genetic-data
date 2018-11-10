#Collect Peptstats features of the Uniprot ids in a gene ids file from the fasta file and pepstats file of the Uniprot ids
#
#input:  a fasta file containing protein sequences e.g. uniprot-mapped_unknown_essentiality_genes.fasta
#	 a pepstats output file containing properties of the protein sequences e.g. uniprot-mapped_unknown_essentiality_genes.pepstats
#output: a data table
#
#pepstats file format:
#
#start of a protein sequence: PEPSTATS of A2AHY8_MOUSE from 1 to 747
#
#some genes encode more than one proteins, the properties of the longest protein are collected
#
#note: A protein can have the same protein name as another protein.
#      A protein has an unique Uniprot entry name.
import genes
import re
import sys
import get_longest_proteins_ids

def main():
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	#fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/uniprot_unknowngenes.fasta'
	#pepstats_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/unknowngenes.pepstats'
	#review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/review_status_obtained_by_mapping_genenames_to_uniprotids.txt'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/unknowngenes.csv'
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/58proteins_not_in_fasta_file.ids'
	#fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/58proteins_not_in_fasta_file.fasta'
	#pepstats_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/58proteins_not_in_fasta_file.pepstats'
	#review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/review_status_obtained_by_mapping_58mgiids_to_uniprotids.txt'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/58proteins_not_in_fasta_file.csv'
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/33proteins_not_in_fasta_file.ids'
	fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/33proteins_not_in_fasta_file.fasta'
	pepstats_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/33proteins_not_in_fasta_file.pepstats'
	review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/review_status_obtained_by_mapping_33uniprotids_to_uniprotids.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/33proteins_not_in_fasta_file.csv'
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_ids.csv'
	#fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable.fasta'
	#pepstats_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable.pepstats'
	#review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/review_status_obtained_by_mapping_mgiids_to_uniprotids.txt'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/pepstats.csv'
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/31ids.csv'
	#fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/31.fasta'
	#pepstats_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/31.pepstats'
	#review_status_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/review_status_obtained_by_mapping_31genenames_to_uniprotids.txt'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/31.csv'
	number_of_features=30 #no. of features to collect from pepstats output file
	proteins_position_length ={} #dictionary: key=Uniprot entry name, value= [position of protein sequence, length of protein sequence]
	genesid_fileL = [line.strip() for line in open(genesid_file)]
	genenames = set()
	proteins_to_collect = {} #proteins whose Uniprot ids are not in fasta file
	hash_ids = {}
	original_genenames = {}
	for line in genesid_fileL[1:len(genesid_fileL)]:
		valsL=line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		original_genenames[genename]=original_genename
		genenames.add(genename)
		uniprotid = valsL[2]
		hash_ids[genename] = uniprotid
	pepstats = [line.strip() for line in open(pepstats_file)]
	fasta_fileL = [line.strip() for line in open(fasta_file)]
	review_statusL = [line.strip() for line in open(review_status_file)]
	fw = open(outfile,"w")
	fw.write("GeneName,UniProt_ID,MW,ProteinLength,IP,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Aliphatic,Aromatic,NonPolar,Polar,Charged,basic,Acidic\n");
	(proteins,uniprotids_entrynames,proteins_of_genes_not_to_collect,proteins_with_no_gene_names,genes_not_collected) = genes.get_proteins(list(genenames),fasta_fileL)#proteins: a hashtable with key=gene_name, value=set of tuples (protein name, UniProt id, UniProt entry name)	
	genesL = list(hash_ids.keys())
	i=1 #ith gene that is writen to file
	for gene in genesL:
		print(gene)
		#get the Uniprot entry name (from fasta file) of the Uniprot id of the gene ids file
		uniprotid = hash_ids[gene]
		properties_of_protein_obtained = False
		if uniprotids_entrynames.get(uniprotid)!=None:
			uniprot_entry_name = uniprotids_entrynames[uniprotid]
			(properties,proteins_position_length) = get_pepstats_properties(proteins_position_length,number_of_features,uniprot_entry_name,pepstats)
			properties_of_protein_obtained = True
			if i<len(genesL):
				fw.write(original_genenames[gene]+','+uniprotid+','+properties+'\n')
			else:
				fw.write(original_genenames[gene]+','+uniprotid+','+properties)
			i+=1
		else:#this Uniprot id is not in fasta file
			proteins_to_collect[gene]=uniprotid
	print(str(len(proteins_to_collect))+' Uniprot id(s) they are not in fasta file')
	####For each protein whose uniprot id is not in the fasta file, collect the pepstats features of the protein using a different Uniprot id in the review status file of that protein
	genesL = list(proteins_to_collect.keys())
	for gene in genesL:
		print(original_genenames[gene]+','+proteins_to_collect[gene])
	#get all uniprot ids of each protein from review status file(a protein can have more than 1 uniprot id, but it has 1 uniprot entry name)
	#(proteins_hash,review_status_hash,proteins_length_hash) = get_longest_proteins_ids.get_proteins(list(genenames),review_statusL)#review_status_hash: key=(uniprot entry name,uniprot id), value=review status of protein
	#print(review_status_hash)
	alluniprotids = get_alluniprotids_of_proteins(review_statusL)#get all uniprot ids of each protein. alluniprotids = [{uniprotid1,uniprotid2},{uniprotid3,uniprotid4,uniprotid5,uniprotid6}...]
	print(alluniprotids)
	i=1
 	for gene in genesL:
		print(gene)
		uniprotid = proteins_to_collect[gene]
		properties_of_protein_obtained = False
		uniprotid_eq=''
		for alluniprotids_of_protein in alluniprotids:#get the other uniprot ids of the protein
			if uniprotid in alluniprotids_of_protein:#the other uniprot ids of the protein are found
				alluniprotids_of_protein.remove(uniprotid)
				for uniprotid2 in list(alluniprotids_of_protein):
					if uniprotids_entrynames.get(uniprotid2)!=None:#a different uniprot id of the protein is in fasta file
						uniprotid_eq = uniprotid2
						uniprot_entry_name = uniprotids_entrynames[uniprotid2]
						(properties,proteins_position_length) = get_pepstats_properties(proteins_position_length,number_of_features,uniprot_entry_name,pepstats)
						properties_of_protein_obtained = True
						del proteins_to_collect[gene]
						break
				break
		if properties_of_protein_obtained == False:
			properties='?'
			j=0
			while j < number_of_features:
				properties +=',?'
				j+=1
		if i<len(genesL):
			if uniprotid_eq !='':
				fw.write(original_genenames[gene]+','+uniprotid_eq+','+properties+'\n')
			else:
				fw.write(original_genenames[gene]+','+uniprotid+','+properties+'\n')
		else:
			if uniprotid_eq !='':
				fw.write(original_genenames[gene]+','+uniprotid_eq+','+properties)
			else:
				fw.write(original_genenames[gene]+','+uniprotid+','+properties)
		i+=1
	fw.close()
	#print uniprot ids of the proteins which have not been collected
	if len(proteins_to_collect) > 0:
		print('After checking all the Uniprot ids of each protein in review status file, '+str(len(proteins_to_collect))+' protein(s) have not been collected as none of their Uniprot ids including the ones below are in fasta file i.e. the fasta file does not contain sequences of these proteins.')
		genesL = list(proteins_to_collect.keys())
		for gene in genesL:
			print(original_genenames[gene]+','+proteins_to_collect[gene])
	else:
		print('Pepstats features of all Unirpot ids are collected.')

def  get_alluniprotids_of_proteins(review_statusL):
	#return a list of sets of all uniprot ids of same proteins. e.g. alluniprotids = [{uniprotid1,uniprotid2},{uniprotid3,uniprotid4,uniprotid5,uniprotid6}...]
	#uniprotid1 and uniprotid2 are ids of same protein. uniprotid3, uniprotid4, uniprotid5 and uniprotid6 refer to another protein.
	alluniprotids=[]
	i=0
	while i < len(review_statusL):
		line = review_statusL[i]
		m = re.match('^ID\s+(.+_.+)\s+[ReviwdUnr]+;*\s*\d+\s*AA.*$',line)	
		if m:#get to a Uniprot entry
			m2_matched = False
			while m2_matched == False:
				line = review_statusL[i]
				m2 = re.match('^AC\s+(.+)$',line)#line is like AC   A6MDD3; or AC   Q8R422; Q8BLT6;
				if m2:
					m2_matched = True
					uniprotidsSet = set()
					uniprotids = m2.group(1)
					uniprotidsL = uniprotids.split(';')
					for uniprotid in uniprotidsL:
						m3 = re.match('^\s*$',uniprotid)
						if m3 == None:#string is not white space
							uniprotid = uniprotid.strip()
							uniprotidsSet.add(uniprotid)
					alluniprotids.append(uniprotidsSet)
				i+=1
		else:
			i+=1
	return alluniprotids
'''
def get_pepstats_features_of_protein(uniprotids_entrynames,uniprotid,proteins_position_length,number_of_features,pepstats):
	if uniprotids_entrynames.get(uniprotid)!=None:
		uniprot_entry_name = uniprotids_entrynames[uniprotid]
		(properties,proteins_position_length) = get_pepstats_properties(proteins_position_length,number_of_features,uniprot_entry_name,pepstats)
	else:#this Uniprot id is not in fasta file
		proteins_to_collect[gene]=uniprotid				
		properties='?'
		j=0
		while j < number_of_features:
			properties +=',?'
			j+=1
	return (properties,proteins_position_length)		
'''
def get_pepstats_properties(proteins_position_length,number_of_features,uniprot_entry_name,pepstats):
#input: a dictionary (key=protein, value=[protein position, protein length])
#	protein name
#	uniprotid of a protein
#       a list with each element a line in a pepstats output file 
#output: properties of the protein
	i=0
	values=''	
	features_collected=0
	(i,proteins_position_length) = get_position_of_sequence(proteins_position_length,uniprot_entry_name,pepstats)
	i+=1	
	while(features_collected < number_of_features):
		line = pepstats[i]		
		#match line: "#Molecular weight = 83265.64  		Residues = 747"
		matchObj = re.match("Molecular\\s+weight\\s+=\\s+([0-9\\.]+)\\s+Residues\\s+=\\s+([0-9]+)\\s*",line)
		#Average Residue Weight  = 111.467 	Charge   = -52.5 
		#matchObj2 = re.match("Average\\s+Residue\\s+Weight\\s+=\\s+[0-9\\.]+\\s+Charge\\s+=\\s+([\\-0-9\\.]+)\\s*",line)
		#match line: "Isoelectric Point = 4.3546"
		matchObj3 = re.match("Isoelectric\\s+Point\\s+=\\s+([0-9\\.]+)\\s*",line)
		#match line: "A = Ala		39		5.221  		0.607"
		matchObj4 = re.match("[A,C-I,K-N,P-T,V,W,Y]\\s=\\s[A-Za-z]+\\s+[0-9]+\\s+([0-9\\.]+)\\s+[0-9\\.]+\\s*",line)
		#match line: Tiny	(A+I+L+V)		190		25.435
		matchObj5 = re.match("Tiny\\s+[\\(\\)\\+A-Z]+\\s+[0-9]+\\s+[0-9\\.]+\\s*",line)
		#match line: Small	(A+I+L+V)		190		25.435
		matchObj6 = re.match("Small\\s+[\\(\\)\\+A-Z]+\\s+[0-9]+\\s+([0-9\\.]+)\\s*",line)
		#match line: Aliphatic	(A+I+L+V)		190		25.435
		#	     Aromatic	(F+H+W+Y)		40		13.468
		#	     Non-polar	(A+C+F+G+I+L+M+P+V+W+Y)	161		54.209
		#	     Polar	(D+E+H+K+N+Q+R+S+T+Z)	52		37.143
		#	     Charged	(B+D+E+H+K+R+Z)		23		16.429
		#	     Basic	(H+K+R)			17		12.143
		#	     Acidic	(B+D+E+Z)		6		 4.286	
		matchObj7 = re.match("[\\-A-Za-z]+\\s+[\\(\\)\\+A-Z]+\\s+[0-9]+\\s+([0-9\\.]+)\\s*",line)
		if matchObj != None:
			mw = matchObj.group(1)
			residues = matchObj.group(2)
			values = mw+','+residues
			features_collected +=2
			i+=1
		elif matchObj3 != None:
			IP = matchObj3.group(1)
			values += ','+IP
			features_collected +=1
			i+=1
		elif matchObj4 != None:
			values += ','+matchObj4.group(1)
			features_collected +=1
			i+=1
		elif matchObj5 != None: #line contains "Tiny"
			i+=1
		elif matchObj6 != None: #line contains "Small"
			i+=1
		elif matchObj7 != None:
			values += ','+matchObj7.group(1)
			features_collected +=1
			i+=1
		else:
			i+=1
		#print (features_collected)
	return (values,proteins_position_length)

'''
		elif matchObj2 != None:
			charge = matchObj2.group(1)
			values += ','+charge
			features_collected +=1
			i+=1
'''

def get_position_of_sequence(proteins_position_length,uniprot_entry_name,pepstats):
	#input: a dictionary (key=Uniprot entry name, value=[protein position, protein length])
	#       protein name
	#       uniprot id
	#       a list with each element a line in a pepstats output file
	#output: position of protein header line in the pepstats output file
	if proteins_position_length.get(uniprot_entry_name)==None:
		position_length = get_protein_position_length(uniprot_entry_name,pepstats)
		proteins_position_length[uniprot_entry_name] = position_length	
		return (position_length[0],proteins_position_length)
	else:
		position_length = proteins_position_length[uniprot_entry_name]
		return (position_length[0],proteins_position_length)

def get_protein_position_length(entry_name,pepstats):
	#input: uniprot entry name of a protein sequence e.g. A2AHY8_MOUSE
	#       a list with each element a line in a pepstats output file
	#output: [position of the protein sequence header line, residues of the protein sequence]
	
	i=0
	position=0#position of header line of the sequence in pepstats list
	residues=''
	protein_found=False
	while(i<len(pepstats)):
		line = pepstats[i]	
		#if protein_found == False:
		#match line: "PEPSTATS of A2AHY8_MOUSE from 1 to 747"
		matchObj = re.match("PEPSTATS of "+entry_name+".+\\s+to\\s+(\d+)",line)
		if matchObj:
			residues = matchObj.group(1)				
			protein_found=True
			position = i		
			i+=1
			break		
		else:
			i+=1
	if protein_found==False:
		print ('Uniprot entry name '+entry_name+' is not found in pepstats output file')
	return [position,int(residues)]	

'''				
def get_length_of_sequence(proteins_position_length,uniprot_entry_name,pepstats):
	#input a dictionary: key=protein, value=[protein position, protein length]
	#      protein name
	#      uniprot id
	#      a list with each element a line in a pepstats output file
	#output: length of protein sequence
	if proteins_position_length.get(uniprot_entry_name)==None:
		proteins_position_length[uniprot_entry_name] = get_protein_position_length(uniprot_entry_name,pepstats)
		position_length = proteins_position_length[uniprot_entry_name]
		#print(protein_name+': '+str(position_length[1]))		
		return position_length[1]
	else:
		position_length = proteins_position_length[uniprot_entry_name]
		#print(protein_name+': '+str(position_length[1]))
		return position_length[1]
'''

if __name__ == "__main__":
	main()
				


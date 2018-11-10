#collect all known PPIs between each longest reviewed protein retrieved and all proteins of other genes
#collect all the PPIs (known and predicted PPIs) between each longest reviewed protein retrieved and all proteins of other genes
#input: i2d.2_3.Public.MOUSE.tab (PPI file containing known and predicted PPIs of all mouse genes from the I2D database)
#	a csv file containing the Uniprot ids of the longest reviewed proteins of the mouse genes retrieved from Uniprot (e.g. new_lethal_new_viable_ids.csv or uniprot-mapped_unknown_essentiality_genes.pepstats)
#output known ppi file
#	kp ppi file
import collectfeatures_ppi

def main():
	ppi_file = '/home/david/Dropbox/datasets/essential genes prediction/i2d.2_3.Public.MOUSE.tab'
	#genesid_file = 'C:/Users/David/Dropbox/datasets/essential genes prediction/test set/new_lethal_new_viable_ids.csv'
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15.csv'
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	#known_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/known_ppi.tab'
	#kp_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/kp_ppi.tab'
	known_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/known_ppi.tab'
	kp_ppi_outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/kp_ppi.tab'
	ppiL = [line.strip() for line in open(ppi_file)]
	genesidL = [line.strip() for line in open(genesid_file)]
	known_ppi_databases = set(['BioGrid_Mouse','BIND_Mouse','Chen_PiwiScreen','DIP_Mouse','I2D-c_Fiona_Mouse,I2D-c_Mouse','IntAct_Mouse','INNATEDB_Mouse','KIM_MYC','MGI','MINT_Mouse','WangEScmplx','WangEScmplxlow','WangEScoIP'])
	(known_ppi,kp_ppi) = collect_ppis_print_statistics(ppiL,genesidL,known_ppi_databases)
	write_ppi_to_files(known_ppi,kp_ppi,known_ppi_outfile,kp_ppi_outfile)

def write_ppi_to_files(known_ppi,kp_ppi,known_ppi_outfile,kp_ppi_outfile):
	all_known_ppi = set()
	#remove duplicate ppis and self ppis from all known ppis
	ks = list(known_ppi.keys())
	for k in ks:
		ppis = known_ppi[k]
		all_known_ppi = all_known_ppi.union(ppis)
	all_known_ppiL = collectfeatures_ppi.remove_duplicate_self_ppi(list(all_known_ppi))
	print('removal done')
	#write known ppis to file		
	fw = open(known_ppi_outfile,'w')
	fw.write('uniprot1\tuniprot2\n')
	for ppi in all_known_ppiL:
		fw.write(ppi[0]+'\t'+ppi[1]+'\n')		
	fw.close()
	print('write to file 1 done')
	#remove duplicate ppis and self ppis from all ppis
	all_ppi = set()
	ks2 = list(kp_ppi.keys())
	for k2 in ks2:
		ppis = kp_ppi[k2]
		all_ppi = all_ppi.union(ppis)
	all_ppiL = collectfeatures_ppi.remove_duplicate_self_ppi(list(all_ppi))
	print('removal done')
	#write predicted ppis to file		
	fw2 = open(kp_ppi_outfile,'w')
	fw2.write('uniprot1\tuniprot2\n')
	for ppi in all_ppiL:
		fw2.write(ppi[0]+'\t'+ppi[1]+'\n')		
	fw2.close()
	print('write to file 2 done')

def add_ppi_to_hashtable(ppi_hash,(p1,p2)):
	if ppi_hash.get(p1)==None:
		ppis = set()
		ppis.add((p1,p2))
		ppi_hash[p1] = ppis					
	else:
		ppis = ppi_hash[p1]
		ppis.add((p1,p2))
		ppi_hash[p1] = ppis
	return ppi_hash

def collect_ppis_print_statistics(ppiL,genesidL,known_ppi_databases):
	#output the following statistics:
	#	no. of longest reviewed proteins with no known PPIs
	#	no. of longest reviewed proteins with known PPIs
	#	minimum no. of known PPIs and maximum no. of known PPIs of longest reviewed proteins
	#	no. of longest reviewed proteins with no predicted PPIs
	#	no. of longest reviewed proteins with predicted PPIs
	#	minimum no. of predicted PPIs and maximum no. of predicted PPIs of longest reviewed proteins
	k = 0 #no. of longest reviewed proteins with 0 known PPIs
	l = 0 #no. of longest reviewed proteins with >0 known PPIs
	m = 0 #no. of longest reviewed proteins with 0 predicted PPIs
	n = 0 #no. of longest reviewed proteins with >0 predicted PPIs
	min_known_ppi = 999999999 #minimum no. of known PPIs
	max_known_ppi = 0 #maximum no. of known PPIs 
	min_kp_ppi = 999999999 #minimum no. of predicted PPIs
	max_kp_ppi = 0 #maximum no. of predicted PPIs
	known_ppi = {} #hashtable: key = uniprot id of a longest reviewed protein, value = set of all known ppis (tupes) of the protein e.g. {(p1,p2),(p2,p3),...} 
	kp_ppi = {} #hashtable: key = uniprot id of a longest reviewed protein, value = set of all ppis (tupes) of the protein e.g. {(p1,p2),(p2,p3),...}
	known_ppi_dbs = set()
	retrieved_longest_reviewed_proteins = collectfeatures_ppi.get_retrieved_longest_reviewed_proteins(genesidL)
	#print(str(retrieved_longest_reviewed_proteins))
	#get all known ppis of each uniprot id in the data set and all ppis of each uniprot id in the data set 
	#loop:
	#	Take a PPI (p1,p2) from ppiL (file containing all the PPIs) 
	#	If p1 is a retrieved longest reviewed protein
	#       Then 
	#	     If (p1,p2) is a known PPI
	#     	     Then values = hashtable1.get(p1)
	#		  values.add((p1,p2))
	#		  hastable1[p1] = values 
	#	     Else values = hashtable2.get(p1)
	#		  values.add((p1,p2))
	#		  hastable2[p1] = values 
	#	Else
	#	     If p2 is a retrieved longest reviewed protein
	#	     Then values = hashtable1.get(p2)
	#		  values.add((p1,p2))
	#		  hastable1[p2] = values 
	#     	     Else values = hashtable2.get(p2)
	#		  values.add((p1,p2))
	#		  hastable2[p2] = values	
	#until ppiL is an empty list
	for known_ppi_db in list(known_ppi_databases):
		known_ppi_db = known_ppi_db.lower()
		known_ppi_dbs.add(known_ppi_db)
	for ppi in ppiL[1:len(ppiL)]:
		ppiL = ppi.split('\t')
		db = ppiL[0]
		p1 = ppiL[1]
		p2 = ppiL[2]
		db = db.lower()
		if p1 in retrieved_longest_reviewed_proteins:#p1 is an uniprot id in the data set
			if db in known_ppi_dbs:#ppi is a known ppi
				known_ppi = add_ppi_to_hashtable(known_ppi,(p1,p2))				
				kp_ppi = add_ppi_to_hashtable(kp_ppi,(p1,p2))
			else:#ppi is a predicted ppi
				kp_ppi = add_ppi_to_hashtable(kp_ppi,(p1,p2))
		else:
			if p2 in retrieved_longest_reviewed_proteins:#p2 is an uniprot id in the data set
				if db in known_ppi_dbs:#ppi is a known ppi
					known_ppi = add_ppi_to_hashtable(known_ppi,(p2,p1))
					kp_ppi = add_ppi_to_hashtable(kp_ppi,(p2,p1))
				else:#ppi is a predicted ppi
					kp_ppi = add_ppi_to_hashtable(kp_ppi,(p2,p1))
	for uniprotid in list(retrieved_longest_reviewed_proteins):
		if known_ppi.get(uniprotid) == None:
			k += 1
		else:
			ppis = known_ppi[uniprotid]
			l += 1
			ppis = known_ppi[uniprotid]
			if len(ppis) < min_known_ppi:
				min_known_ppi = len(ppis)
			if len(ppis) > max_known_ppi:
				max_known_ppi = len(ppis)
	for uniprotid in list(retrieved_longest_reviewed_proteins):
		if kp_ppi.get(uniprotid) == None:
			m += 1
		else:
			ppis = kp_ppi[uniprotid]
			n += 1
			ppis = kp_ppi[uniprotid]
			if len(ppis) < min_kp_ppi:
				min_kp_ppi = len(ppis)
			if len(ppis) > max_kp_ppi:
				max_kp_ppi = len(ppis)	
	print(str(k)+' longest reviewed proteins have 0 known PPIs')			
	print(str(l)+' longest reviewed proteins have >0 known PPIs')
	print('minimum no. of known PPIs of a longest reviewed protein: '+str(min_known_ppi))
	print('maximum no. of known PPIs of a longest reviewed protein: '+str(max_known_ppi))
	print(str(m)+' longest reviewed proteins have 0 KP PPIs')
	print(str(n)+' longest reviewed proteins have >0 KP PPIs')
	print('minimum no. of all PPIs of a longest reviewed protein: '+str(min_kp_ppi))
	print('maximum no. of all PPIs of a longest reviewed protein: '+str(max_kp_ppi))
	return (known_ppi,kp_ppi)
	
if __name__ == "__main__":
	main()

#get the duplicate ids of a file if any

import sys

def main():
	id_indx=0
	#infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_all_uniprotids.csv'	
	#infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_ensembleids.csv'
	
	#infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	#infile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals_non-mouse_genes_removed_and_duplicates_removed.csv'	
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/human essential genes/genes_of_different_essentiality_ensembleids.csv'
	infile = '/home/david/Desktop/result7.csv'
	infileL = [line.strip() for line in open(infile)]
	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	idsReadSoFar = set()#set of ids in lower cases in infile
	duplicate_ids_infile = {}#hashtable to store duplicate ids: key=id, value= number of duplicates of id (first occurrance of id is not counted as a duplicate)
	ids_infile = []
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_id = valsL[id_indx]
		id_lower = original_id.lower()
		ids_infile.append(id_lower)
		if id_lower in idsReadSoFar:
			if duplicate_ids_infile.get(id_lower) != None:
				k = duplicate_ids_infile.get(id_lower)
				k += 1
				duplicate_ids_infile[id_lower] = k
			else:
				duplicate_ids_infile[id_lower] = 1	
		else:
			idsReadSoFar.add(id_lower)
	print('total no. of gene ids including duplicate ids: '+str(len(ids_infile)))
	print('total no. of unique gene ids: '+str(len(idsReadSoFar)))
	if duplicate_ids_infile != {}:
		print(str(len(duplicate_ids_infile))+' genes have duplicate ids.')
		ids = list(duplicate_ids_infile.keys())
		k = 0
		for id_lower in ids:
			print(id_lower+' has '+str(duplicate_ids_infile[id_lower])+' duplicate(s) in infile.')
			k += duplicate_ids_infile[id_lower]
		print('total no. of duplicate ids: '+str(k))
	else:
		print('There are no duplicate gene ids.')
	print('The first occurrence of a gene id is not counted as a duplicate gene id.')

if __name__ == "__main__":
        main()


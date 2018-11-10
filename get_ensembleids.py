#Get Ensemble ids from Ensemble and MGI databases as Ensemble and MGI databases may contain Ensemble ids of different gene names.
#input: merged_gene_protein_features_with_ids_missing_vals.csv 
#	genenames_ensembleids_from_ensemble.txt
#	genenames_ensembleids_from_mgi.txt
#output: genenames_ensembleids.txt
#
#format of genenames_ensembleids_from_ensemble.txt
#
#MGI symbol,Ensembl Gene ID
#mt-Cytb,ENSMUSG00000064370
#mt-Nd6,ENSMUSG00000064368
#mt-Nd5,ENSMUSG00000064367
#
#format of genenames_ensembleids_from_mgi.txt
#
#Input	Input Type	MGI Gene/Marker ID	Ensembl ID
#Plekhg2	current symbol	MGI:2141874	ENSMUSG00000037552
#1700006E09Rik	current symbol	MGI:1922687	ENSMUSG00000010841
#Cers4	current symbol	MGI:1914510	ENSMUSG00000008206
#Cers5	current symbol	MGI:1919199	ENSMUSG00000023021

import sys

def main():
	infile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals_non-mouse_genes_removed_and_duplicates_removed.csv'
	infile2='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_ensembleids_from_ensemble.txt'
	infile3='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_ensembleids_from_mgi.txt'
	outfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_ensembleids2.csv'

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

	hash_ensembleids_infile = {}
	hash_ensembleids_infile2 = {}
	hash_ensembleids_infile3 = {}
	for line in infile2L[1:len(infile2L)]:
		vals = line.split(',')
		genename = vals[0]
		genename = genename.lower()
		ensembleid = vals[1]
		if hash_ensembleids_infile.get(genename) == None:
			if ensembleid != '':
				hash_ensembleids_infile[genename] = ensembleid
		if hash_ensembleids_infile2.get(genename) == None:
			if ensembleid != '':
				hash_ensembleids_infile2[genename] = ensembleid
		#else:
		#	print(genename+' in '+infile2+' has more than one Ensemble id. '+hash_ensembleids_infile2[genename]+' is chosen.')
	for line in infile3L[1:len(infile3L)]:
		vals = line.split('\t')
		genename = vals[0]
		genename = genename.lower()
		input_type = vals[1]
		if input_type == 'current symbol' or input_type == 'old symbol' or input_type == 'synonym' or input_type == 'related synonym':
			if len(vals) == 4:
				ensembleid = vals[3]
				if hash_ensembleids_infile.get(genename) == None:
					hash_ensembleids_infile[genename] = ensembleid
				if hash_ensembleids_infile3.get(genename) == None:
					hash_ensembleids_infile3[genename] = ensembleid
				#else:
				#	print(genename+' in '+infile2+' has more than one Ensemble id. '+hash_ensembleids_infile3[genename]+' is chosen.')
			#else:
			#	print(genename+' does not have Ensemble id.')#e.g. Gm10427	current symbol	MGI:3642375
	#write Ensemble ids to file
	t = len(infileL)-1 #total no. of gene names
	j = 0 #no. of gene names with Ensemble ids
	k = 0 #no. of gene names with no Ensemble ids in MGI and Ensemble
	fw = open(outfile,'w')
	fw.write('GeneName,Ensemble_ID\n')
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		genename2 = genename.lower()		
		if hash_ensembleids_infile.get(genename2)!= None:
			j += 1
			fw.write(genename+','+hash_ensembleids_infile.get(genename2)+'\n')
		else:
			fw.write(genename+',none\n')
			k += 1
	fw.close()
	print('total no. of gene names in infile: '+str(t))
	print(str(j)+' gene names have Ensemble ids.')
	print(str(k)+' gene names do not have Ensemble ids in MGI and Ensemble.')
	print('MGI has '+str(len(hash_ensembleids_infile3))+' Ensemble ids.')
	print('Ensemble has '+str(len(hash_ensembleids_infile2))+' Ensemble ids.')

if __name__ == "__main__":
	main()

	

	


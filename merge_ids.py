#merge gene names, MGI ids, Uniprot ids, Ensemble ids and Unigene ids in different files into one file by gene names
#input: id file1,
#	id file2,
#	id file3,
#	etc
#output: a merged id file

import re
import sys

def main():
	infile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals_non-mouse_genes_removed_and_duplicates_removed.csv'
	infile2='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_ensembleids.csv'
	infile3='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids.csv'
	infile4='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_all_uniprotids.csv'
	outfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'

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

	hash_ensembleids = {}
	hash_mgiids = {}
	hash_uniprotids = {}
	for line in infile2L[1:len(infile2L)]:
		m = re.match('^(.+)[\\t,]+(.+)$',line)#match either tab format or csv format			
		if m:
			genename = m.group(1)
			ensembleid = m.group(2)
			hash_ensembleids[genename]=ensembleid
	for line in infile3L[1:len(infile3L)]:
		m = re.match('^(.+)[\\t,]+(.+)$',line)#match either tab format or csv format			
		if m:
			genename = m.group(1)
			mgiid = m.group(2)
			hash_mgiids[genename]=mgiid
	for line in infile4L[1:len(infile4L)]:
		vals = line.split(',')
		genename = vals[0]
		uniprotid = vals[1]
		hash_uniprotids[genename] = uniprotid
	fw=open(outfile,'w')
	fw.write('GeneName,MGI_ID,Uniprot_ID,Ensemble_ID\n')
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		genename = valsL[0]
		if hash_mgiids.get(genename)!=None:
			mgiid = hash_mgiids[genename]
		else:
			mgiid = 'none'
		if hash_uniprotids.get(genename)!=None:
			uniprotid = hash_uniprotids[genename]
		else:
			uniprotid = 'none'	
		if hash_ensembleids.get(genename)!=None:
			ensembleid = hash_ensembleids[genename]
		else:
			ensembleid = 'none'
		fw.write(genename+','+mgiid+','+uniprotid+','+ensembleid+'\n')
	fw.close()

if __name__ == "__main__":
	main()

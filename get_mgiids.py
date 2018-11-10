#Get mgi ids of genes from MGI, Ensemble and Uniprot databases as these databases may contain MGI ids of different gene names.
#input: 
#	merged_gene_protein_features_with_ids_missing_vals.csv (contain gene names of the genes)
#	genenames_mgiids_from_mgi.txt (downloaded from MGI)
#	genenames_mgiids_from_ensemble.txt (downloaded from Ensemble)
#	genenames_mgiids_from_uniprot.txt (downloaded from Uniprot)
#output: a csv file containing gene names and their MGI ids
#
#format of genenames_mgiids_from_mgi.txt
#
#Input	Input Type	MGI Gene/Marker ID	Symbol	Name	Feature Type
#Plekhg2	current symbol	MGI:2141874	Plekhg2	pleckstrin homology domain containing, family G (with RhoGef domain) member 2	protein coding gene
#1700006E09Rik	current symbol	MGI:1922687	1700006E09Rik	RIKEN cDNA 1700006E09 gene	protein coding gene
#Cers4	current symbol	MGI:1914510	Cers4	ceramide synthase 4	protein coding gene
#Cers5	current symbol	MGI:1919199	Cers5	ceramide synthase 5	protein coding gene
#Cers5	human synonym	MGI:2442564	Cers6	ceramide synthase 6	protein coding gene
#
#format of genenames_mgiids_from_ensemble.txt
#
#MGI symbol,MGI ID
#mt-Nd1,MGI:101787
#mt-Nd2,MGI:102500
#mt-Co1,MGI:102504
#mt-Co2,MGI:102503
#mt-Atp8,MGI:99926
#
#format of genenames_mgiids_from_uniprot.txt
#
#ID   D3YY99_MOUSE            Unreviewed;       205 AA.
#AC   D3YY99;
#DT   20-APR-2010, integrated into UniProtKB/TrEMBL.
#DT   20-APR-2010, sequence version 1.
#DT   24-JUN-2015, entry version 35.
#DE   SubName: Full=Pleckstrin homology domain-containing family G member 2 {ECO:0000313|Ensembl:ENSMUSP00000118217};
#DE   Flags: Fragment;
#...
#DR   MGI; MGI:2141874; Plekhg2.
#
#note: if a gene name has no MGI id in MGI, Ensemble and Uniprot, its MGI id is 'none' for the MGI_ID feature.

import sys
import re

def main():
	infile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals_non-mouse_genes_removed_and_duplicates_removed.csv'
	infile2='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_from_mgi.txt'
	infile3='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_from_ensemble.txt'
	infile4 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_from_uniprot.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids.csv'

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

	hash_mgiids_infile = {}
	hash_mgiids_infile2 = {}
	hash_mgiids_infile3 = {}
	hash_mgiids_infile4 = {}
	for line in infile2L[1:len(infile2L)]:
		vals = line.split('\t')
		genename = vals[0]
		genename = genename.lower()
		input_type = vals[1]
		if input_type == 'current symbol' or input_type == 'old symbol' or input_type == 'synonym' or input_type == 'related synonym':
			mgiid = vals[2]
			if hash_mgiids_infile.get(genename) == None:
				if mgiid != '':
					hash_mgiids_infile[genename] = mgiid
		#	else:
		#		print(genename+' in '+infile2+' has more than one MGI id. '+hash_mgiids_infile2[genename]+' is chosen.')
			if hash_mgiids_infile2.get(genename) == None:
				if mgiid != '':
					hash_mgiids_infile2[genename] = mgiid

	for line in infile3L[1:len(infile3L)]:
		vals = line.split(',')
		genename = vals[0]
		genename = genename.lower()
		mgiid = vals[1]
		if hash_mgiids_infile.get(genename) == None:
			if mgiid != '':
				hash_mgiids_infile[genename] = mgiid
		#else:
		#	print(genename+' in '+infile3+' has more than one MGI id. '+hash_mgiids_infile3[genename]+' is chosen.')
		if hash_mgiids_infile3.get(genename) == None:
			if mgiid != '':
				hash_mgiids_infile3[genename] = mgiid

	for line in infile4L:
		m = re.match('^DR\s+MGI;\s+(MGI:\d+);\s+(.+)$',line)#DR   MGI; MGI:103007; Epb4.1l4a.
		if m:
			genename = m.group(2)
			genename = genename.strip('.')
			genename = genename.lower()
			mgiid = m.group(1)
			if hash_mgiids_infile.get(genename) == None:
				hash_mgiids_infile[genename] = mgiid
		#	else:
		#		print(genename+' in '+infile4+' has more than one MGI id. '+hash_mgiids_infile4[genename]+' is chosen.')
			if hash_mgiids_infile4.get(genename) == None:
				if mgiid != '':
					hash_mgiids_infile4[genename] = mgiid

	#write MGI ids to file
	t = len(infileL)-1 #total no. of gene names
	j = 0 #no. of gene names with MGI ids
	k = 0 #no. of gene names with no MGI ids in MGI, Ensemble and Uniprot
	fw = open(outfile,'w')
	fw.write('GeneName,MGI_ID\n')
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		genename2 = genename.lower()		
		if hash_mgiids_infile.get(genename2)!= None:
			j += 1
			fw.write(genename+','+hash_mgiids_infile.get(genename2)+'\n')
		else:
			fw.write(genename+',none\n')
			k += 1
	fw.close()
	print('total no. of gene names in infile: '+str(t))
	print(str(j)+' gene names have MGI ids.')
	print(str(k)+' gene names do not have MGI ids in MGI, Ensemble and Uniprot.')
	print('MGI has '+str(len(hash_mgiids_infile2))+' MGI ids.')
	print('Ensemble has '+str(len(hash_mgiids_infile3))+' MGI ids.')
	print('Uniprot has '+str(len(hash_mgiids_infile4))+' MGI ids.')

if __name__ == "__main__":
	main()


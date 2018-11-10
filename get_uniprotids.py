#get uniprot ids from genenames_uniprotids.csv and genenames_uniprotids2.csv
#
#format of genenames_uniprotids.csv and genenames_uniprotids2.csv
#
#GeneName,UniProtID,UniprotEntryName,ReviewStatus,ProteinLength,NumberofLongestProteins
#Pycrl,Q9DCC4,P5CR3_MOUSE,Reviewed,274,1
#Gpr182,G3X9R9,G3X9R9_MOUSE,Unreviewed,398,1
#Scrt2,Q8BTH6,SCRT2_MOUSE,Reviewed,311,1
#Neurl3,Q8CJC5,NEUL3_MOUSE,Reviewed,254,1

import sys

def main():
	
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids2.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_all_uniprotids0.csv'
	
	'''
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids3.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_uniprotids4.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_all_uniprotids2.csv'
	'''
	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	hash_genenames_uniprotids = {}
	hash_genenames_uniprotids2 = {}
	hash_genenames_uniprotids3 = {}
	hash_original_genenames = {}
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		uniprotid = valsL[1]
		genename = original_genename.lower()
		hash_genenames_uniprotids[genename] = uniprotid
		hash_original_genenames[genename] = original_genename
		hash_genenames_uniprotids2[genename] = uniprotid
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		uniprotid = valsL[1]
		genename = original_genename.lower()
		hash_genenames_uniprotids[genename] = uniprotid
		hash_original_genenames[genename] = original_genename
		hash_genenames_uniprotids3[genename] = uniprotid
	print(str(len(hash_genenames_uniprotids))+' has Uniprot ids in total.')
	print('infile has '+str(len(hash_genenames_uniprotids2))+' ids.')
	print('infile2 has '+str(len(hash_genenames_uniprotids3))+' ids.')	
	fw = open(outfile,'w')
	fw.write('GeneName,Uniprot_ID\n')
	ks = list(hash_genenames_uniprotids.keys())
	for k in ks:
		fw.write(hash_original_genenames[k]+','+hash_genenames_uniprotids[k]+'\n')
	fw.close()

if __name__ == "__main__":
        main()

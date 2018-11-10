#Get Unigene ids from Uniprot mapping files and write the Unigene ids into a new file that maps uniprot ids to unigene ids
#input: uniprotids_unigeneids.txt
#	uniprotids2_unigeneids.txt
#output: all_uniprotids_unigeneids.txt

import sys
import re

def main():
	infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/uniprotids_unigeneids.txt'
	infile2='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/uniprotids2_unigeneids.txt'
	outfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/all_uniprotids_unigeneids.txt'

	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	
	uniprotids_unigeneids = set()
	k = 0
	for line in infileL[1:len(infileL)]:
		m = re.match('^(\w+)\t(Mm.\d+)$',line)
		if m:
			k += 1
			uniprotid = m.group(1)
			unigeneid = m.group(2)
			uniprotids_unigeneids.add(uniprotid+'\t'+unigeneid)
	l = 0
	for line in infile2L[1:len(infile2L)]:
		m = re.match('^(\w+)\t(Mm.\d+)$',line)
		if m:
			l += 1
			uniprotid = m.group(1)
			unigeneid = m.group(2)
			uniprotids_unigeneids.add(uniprotid+'\t'+unigeneid)
	fw=open(outfile,'w')
	fw.write('From\tTo\n')
	uniprotids_unigeneidsL = list(uniprotids_unigeneids)	
	for uniprotid_unigeneid in uniprotids_unigeneidsL:
		fw.write(uniprotid_unigeneid+'\n')
	fw.close()
	print('total no. of uniprot id to unigene id mappings: '+str(len(uniprotids_unigeneids)))
	print('infile has '+str(k)+' uniprot id to unigene id mappings.')
	print('infile2 has '+str(l)+' uniprot id to unigene id mappings.')
	
if __name__ == "__main__":
	main()




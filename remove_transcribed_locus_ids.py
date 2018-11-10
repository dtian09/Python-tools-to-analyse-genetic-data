#Remove the Unigene ids of the transcribed locus from a mapping file: 
#	uniprot id to unigene id file or MGI id to Unigene id file or gene name to Unigene id file or Ensemble id to unigene id file
#
#format of a transcribed locus in Mm.data:
#
#ID          Mm.393562
#TITLE       Transcribed locus
#
#ID          Mm.441460
#TITLE       Transcribed locus, strongly similar to NP_001026562.1 ATP5B gene product [Gallus gallus]
#
#format of a gene in Mm.data
#
#ID          Mm.272264
#TITLE       Syntaxin 3
#
#GENE        Stx3
#
import sys
import re

def main():
	'''
	mappingfile='c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_unigeneids.csv'
	Mmfile='c:/Users/David/Dropbox/datasets/essential genes prediction/Mm.data/Mm.data'
	outfile='c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_unigeneids_no_transcribed_locus.csv'
	
	mappingfile='c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/mgiids_unigeneids.csv'
	Mmfile='c:/Users/David/Dropbox/datasets/essential genes prediction/Mm.data/Mm.data'
	outfile='c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/mgiids_unigeneids_no_transcribed_locus.csv'
	
	mappingfile='c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ensembleids_unigeneids.csv'
	Mmfile='c:/Users/David/Dropbox/datasets/essential genes prediction/Mm.data/Mm.data'
	outfile='c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ensembleids_unigeneids_no_transcribed_locus.csv'
	'''
	mappingfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/all_uniprotids_unigeneids.txt'
	Mmfile='/home/david/Dropbox/datasets/essential genes prediction/Mm.data/Mm.data'
	outfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/all_uniprotids_unigeneids_no_transcribed_locus.txt'
	
	mappingfileL = [line.strip() for line in open(mappingfile)]
	if len(mappingfileL)==0:
		print(mappingfile+" is empty.")
		sys.exit(-1)
	
	MmfileL = [line.strip() for line in open(Mmfile)]
	if len(MmfileL)==0:
		print(Mmfile+" is empty.")
		sys.exit(-1)

	unigeneids_in_Mmfile = set() #all the Unigene ids in Mm.data
	transcribed_locus = set()
	for line in MmfileL[1:len(MmfileL)]:
		m = re.match('^ID\s+(Mm.\d+)$',line)
		if m:
			unigeneid = m.group(1)
			unigeneids_in_Mmfile.add(unigeneid)
		else:
			m2 = re.match('^TITLE\s+Transcribed\s+locus.*$',line)
			if m2:
				transcribed_locus.add(unigeneid)
	unigeneids_not_in_Mmfile = set()
	transcribed_locus_in_mappingfile = set()
	fw=open(outfile,'w')
	fw.write(mappingfileL[0]+'\n')
	for line in mappingfileL[1:len(mappingfileL)]:
		m = re.match('^.+[\\t,]+(Mm.\d+)$',line)#match either tab format or csv format			
		if m:
			unigeneid = m.group(1)
			if unigeneid not in transcribed_locus:
				fw.write(line+'\n')
			else:
				transcribed_locus_in_mappingfile.add(unigeneid)
			if unigeneid not in unigeneids_in_Mmfile:
				unigeneids_not_in_Mmfile.add(unigeneid)
	fw.close()
	print(str(len(transcribed_locus_in_mappingfile))+' Unigene ids in the mapping file: '+mappingfile+' are transcribed locuses and removed from the mapping file.')
	#for unigeneid in list(transcribed_locus_in_mappingfile):
	#	print(unigeneid)
	#Check if any unigene ids in the mapping file are not in Mm.data file. 
	#If an unigene id is not in Mm.data, it is unknown whether the unigene id is a protein coding gene or not.
	if len(unigeneids_not_in_Mmfile) > 0:
		print(str(len(unigeneids_not_in_Mmfile))+' Unigene ids are not in Mm.data.')
		#for unigeneid in list(unigeneids_not_in_Mmfile):
		#	print(unigeneid)
		priint('It is unknown whether these Unigene ids are protein coding genes or not.')

if __name__ == "__main__":
	main()
	

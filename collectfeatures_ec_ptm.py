#collect EC number and post-translational modification features: Glycoprotein, phosphoprotein, Acetylation and transcription of unknown genes from the keywords field under the Miscellaneous field on Uniprot
#
#input: a ids file containing Uniprot ids of the proteins whose EC numbers and keywords features are to be collected
#	a file containing EC numbers and the post-translational modification features in tab delimited format
#output: a data file containing EC numbers and keywords features in csv format
#
#file format:
#
#Entry	EC number	Keywords
#Q9Z110	2.7.2.11; 1.2.1.41	ATP-binding; Alternative splicing; Amino-acid biosynthesis; Complete proteome; Kinase; Membrane; Mitochondrion; Mitochondrion inner membrane; Multifunctional enzyme; NADP; Nucleotide-binding; Oxidoreductase; Proline biosynthesis; Reference proteome; Transferase
#notes: EC numbers are merged into 6 groups 1, 2, 3, 4, 5, 6 such that 1.2.3.5, 1.4._._ are merged to 1, 2.3.1.5, 2._._._ are merged to 2 etc
#ENSMUSG00000064341,NADH dehydrogenase (ubiquinone) activity,molecular_function
import sys
import re
import collectfeatures_ppi

def main():
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	ec_keywords_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/EC and keywords/uniprot-unknown genes EC_keywords.tab'
	out_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/EC and keywords/EC_keywords.csv'
	ec_keywordsL = [line.strip() for line in open(ec_keywords_file)]
	genesidL = [line.strip() for line in open(genesid_file)]
	uniprotids = collectfeatures_ppi.get_uniprotids(genesidL)
	fw = open(out_file,"w")
	(EC_hash,keywords_hash) = get_EC_keywords(ec_keywordsL[1:len(ec_keywordsL)])
	printdata(uniprotids,EC_hash,keywords_hash,fw)
	fw.close()
	
def get_EC_keywords(data):
	EC_hash = {} #key=uniprot id, value = set of EC groups of a protein
	keywords_hash = {} #key=uniprot id, value = set of keywords of a protein
	missing_ec = 0
	missing_keywords = 0
	for line in data:
		l = line.split('\t')
		#print (str(len(l)))
		uniprotid = l[0]
		if len(l)== 1:#missing EC numbers and missing keywords
			EC_groups = set()
			keywords = set()
			missing_ec += 1
			missing_keywords += 1
			#print(line+' has missing EC number and keywords')
		else:
			EC_nums = l[1]
			keywords = l[2]
			EC_nums = EC_nums.strip()
			keywords = keywords.strip()
			if missing(EC_nums):
				EC_groups = set()
				missing_ec += 1
				#print('EC number missing: '+line)
			else:
				EC_groups = map_EC_nums_to_EC_groups(EC_nums)
			if missing(keywords):
				keywords = set()
				missing_keywords += 1
				#print('keywords missing: '+line)
			else:
				keywords = get_keywords(keywords)		
		EC_hash[uniprotid] = EC_groups
		keywords_hash[uniprotid] = keywords
	print('no. of proteines with missing EC numbers: '+str(missing_ec))
	print('no. of proteines with missing keywords: '+str(missing_keywords))
	return (EC_hash,keywords_hash)

def missing(value):
	m = re.match("^\s*$",value)
	if m:
		return True
	else:
		return False

def map_EC_nums_to_EC_groups(EC_nums):
	#input: a string of EC numbers of a protein
	#output: a set of EC groups of the EC numbers
	#
	#EC numbers format: 2.7.2.11; 4.2.1.41
	
	EC_groups = set()
	EC_numsL= EC_nums.split(';')
	for EC_num in EC_numsL:
		EC_num = EC_num.strip()
		m1 = re.match("1(\.)(.*)",EC_num)
		m2 = re.match("2(\.)(.*)",EC_num)
		m3 = re.match("3(\.)(.*)",EC_num)
		m4 = re.match("4(\.)(.*)",EC_num)
		m5 = re.match("5(\.)(.*)",EC_num)
		m6 = re.match("6(\.)(.*)",EC_num)
		if m1:
			EC_group = 1
		elif m2:
			EC_group = 2
		elif m3:
			EC_group = 3
		elif m4:
			EC_group = 4
		elif m5:
			EC_group = 5
		elif m6:
			EC_group = 6
		else:
			#print(EC_num+' does not map to any EC group.')
			EC_group='none'
			#sys.exit(-1)
		if EC_group != 'none':
			EC_groups.add(EC_group)
	return EC_groups

def get_keywords(keywordsStr):
	#input: a string of keywords of a protein
	#output: a set of the keywords: Glycoprotein, Phosphoprotein, Acetylation and Transcription
	#
	#format of keywords field: Chaperone; Chromatin regulator; Complete proteome; Nucleus; Phosphoprotein; Reference proteome; Transcription; Transcription regulation
	keywordsSet = set()
	keywordsL= keywordsStr.split(';')
	for keyword in keywordsL:
		keyword = keyword.strip().lower()
		if keyword == 'glycoprotein':
			keywordsSet.add('glycoprotein')
		elif keyword == 'phosphoprotein':
			keywordsSet.add('phosphoprotein')
		elif keyword == 'acetylation':
			keywordsSet.add('acetylation')
		elif keyword == 'transcription':
			keywordsSet.add('transcription')
	return keywordsSet

def printdata(uniprotids,EC_hash,keywords_hash,fw):
	fw.write('UniProt_ID,1.-.-.-,2.-.-.-,3.-.-.-,4.-.-.-,5.-.-.-,6.-.-.-,GlycoProtein,Phosphoprotein,Acetylation,Transcription\n')
	for uniprotid in list(uniprotids):
		if EC_hash.get(uniprotid)==None:#Uniprot does not have an entry for this Uniprot id.
			EC_groupsStr = '0,0,0,0,0,0,'
		else:
			EC_groups = EC_hash[uniprotid]
			if len(EC_groups) == 0:#Uniprot has an entry for this Uniprot id, but the EC numbers are missing.
				EC_groupsStr = '0,0,0,0,0,0,'
			else:
				if 1 in EC_groups:
					EC_1='1,'
				else:
					EC_1='0,'
				if 2 in EC_groups:
					EC_2='1,'
				else:
					EC_2='0,'
				if 3 in EC_groups:
					EC_3='1,'
				else:
					EC_3='0,'
				if 4 in EC_groups:
					EC_4='1,'
				else:
					EC_4='0,'
				if 5 in EC_groups:
					EC_5='1,'
				else:
					EC_5='0,'
				if 6 in EC_groups:
					EC_6='1,'
				else:
					EC_6='0,'
				EC_groupsStr = EC_1+EC_2+EC_3+EC_4+EC_5+EC_6
		if keywords_hash.get(uniprotid)==None:#Uniprot does not have an entry for this Uniprot id.
			keywordsStr='0,0,0,0'
		else:
			keywords = keywords_hash[uniprotid]
			if len(keywords) == 0:#Uniprot has an entry for this Uniprot id, but the keywords are missing.
				keywordsStr='0,0,0,0'
			else:
				if 'glycoprotein' in keywords:
					gly='1,'
				else:
					gly='0,'
				if 'phosphoprotein' in keywords:
					pho='1,'
				else:
					pho='0,'
				if 'acetylation' in keywords:
					ace='1,'
				else:
					ace='0,'
				if 'transcription' in keywords:
					tra='1'
				else:
					tra='0'
				keywordsStr = gly+pho+ace+tra			
		fw.write(uniprotid+','+EC_groupsStr+keywordsStr+'\n')

if __name__ == "__main__":
	main()

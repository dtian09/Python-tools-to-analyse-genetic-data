#add age feature to the unknowngenes.csv before the class attribute.
def main():
	genesfile ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes.csv'
	agefile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/age/age.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes2.csv'
	genesL = [line.strip() for line in open(genesfile)]
	ageL = [line.strip() for line in open(agefile)]
	age_hash = {}
	genes_hash = {}
	idshash = {}#contains ensemble ids and gene names of genesfile
	for line in ageL[1:len(ageL)]:
		vals = line.split(',')
		ensembleid = vals[0]
		age_hash[ensembleid]=vals[1]
	for line in genesL[1:len(genesL)]:
		vals = line.split(',')
		genename = vals[0]
		ensembleid = vals[3]
		if idshash.get(ensembleid)==None:
			genenames = set()
			genenames.add(genename)
			idshash[ensembleid] = genenames
		else:
			genenames = idshash[ensembleid] 
			genenames.add(genename)
			idshash[ensembleid] = genenames
	for line in genesL[1:len(genesL)]:
		vals = line.split(',')
		genename = vals[0]
		genes_hash[genename] = vals
	featuresnames = genesL[0]
	add_age_feature(featuresnames,outfile,idshash,genes_hash,age_hash)

def add_age_feature(featuresnames,outfile,idshash,genes_hash,age_hash):
	fw = open(outfile,'w')
	fnamesL = featuresnames.split(',')
	fnames=''
	for fname in fnamesL[0:len(fnamesL)-1]:
		fnames += ','+fname
		fnames = fnames.strip(',')	
	fw.write(fnames+',Age'+','+fnamesL[len(fnamesL)-1]+'\n')
	ensembleids = idshash.keys()
	for ensembleid in ensembleids:
		genenames = idshash[ensembleid]
		for genename in genenames:
			vals = genes_hash[genename]
			if age_hash.get(ensembleid) != None:
				age = age_hash[ensembleid]
			else:
				age = '?'
			gene=''
			for val in vals[0:len(vals)-1]:
				gene += ','+val
			gene += ','+age
			gene += ','+vals[len(vals)-1]#add class
			gene = gene.strip(',')
			fw.write(gene+'\n')
	fw.close()

if __name__ == "__main__":
	main()		


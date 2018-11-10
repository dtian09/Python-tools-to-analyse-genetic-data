#Get the gene ids of the genes in balanced1307_missing_values.csv from training_all_genesinfo.csv
import sys

def main():
	infile='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'#known gene set
	infile2 ='/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_missing_values_as_-1.csv'#train set
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/balanced1307_missing_values_as_-1_ids.csv'
	infileL = [line.strip() for line in open(infile)]
	infile2L = [line.strip() for line in open(infile2)]
	hash1={}
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		genename = valsL[0]
		mgiid = valsL[1]
		uniprotid = valsL[2]
		ensembleid = valsL[3]
		vals = ''
		for val in valsL[4:len(valsL)]:
			vals += ','+val
		vals = vals.strip(',')
		hash1[vals]=[genename,mgiid,uniprotid,ensembleid]
	#print(hash1)
	infile2S = set()
	for line in infile2L[1:len(infile2L)]:
		infile2S.add(line)
	fw = open(outfile,'w')
	fw.write(infileL[0]+'\n')
	for gene in infile2S:
		if hash1.get(gene)!=None:
			idsL = hash1[gene]
			print('got')
		else:
			print(gene)
			sys.exit()
		genename = idsL[0]
		mgiid = idsL[1]
		uniprotid = idsL[2]
		ensembleid = idsL[3]
		instance = ''
		instance = genename+','+mgiid+','+uniprotid+','+ensembleid+','+gene
		fw.write(instance+'\n')
	fw.close()
			
if __name__ == "__main__":
        main()


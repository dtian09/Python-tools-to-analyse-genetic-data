#keep the first occurrance of each gene and remove duplicate genes from a genes file
def main():
	infile = '/home/david/Desktop/result7.csv'
	outfile = '/home/david/Desktop/result8.csv'
	fw=open(outfile,'w')
	infileL = [line.strip() for line in open(infile)]
	genesReadSoFar = set()
	geneSet = set()
	n=0
	for gene in infileL[1:len(infileL)]:
		genename = gene.split(',')[0] 
		if genename in genesReadSoFar:
			n += 1
		else:
			genesReadSoFar.add(genename)
			geneSet.add(gene)
	fw.write(infileL[0]+'\n')
	for gene in geneSet:
		fw.write(gene+'\n')
	fw.close()
	print('infile has '+str(len(infileL)-1)+' genes')
	print('no. of duplicate genes: '+str(n))
	print('outfile has '+str(len(geneSet))+' genes')

if __name__ == "__main__":
	main()




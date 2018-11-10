#get mutant gene names
import re

def main():
	infile = '/home/david/Dropbox/datasets/essential genes prediction/mutant lethal genes/MGIalleleQuery_20160316_125157.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/mutant lethal genes/mutant_genenames'
	infileL = [line.strip() for line in open(infile)]
	fw = open(outfile,'w')
	genes = set()
	for line in infileL[1:len(infileL)]:
		valsL = line.split('\t')
		genename = valsL[1]
		m = re.match('^([\w]+)<.+>$',genename)
		if m:
			genename = m.group(1)
		genes.add(genename)
	for gene in genes:
		print(gene)
		fw.write(gene+'\n')				
	fw.close()
	print('no. of genes: '+str(len(genes)))

if __name__ == "__main__":
	main()

	


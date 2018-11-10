#compute length of genes by (end bp - start bp)+1.
#input: gene_start_end_bps.csv, a csv file containing start bp and end bp of genes
#output: a csv file containing length of genes
#
#format of csv file
#MGI ID,Ensembl Gene ID,Gene Start (bp),Gene End (bp)
#MGI:1329010,ENSMUSG00000000204,83175186,83190221
#MGI:1344333,ENSMUSG00000000266,140664599,140767715
#MGI:107506,ENSMUSG00000000983,83709015,83711348
#MGI:1339467,ENSMUSG00000001025,90612882,90624181
#MGI:2442218,ENSMUSG00000001053,51643063,51650842

import sys

def main():
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/test set/bps.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/gene_length.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/bps.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/bps1013.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length1013.csv'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/bps10.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length10.csv'
	genesL = [line.strip() for line in open(infile)]
	get_length(genesL[1:len(genesL)],outfile)
	
def get_length(genesL,outfile):
	hash_length = {} #key = (MGI id, ensemble id), value = length of gene
	for gene in genesL:
		geneL = gene.split(',')
		mgi_id = geneL[0]
		ensemble_id = geneL[1]
		bp_start = geneL[2]
		bp_end = geneL[3]
		if mgi_id == '' or mgi_id == '?' or mgi_id == 'none':
			print('missing mgi id')
			mgi_id = 'none'
		if ensemble_id == '' or ensemble_id == '?' or ensemble_id == 'none':
			print('missing ensemble id')
			ensemble_id = 'none'
		length = int(bp_end) - int(bp_start) + 1
		hash_length[(mgi_id,ensemble_id)] = length
		if length <= 0:
			print(mgi_id+' ,'+ensemble_id+' : bp end <= bp start; bp end: '+bp_end+', bp start: '+bp_start)
	fw = open(outfile,'w')
	fw.write("MGI_ID,Ensemble_ID,GeneLength\n")
	ks = list(hash_length.keys())
	for k in ks:
		mgi_id = k[0]
		ensemble_id = k[1]
		fw.write(mgi_id+','+ensemble_id+','+str(hash_length[k])+'\n')
	fw.close()

if __name__ == "__main__":
	main()	

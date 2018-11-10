#compute intron length: the sum of length of all introns of a gene by subtracting the total exon length from the gene length
#input: gene_length.csv
#	exon_length.csv
#output: a csv file containing intron length of genes
#
# format of gene_length.csv
#
#MGI_ID,Ensemble_ID,GeneLength
#MGI:3583955,ENSMUSG00000074639,10844
#MGI:2145430,ENSMUSG00000038042,47701
#MGI:1917946,ENSMUSG00000031179,5231
#MGI:1920149,ENSMUSG00000068205,1997658
#MGI:1924781,ENSMUSG00000058589,1099792
#
#format of exon_length.csv
#
#Ensemble_ID,ExonLength
#ENSMUSG00000002012,1536
#ENSMUSG00000047613,1098
#ENSMUSG00000049598,1801
#ENSMUSG00000029518,2820
#ENSMUSG00000021715,2071
import sys

def main():
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/test set/G6b_gene_length.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/test set/G6b_exon_length.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/G6b_intron_length.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/intron_length.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length1013.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length1013.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/intron_length1013.csv'	
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length10.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length10.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/intron_length10.csv'	
	geneslengthL = [line.strip() for line in open(infile)]
	exonslengthL = [line.strip() for line in open(infile2)]
	fw=open(outfile,'w')
	fw.write('MGI_ID,Ensemble_ID,IntronLength\n')
	hash_exons_length = exons_length(exonslengthL[1:len(exonslengthL)])
	(hash_genes_length,hash_ensembleid_mgiid) = genes_length(geneslengthL[1:len(geneslengthL)])
	ensemble_ids = list(hash_exons_length.keys())
	for ensemble_id in ensemble_ids:
		if hash_genes_length.get(ensemble_id) != None:
			gene_length = hash_genes_length[ensemble_id]
		else:
			print(ensemble_id+' is missing in hash_genes_length.')
			gene_length = '?'
		if hash_exons_length.get(ensemble_id) != None:
			exon_length = hash_exons_length[ensemble_id]
		else:
			print(ensemble_id+' is missing in hash_exons_length.')
			exon_length = '?'
		if gene_length != '?' and exon_length != '?':
			intron_length = gene_length - exon_length
		else:
			print('gene length or exon length is missing.')
			intron_length = '?'
		mgi_id = hash_ensembleid_mgiid[ensemble_id]
		fw.write(mgi_id+','+ensemble_id+','+str(intron_length)+'\n')
	fw.close()

def exons_length(exonslengthL):
	hash_exons_length = {}#hash_exons_length: key = ensemble gene id, value = exon length
	for line in exonslengthL:
		valsL = line.split(',')
		gene_id = valsL[0]
		exon_length = valsL[1]
		if exon_length == '?':
			hash_exons_length[gene_id] = exon_length
		else:
			hash_exons_length[gene_id] = int(exon_length)
	return hash_exons_length

def genes_length(geneslengthL):
	hash_genes_length = {}#hash_genes_length: key = ensemble gene id, value = gene length
	hash_ensembleid_mgiid = {}#key = ensemble id of gene, value = mgi id of gene
	for line in geneslengthL:
		valsL = line.split(',')
		mgi_id = valsL[0]
		ensemble_id = valsL[1]
		gene_length = valsL[2]
		if gene_length == '?':
			hash_genes_length[ensemble_id] = '?'
		else:
			hash_genes_length[ensemble_id] = int(gene_length)
		hash_ensembleid_mgiid[ensemble_id] = mgi_id
	return (hash_genes_length,hash_ensembleid_mgiid)

if __name__ == "__main__":
	main()	

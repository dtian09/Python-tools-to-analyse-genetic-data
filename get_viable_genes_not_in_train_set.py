#Get the viable genes from the known gene set which are not in the training set
import sys

def main():
	'''	
	infile='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo_no_ambiguous_genes.csv'#known gene set
	infile2 ='/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/new_lethal_new_viable_balanced_missing_vals.csv'#train set
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/viable_genes_test_set.csv'
	'''
	infile='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'#known gene set
	infile2 ='/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_missing_values.csv'#train set
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/viable_genes_test_set2_2.csv'
	infileL = [line.strip() for line in open(infile)]
	infile2L = [line.strip() for line in open(infile2)]
	train_set = set()
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		genename = valsL[0]
		genename = genename.lower()
		train_set.add(genename)
	#n = 2149
	k = 0 #no. of viable genes which have been retrieved from the known gene set
	fw = open(outfile,'w')
	fw.write(infileL[0]+'\n')
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		gene_class = vals[len(vals)-1]
		#if k == 2149:
		#	break
		if genename.lower() not in train_set and gene_class == 'Viable':
			fw.write(line+'\n')
			k += 1
	fw.close()
	print(str(k)+' viable genes have been retrieved from the known gene set.')
			
if __name__ == "__main__":
        main()


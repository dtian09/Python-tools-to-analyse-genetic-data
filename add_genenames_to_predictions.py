#add gene names to predictions of unknown genes output by classifiers
#input: merged_gene_protein_features_with_ids.csv (unknown genes data set with gene names and ids)
#	prediction.chimerge_significance_level0.95_random_forest_ROC_area=0.998 (prediction file)
#output: prediction_with_gene_names.chimerge_significance_level0.95_random_forest_ROC_area=0.998
#

import sys
import re

def main():
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_44_features_data.random_forest_230_trees'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_44_features_data_with_gene_names.random_forest_230_trees'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_48_features_data.random_forest'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_48_features_data_with_gene_names.random_forest'
	infileL = [line.strip() for line in open(infile)]
	predictL = [line.strip() for line in open(infile2)]
	hash_genenames = get_genenames(infileL)#key=instance no., value=gene name of the instance no.
	add_gene_names(hash_genenames,predictL,outfile)

def get_genenames(infileL):
	instance_no = 1
	hash_genenames={}
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		hash_genenames[instance_no]=genename
		instance_no += 1
	return hash_genenames

def add_gene_names(hash_genenames,predictL,outfile):
	fw = open(outfile,'w')
	i=0
	while i < len(predictL):
		line = predictL[i]
		m = re.match('^(inst#\s+actual\s+predicted\s+error\s+prediction)$',line)
		if m:
			fw.write('gene names  '+line+'\n')			
			break
		else:
			fw.write(line+'\n')
			i +=1	
	while i < len(predictL):
		line = predictL[i]	
		m2 = re.match('^([\d]+)(.+)$',line)#1        1:?   2:Viable       0.617 
		if m2:
			instance_no = int(m2.group(1))
			genename = hash_genenames[instance_no]
			fw.write(genename+'  '+line+'\n')
			i += 1
		else:
			i += 1
	fw.close()

if __name__ == "__main__":
	main()



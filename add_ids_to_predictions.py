#add gene names, MGI Ids and ensemble Ids to predictions of unknown genes output by classifiers
#input: merged_gene_protein_features_with_ids.csv (unknown genes data set with gene names and ids)
#	prediction.chimerge_significance_level0.95_random_forest_ROC_area=0.998 (prediction file)
#output: prediction_with_gene_names.chimerge_significance_level0.95_random_forest_ROC_area=0.998
#
import sys
import re

def main():
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids_unigeneids.csv'
	infile2 = '/home/david/Dropbox/working papers/gene essentiality predictions/supplementary tables/Supplementary Table 4.csv'
	outfile = '/home/david/Dropbox/working papers/gene essentiality predictions/supplementary tables/Supplementary Table 4_ids.csv'
	infileL = [line.strip() for line in open(infile)]
	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	predictL = [line.strip() for line in open(infile2)]
	if len(predictL)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	#hash_genenames_mgiids: key=gene name, value=mgi id
	#hash_genenames_ensembleids: key=gene name, value = ensemble id
	(hash_genenames_mgiids,hash_genenames_ensembleids) = get_genenames_mgiids_ensembleids(infileL)
	fw = open(outfile,'w')
	fw.write('gene names,MGI ids,Ensemble ids,predicted classes,probabilities (confidences) of predicted classes\n')
	for line in predictL[1:len(predictL)]:
		vals = line.split(',')
		genename = vals[0]
		pred_class = vals[1]
		prob = vals[2]
		if hash_genenames_mgiids.get(genename) != None:
			mgiid = hash_genenames_mgiids[genename]
		else:
			mgiid = 'none'
		if hash_genenames_ensembleids.get(genename) != None:
			ensembleid = hash_genenames_ensembleids[genename]
		else:
			ensembleid = 'none'
		line2 = genename+','+mgiid+','+ensembleid+','+pred_class+','+prob+'\n'
		fw.write(line2)
	fw.close()

def get_genenames_mgiids_ensembleids(infileL):
	hash_genenames_mgiids = {}
	hash_genenames_ensembleids = {}
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		ensembleid = vals[3]
		hash_genenames_mgiids[genename] = mgiid
		hash_genenames_ensembleids[genename] = ensembleid
	return (hash_genenames_mgiids,hash_genenames_ensembleids)

'''
def add_genenames_mgiids_ensembleids(hash_ids,predictL,outfile):
	# inst#     actual  predicted error prediction
	fw = open(outfile,'w')
	i=0
	while i < len(predictL):
		line = predictL[i]
		#if i == 4:
		#	print(line)
		m = re.match('^(inst#\s+actual\s+predicted\s+error\s+prediction)$',line)
		if m:
			#fw.write('gene name  '+m.group(1)+'\n')
			fw.write('GeneName\tMGI_ID\tEnsemble_ID\tinst#\tactual\tpredicted\terror\tprediction\n')			
			#print(line)
			break
		else:
			fw.write(line+'\n')
			#print(line+'\n')
			i +=1
	while i < len(predictL):
		line = predictL[i]
		m2 = re.match('^([\d]+).+$',line)
		m2 = re.match('^([\d]+)\s+(\d+\:\?)\s+(\d+\:[LethalVib]+)\s+([.\d]+)$',line)#matche line: '1        1:?   2:Viable       0.617'
		if m2:
			instance_no = int(m2.group(1))
			actual_class = m2.group(2)
			predicted_class = m2.group(3)
			confidence = m2.group(4)
			(genename,mgiid,ensembleid) = hash_ids[instance_no]
			fw.write(genename+'\t'+mgiid+'\t'+ensembleid+'\t'+str(instance_no)+'\t'+actual_class+'\t'+predicted_class+'\t \t'+confidence+'\n')
			i += 1
		else:
			i += 1
	fw.close()
'''

if __name__ == "__main__":
	main()



#Use the predictions of the mouse genes as the predictions of the human genes
#
#format of prediction file
#
#=== Predictions on test data ===
#
#gene names  inst#     actual  predicted error prediction
#Plekhg2  1        1:?   1:Lethal       0.557
#1700006E09Rik  2        1:?   2:Viable       0.826
#Cers4  3        1:?   2:Viable       0.622
#Cers5  4        1:?   1:Lethal       0.587
#Plekhg6  5        1:?   2:Viable       0.639
#Plekhg4  6        1:?   2:Viable       0.661
#
#orthologs:
#human gene	mouse gene
#ZNF721		Zfp721		
#ZNF140		Zfp140
#SLC5A4		Slc5a4a
#OR51G1		Olfr578
#C4orf21	Zgrf1
#ABCB1		Abcb1a
#C19orf45 	1700019B03Rik
#KIAA1462	9430020K01Rik
#ZNF793		Zfp793
#BOD1L1		Bod1l
#CES3		Ces3b
#OR6B3		Olfr1414

import sys
import re

def main():
	hash_predictions = {} #key = gene name, value = (predicted class, confidence of prediction)
	hash_genes_predicted_wrong = {} #key = gene name, value = (predicted_class,confidence,human_gene_actual_class,mouse_gene_actual_class,actual_classes_mismatch)
	n = 0 #no. of correct predictions of human genes
	t = 0 #no. of human genes with the same names as the mouse genes in prediction_file
	human_gene_actual_class = 'Lethal'
	human_genes_in_prediction_file = set()
	human_genes_not_in_prediction_file = set()
	orthologs_in_prediction_file = 0 #no. of orthologs in the prediction file
	orthologs = 0 #no. of orthologs which are correctly predicted
	orthologsL = []
	hash_orthologs = {}
	hash_known_genes = {}
	'''#human viable genes files 
	genesfile = 'c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes_not_in_train_set_and_test_set.csv'
	genesfile2 = 'c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes2_not_in_train_set_and_test_set.csv'
	testfile = 'c:/Users/David/Dropbox/datasets/essential genes prediction/test set/103 features data/new_lethal_new_viable_genes_not_in_train_set_missing_vals.csv'
	hash_orthologs = {}
	orthologsL = [('ZNF721','Zfp721'),('ZNF140','Zfp140'),('SLC5A4','Slc5a4a'),('OR51G1','Olfr578'),('C4orf21','Zgrf1'),('ABCB1','Abcb1a'),('C19orf45','1700019B03Rik'),('KIAA1462','9430020K01Rik'),('ZNF793','Zfp793'),('BOD1L1','Bod1l'),('CES3','Ces3b'),('OR6B3','Olfr1414')]
	genesL = [line.strip() for line in open(genesfile)]

	if len(genesL)==0:
		print(genesfile+" is empty.")
		sys.exit(-1)

	genes2L = [line.strip() for line in open(genesfile2)]

	if len(genes2L)==0:
		print(genesfile2+" is empty.")
		sys.exit(-1)

	testfileL = [line.strip() for line in open(testfile)]

	if len(testfileL)==0:
		print(testfile+" is empty.")
		sys.exit(-1)
	testfile_genenames = set()
	for line in testfileL[1:len(testfileL)]:
		vals = line.split(',')
		genename = vals[0]
		testfile_genenames.add(genename.lower())
	'''
	#human essential genes file
	outfile = sys.argv[1]
	essentialgenesfile = '/home/david/Dropbox/datasets/essential genes prediction/human essential genes/human essential genes.csv' #from paper "Analysis of protein-coding genetic 1 variation in 60,706 humans"
	prediction_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/random forest with AUC 0_86/prediction_reduced_data_with_non-mouse_genes_removed_and_duplicates_removed.random_forest_230_trees'
	knowngenesfile = '/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'
	genesL = [line.strip() for line in open(essentialgenesfile)]
	if len(genesL)==0:
		print(genesfile+" is empty.")
		sys.exit(-1)
	predL = [line.strip() for line in open(prediction_file)]
	if len(predL)==0:
		print(prediction_file+" is empty.")
		sys.exit(-1)
	knowngenesL = [line.strip() for line in open(knowngenesfile)]
	###store essentialites of genes in the known genes set into hash table
	for line in knowngenesL[1:len(knowngenesL)]:
		vals = line.split(',')
		genename = vals[0]
		ensembleid = vals[3]
		essentiality = vals[len(vals)-1]
		hash_known_genes[genename.lower()] = (essentiality,ensembleid)	
	###store predictions into hash table
	k = 0 #no. of genes with more than predictions
	for line in predL:
		prediction_obtained = False
		#match line: Plekhg2	1:Lethal	0.557
		m = re.match('^([\(\)\?.\-\w\s]+)\s+\d+\:([LethalVib]+)\s+([.\d]+)$',line)
		if m:
			genename = m.group(1)
			prediction = m.group(2)
			confidence = m.group(3)
			prediction_obtained = True
		else:   
			prediction_obtained = False
		if prediction_obtained == True:
			if hash_predictions.get(genename) == None:
				hash_predictions[genename.lower()] = (prediction,confidence)
			else:
				k += 1
				print('There are duplicate predictions of '+genename+'. Its first prediction is chosen.') 
				print('The first prediction: '+str(hash_predictions[genename.lower()]))
				print('Duplicate prediction: '+'('+prediction+','+confidence+')')
	print(str(k)+' genes have numerous predictions in prediction_file.')
	print('no. of predictions: '+str(len(hash_predictions)))
	###read orthologs
	'''
	for pair in orthologsL:
		human_gene = pair[0].lower()
		mouse_gene = pair[1].lower()
		hash_orthologs[human_gene] = mouse_gene
	#print(hash_orthologs)
	#get the orthologs which are in the prediction file but not in test set and training set
	#if a gene name is in the prediction file and it is not in the training set, it may be in the test set
	ks = list(hash_orthologs.keys())
	prediction_file_genenames = list(hash_predictions.keys())
	for k in ks:
		if hash_orthologs[k] in testfile_genenames: #delete the orthologs in the test set
			del hash_orthologs[k]
		if hash_orthologs[k] not in prediction_file_genenames:#delete the orthologs in the training set
			del hash_orthologs[k]
	print('orthologs in prediction file but not in the training set and the test set: '+str(len(hash_orthologs)))
	
	'''
	###get predictions of genes and orthologs
	(t,n,genes_of_different_essentialities,human_genes_in_known_gene_set) = get_predictions(genesL,hash_known_genes,hash_orthologs,hash_predictions,human_gene_actual_class,hash_genes_predicted_wrong,human_genes_in_prediction_file,human_genes_not_in_prediction_file,orthologs_in_prediction_file,t,n)
	print('total no. of human essential genes: '+str(len(genesL)))
	print('t: '+str(t)+' n: '+str(n))
	#get_predictions(genesL2,hash_known_genes,hash_orthologs,hash_predictions,human_gene_actual_class,hash_genes_predicted_wrong,human_genes_in_prediction_file,human_genes_not_in_prediction_file,orthologs_in_prediction_file,n,t)
	print(str(len(human_genes_in_prediction_file))+' human '+human_gene_actual_class+' genes have the same names as the mouse genes in the prediction_file.')
	print(str(len(human_genes_not_in_prediction_file))+' human '+human_gene_actual_class+' genes have different names to the mouse genes in the prediction_file.')
	print(str(len(human_genes_in_known_gene_set))+' human '+human_gene_actual_class+' genes are in known essentiality genes set.')
	print(str(len(genes_of_different_essentialities))+' '+human_gene_actual_class+' human genes have different essentiality to known essentiality gene set.')
	#write the human essential genes which are viable mouse genes in a file
	fw = open(outfile,'w')
	l = list(genes_of_different_essentialities)
	for genename_ensembleid in l:
		genename = genename_ensembleid[0]
		ensembleid = genename_ensembleid[1]
		fw.write(genename+','+ensembleid+'\n')
	fw.close()
	###print results of human genes prediction
	print('Using the predicted classes of mouse genes as the predicted classes of human genes having the same names as the mouse genes')
	print('no. of correct predictions: '+str(n))
	print('prediction accuracy on '+str(t)+' human '+human_gene_actual_class+' genes: '+str(float(n)/float(t)))
	if orthologsL != []:
		print('There are '+str(len(orthologsL))+' pairs of orthologs in total.')
		print('The prediction file contains '+str(orthologs_in_prediction_file)+' orthologs.')
		print(str(orthologs)+' orthologs of '+str(orthologs_in_prediction_file)+' orthologs in the prediction file are predicted correctly.')
	
def get_predictions(genesL,hash_known_genes,hash_orthologs,hash_predictions,human_gene_actual_class,hash_genes_predicted_wrong,human_genes_in_prediction_file,human_genes_not_in_prediction_file,orthologs_in_prediction_file,t,n):
	human_genes_in_known_gene_set = set()
	#essentiality_mismatch = 0
	genes_of_different_essentialities = set()
	for gene in genesL:
		genename = gene
		#print(genename)
		if hash_orthologs != {}:
			if hash_orthologs.get(genename.lower()) != None:
				mouse_gene = hash_orthologs[genename.lower()]		
				genename = mouse_gene.lower()
				is_ortholog = True
			else:
				is_ortholog = False
		else:
			is_ortholog = 'no_orthologs'
		if hash_predictions.get(genename.lower()) != None:
			t += 1
			if is_ortholog == True:
				orthologs_in_prediction_file += 1
			(predicted_class,confidence) = hash_predictions[genename.lower()]
			if predicted_class == human_gene_actual_class:
				n += 1
				if is_ortholog == True:
					orthologs += 1	
			else:
				hash_genes_predicted_wrong[genename.lower()] = (predicted_class,confidence,human_gene_actual_class)
			human_genes_in_prediction_file.add(genename)
		else:
			human_genes_not_in_prediction_file.add(genename)
			if hash_known_genes.get(genename.lower()) != None:
				human_genes_in_known_gene_set.add(genename)
				(essentiality,ensembleid) = hash_known_genes.get(genename.lower())
				if essentiality != 'Lethal':
					#essentiality_mismatch += 1
					genes_of_different_essentialities.add((genename,ensembleid))
	return (t,n,genes_of_different_essentialities,human_genes_in_known_gene_set)

if __name__ == "__main__":
	main()	


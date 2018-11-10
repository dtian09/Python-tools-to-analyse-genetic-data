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
	genesfile = 'c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes_not_in_train_set_and_test_set.csv'
	genesfile2 = 'c:/Users/David/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes2_not_in_train_set_and_test_set.csv'	
	testfile = 'c:/Users/David/Dropbox/datasets/essential genes prediction/test set/103 features data/new_lethal_new_viable_genes_not_in_train_set_missing_vals.csv'
	prediction_file = 'c:/Users/David/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/random forest with AUC 0_86/prediction_reduced_data_with_ids.random_forest_230_trees'
	
	hash_predictions = {} #key = gene name, value = (predicted class, confidence of prediction)
	hash_mouse_gene_actual_classes = {} #key = gene name, value = actual class of the mouse gene with this gene name
	hash_genes_predicted_wrong = {} #key = gene name, value = (predicted_class,confidence,human_gene_actual_class,mouse_gene_actual_class,actual_classes_mismatch)
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

	predL = [line.strip() for line in open(prediction_file)]

	if len(predL)==0:
		print(prediction_file+" is empty.")
		sys.exit(-1)

	testfile_genenames = set()
	for line in testfileL[1:len(testfileL)]:
		vals = line.split(',')
		genename = vals[0]
		testfile_genenames.add(genename.lower())
	###read data (actual classes and predictions) into hash tables
	k = 0 #no. of genes with more than predictions
	for line in predL:
		prediction_obtained = False
		#m = re.match('^([\(\)\?.\-\w\s]+)\s+([MGI\?\:\d]+)\s+([ENSMUSG\?\d]+)\s+\d+\s+\d+\:\?\s+\d+\:([LethalVib]+)\s+([.\d]+)$',line)#match line: 'Bckdk	MGI:1276121	ENSMUSG00000030802	16	1:?	1:Lethal	 	0.613'	
		m = re.match('^([\(\)\?.\-\w\s]+)\s+([MGI\?\:\d]+)\s+([ENSMUSG\?\d]+)\s+\d+\s+(\d+\:[^\s]*)\s+\d+\:([LethalVib]+)\s+([.\d]+)$',line)#match line: 'Bckdk	MGI:1276121	ENSMUSG00000030802	16	1:?	1:Lethal	 	0.613'	
		if m:
			genename = m.group(1)
			#mgi_id = m.group(2)
			#ensemble_id = m.group(3)
			mouse_gene_actual_class = m.group(4)
			prediction = m.group(5)
			confidence = m.group(6)
			hash_mouse_gene_actual_classes[genename.lower()] = mouse_gene_actual_class
			prediction_obtained = True
		else:   
			prediction_obtained = False
		#	print('This line does not match patterns: '+line)
		if prediction_obtained == True:
			if hash_predictions.get(genename) == None:
				hash_predictions[genename.lower()] = (prediction,confidence)
			else:
				k += 1
				print('There are duplicate predictions of '+genename+'. Its first prediction is chosen.') 
				print('The first prediction: '+str(hash_predictions[genename.lower()]))
				print('Duplicate prediction: '+'('+prediction+','+confidence+')')
		#else:
		#	print('Prediction is not obtained from this line: '+line)
	print(str(k)+' genes have numerous predictions in prediction_file.')
	###read orthologs
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
	n = 0 #no. of correct predictions of human genes
	t = 0 #no. of human genes with the same names as the mouse genes in prediction_file
	human_gene_actual_class = 'Viable'
	human_genes_in_prediction_file = set()
	human_genes_not_in_prediction_file = set()
	orthologs_in_prediction_file = 0 #no. of orthologs in the prediction file
	orthologs = 0 #no. of orthologs which are correctly predicted
	###get predictions of genes and orthologs
	for gene in genesL:
		#valsL = gene.split(',')
		#genename = valsL[0]
		#mgi_id = valsL[2]
		#ensemble_id = valsL[3]
		genename = gene
		if hash_orthologs.get(genename.lower()) != None:
			mouse_gene = hash_orthologs[genename.lower()]		
			genename = mouse_gene.lower()
			is_ortholog = True
		else:
			is_ortholog = False
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
				mouse_gene_actual_class = hash_mouse_gene_actual_classes[genename.lower()]
				if human_gene_actual_class != mouse_gene_actual_class:
					actual_classes_mismatch = '1'
				else:
					actual_classes_mismatch = '0'
				hash_genes_predicted_wrong[genename.lower()] = (predicted_class,confidence,human_gene_actual_class,mouse_gene_actual_class,actual_classes_mismatch)
			human_genes_in_prediction_file.add(genename)
		else:
			human_genes_not_in_prediction_file.add(genename)
	for gene in genes2L:
		#valsL = gene.split('\t')
		#genename = valsL[0]
		#mgi_id = valsL[2]
		#ensemble_id = valsL[3]
		genename = gene
		if hash_orthologs.get(genename.lower()) != None:
			mouse_gene = hash_orthologs[genename.lower()]		
			genename = mouse_gene.lower()
			is_ortholog = True
		else:
			is_ortholog = False
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
				mouse_gene_actual_class = hash_mouse_gene_actual_classes[genename.lower()]
				if human_gene_actual_class != mouse_gene_actual_class:
					actual_classes_mismatch = '1'
				else:
					actual_classes_mismatch = '0'
				hash_genes_predicted_wrong[genename.lower()] = (predicted_class,confidence,human_gene_actual_class,mouse_gene_actual_class,actual_classes_mismatch)
			human_genes_in_prediction_file.add(genename)
		else:
			human_genes_not_in_prediction_file.add(genename)
	print(str(len(human_genes_in_prediction_file))+' human viable genes have the same names as the mouse genes in the prediction_file.')
	'''
	for genename in list(human_genes_not_in_prediction_file):
		print(genename)
	'''
	print(str(len(human_genes_not_in_prediction_file))+' human viable genes have different names to the mouse genes in the prediction_file.')
	'''
	for genename in list(human_genes_not_in_prediction_file):
		print(genename)
	'''
	###print results of human viable genes prediction
	print('Using the predicted classes of mouse genes as the predicted classes of human genes having the same names as the mouse genes')
	print('no. of correct predictions: '+str(n))
	print('prediction accuracy on '+str(t)+' human viable genes: '+str(float(n)/float(t)))
	print('There are '+str(len(orthologsL))+' pairs of orthologs in total.')
	print('The prediction file contains '+str(orthologs_in_prediction_file)+' orthologs.')
	print(str(orthologs)+' orthologs of '+str(orthologs_in_prediction_file)+' orthologs in the prediction file are predicted correctly.')
	print('human genes predicted wrong:')
	print('GeneName\tPredicted Class\tConfidence\tHuman Gene Actual Class\tMouse Gene Actual Class\tActual Classes Mismatch\n')
	ks = list(hash_genes_predicted_wrong.keys())
	for k in ks:
		(predicted_class,confidence,human_gene_actual_class,mouse_gene_actual_class,actual_classes_mismatch) = hash_genes_predicted_wrong[k]
		print(k+'\t'+predicted_class+'\t'+confidence+'\t'+human_gene_actual_class+'\t'+mouse_gene_actual_class+'\t'+actual_classes_mismatch+'\n')
	###print results of human lethal (essential) genes prediction
	
if __name__ == "__main__":
	main()	


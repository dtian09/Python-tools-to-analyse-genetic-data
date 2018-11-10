#get the lethal, viable, predicted lethal, predicted viable genes of each chromosome
#file format of infile, infile2, infile3: 
#
#MGI symbol,Ensembl Gene ID,Chromosome Name
#2700049A03Rik,ENSMUSG00000034601,12
#A3galt2,ENSMUSG00000028794,4
#A4galt,ENSMUSG00000047878,15
#A4gnt,ENSMUSG00000037953,9

import sys
import re

def main():
	infile = '/home/david/Dropbox/working papers/gene essentiality predictions/percentages of lethal and viable genes/known_genes_chromosomes.csv'
	infile2 = '/home/david/Dropbox/working papers/gene essentiality predictions/percentages of lethal and viable genes/new_lethal_new_viable_genes_chromosomes.csv'
	infile3 = '/home/david/Dropbox/working papers/gene essentiality predictions/percentages of lethal and viable genes/unknown_genes_chromosomes.csv'
	infile4 = '/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'#infile4 contains labels of genes
	infile5 = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/new_lethal_new_viable_genes_not_in_train_set_missing_vals.csv'#infile5 contains labels of genes
	infile6 =  '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_48_features_data_with_gene_names.random_forest.space'#infile6 contains the predicted labels of genes
	infileL = [line.strip() for line in open(infile)]
	infileL2 = [line.strip() for line in open(infile2)]
	infileL3 = [line.strip() for line in open(infile3)]
	infileL4 = [line.strip() for line in open(infile4)]
	infileL5 = [line.strip() for line in open(infile5)]
	infileL6 = [line.strip() for line in open(infile6)]
	
	hash_genes = {} #key=chromosome name, value = (set of known lethal genes of chromosome,set of known viable genes of chromosome,set of predicted lethal genes of chromosome,set of predicted viable genes of chromosome)
	hash_genes = get_genes(hash_genes,infileL,infileL2,infileL3,infileL4,infileL5,infileL6)
	#percentages of lethal, viable, predicted lethal, predicted viable of each chromosome
	print_in_table_format(hash_genes)
	'''
	print('no. of chromosomes: '+str(len(hash_genes)))
	chromosomes = list(hash_genes.keys())
	for chromosome in chromosomes:
		(lethal,viable,predicted_lethal,predicted_viable) = hash_genes[chromosome]
		
		#print(chromosome)
		#print('lethal: '+str(lethal)+'\n')
		#print('viable: '+str(viable)+'\n')
		#print('predicted lethal: '+str(predicted_lethal)+'\n')
		#print('predicted viable: '+str(predicted_viable)+'\n')
		
		n_lethal = len(lethal)
		n_viable = len(viable)
		n_pred_lethal = len(predicted_lethal)
		n_pred_viable = len(predicted_viable)
		n = n_lethal + n_viable + n_pred_lethal + n_pred_viable
		percent_lethal = float(n_lethal) / float(n) * 100
		percent_viable = float(n_viable) / float(n) * 100
		percent_pred_lethal = float(n_pred_lethal) / float(n) * 100
		percent_pred_viable = float(n_pred_viable) / float(n) * 100		
		print(chromosome+': total no. of genes:'+str(n)+', lethal:'+str(len(lethal))+', viable:'+str(len(viable))+', predicted lethal:'+str(len(predicted_lethal))+', predicted viable:'+str(len(predicted_viable))+'\n')
		print(chromosome+': lethal:'+str(percent_lethal)+'%, viable:'+str(percent_viable)+'%, predicted lethal:'+str(percent_pred_lethal)+'%, predicted viable:'+str(percent_pred_viable)+'%\n')
	'''	

def print_in_table_format(hash_genes):
	i = 0
	chromosomes = list(hash_genes.keys())
	for chromosome in chromosomes:
		if i == 0:
			print(','+chromosome+','),
			i += 1
		else:
			print(chromosome+','),
	print('')
	print('lethal,'),
	for chromosome in chromosomes:
		(lethal,viable,predicted_lethal,predicted_viable) = hash_genes[chromosome]
		n_lethal = len(lethal)
		n_viable = len(viable)
		n_pred_lethal = len(predicted_lethal)
		n_pred_viable = len(predicted_viable)
		n = n_lethal + n_viable + n_pred_lethal + n_pred_viable
		percent_lethal = float(n_lethal) / float(n) * 100
		print(str(percent_lethal)+','),
	print('')
	print('viable,'),
	for chromosome in chromosomes:
		(lethal,viable,predicted_lethal,predicted_viable) = hash_genes[chromosome]
		n_lethal = len(lethal)
		n_viable = len(viable)
		n_pred_lethal = len(predicted_lethal)
		n_pred_viable = len(predicted_viable)
		n = n_lethal + n_viable + n_pred_lethal + n_pred_viable
		percent_viable = float(n_viable) / float(n) * 100
		print(str(percent_viable)+','),
	print('')
	print('predicted lethal,'),
	for chromosome in chromosomes:
		(lethal,viable,predicted_lethal,predicted_viable) = hash_genes[chromosome]
		n_lethal = len(lethal)
		n_viable = len(viable)
		n_pred_lethal = len(predicted_lethal)
		n_pred_viable = len(predicted_viable)
		n = n_lethal + n_viable + n_pred_lethal + n_pred_viable
		percent_pred_lethal = float(n_pred_lethal) / float(n) * 100
		print(str(percent_pred_lethal)+','),
	print('')
	print('predicted viable,'),
	for chromosome in chromosomes:
		(lethal,viable,predicted_lethal,predicted_viable) = hash_genes[chromosome]
		n_lethal = len(lethal)
		n_viable = len(viable)
		n_pred_lethal = len(predicted_lethal)
		n_pred_viable = len(predicted_viable)
		n = n_lethal + n_viable + n_pred_lethal + n_pred_viable
		percent_pred_viable = float(n_pred_viable) / float(n) * 100
		print(str(percent_pred_viable)+','),
	print('')

def get_genes(hash_genes,infileL,infileL2,infileL3,infileL4,infileL5,infileL6):
	hash_genes1 = {}#key=chromosome name, value = set of genes of the chromosome
	hash_genes2 = {} #key=chromosome name, value = (set of known lethal genes of chromosome,set of known viable genes of chromosome,set of predicted lethal genes of chromosome,set of predicted viable genes of chromosome)
	#get all the chromosomes and the genes of each chromosome
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		chromosome = vals[2]
		if hash_genes1.get(chromosome)!=None:
			genes = hash_genes1[chromosome]
			genes.add(genename)
			hash_genes1[chromosome] = genes
		else:
			genes = set()
			genes.add(genename)
			hash_genes1[chromosome] = genes
	for line in infileL2[1:len(infileL2)]:
		vals = line.split(',')
		genename = vals[0]
		chromosome = vals[2]
		if hash_genes1.get(chromosome)!=None:
			genes = hash_genes1[chromosome]
			genes.add(genename)
			hash_genes1[chromosome] = genes
		else:
			genes = set()
			genes.add(genename)
			hash_genes1[chromosome] = genes
	for line in infileL3[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		chromosome = vals[2]
		if hash_genes1.get(chromosome)!=None:
			genes = hash_genes1[chromosome]
			genes.add(genename)
			hash_genes1[chromosome] = genes
		else:
			genes = set()
			genes.add(genename)
			hash_genes1[chromosome] = genes
	'''
	chromosomesL = list(hash_genes1.keys())
	for chromosome in chromosomesL:
		print(chromosome+': '+str(hash_genes1[chromosome])+'\n')
	'''
	#get the labels of the genes
	lethal = set()
	viable = set()
	predicted_lethal = set()
	predicted_viable = set()
	for line in infileL4[1:len(infileL4)]:
		vals = line.split(',')	
		genename = vals[0]
		label = vals[len(vals)-1]
		if label == 'Lethal':
			lethal.add(genename)
		else:
			viable.add(genename)
	for line in infileL5[1:len(infileL5)]:
		vals = line.split(',')	
		genename = vals[0]
		label = vals[len(vals)-1]
		if label == 'Lethal':
			lethal.add(genename)
		else:
			viable.add(genename)
	for line in infileL6[1:len(infileL6)]:
		vals = line.split(',')
		genename = vals[0]
		label = vals[3]
		if label == '1:Lethal':
			predicted_lethal.add(genename)
		else:
			predicted_viable.add(genename)
	#print('lethal genes: '+str(lethal))
	#print('viable genes: '+str(viable))
	#print('predicted lethal genes: '+str(predicted_lethal))
	#get the lethal, viable, predicted lethal, predicted viable genes of each chromosome
	chromosomesL = list(hash_genes1.keys())
	for chromosome in chromosomesL:
		lethal_of_chromosome = set()
		viable_of_chromosome = set()
		predicted_lethal_of_chromosome = set()
		predicted_viable_of_chromosome = set()
		genes = hash_genes1[chromosome]
		genesL = list(genes)
		for gene in genesL:
			if gene in lethal:
				lethal_of_chromosome.add(gene)
			elif gene in viable:
				viable_of_chromosome.add(gene)
			elif gene in predicted_lethal:
				predicted_lethal_of_chromosome.add(gene)
			else:
				predicted_viable_of_chromosome.add(gene)
		hash_genes2[chromosome] = (lethal_of_chromosome,viable_of_chromosome,predicted_lethal_of_chromosome,predicted_viable_of_chromosome)
	return hash_genes2	
				
if __name__ == "__main__":
	main()


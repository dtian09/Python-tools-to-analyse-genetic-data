#get prediction of newly known lethal genes
#input: newly known lethal genes csv file
#	known genes dataset csv file
#	classifier prediction file
#output: prediction of newly known lethal genes file
#
#output file format:
#
#gene name,in training file, prediction
#g1, y, lethal
#g2, n, viable
#...
#percentage of correct prediction

import sys
import re

def main():
        '''
        infile = './unknown essentiality genes/newly known lethal genes/newly_known_lethal_genes.csv'
        infile2 = './known essentiality genes/training_all_genesinfo.csv'
        infile3 = './unknown essentiality genes/predictions/prediction_with_gene_names.chimerge_significance_level0.95_random_forest_ROC_area=0.998'
        outfile =  './unknown essentiality genes/newly known lethal genes/prediction_of_lethal_genes.chimerge_significance_level0.95_random_forest_ROC_area=0.998'
        '''
        #infile3 = './unknown essentiality genes/predictions/prediction_with_gene_names.chimerge_significance_level0.95_c45_ROC_area=0.962'
        #outfile =  './unknown essentiality genes/newly known lethal genes/prediction_of_lethal_genes.chimerge_discretized_sig_level0.95_c45_ROC_area=0.962'
        #infile3 = './unknown essentiality genes/predictions/prediction_with_gene_names.chimerge_significance_level0.95_boolean_transformed_linear_svm_ROC_area=0.949'
        #outfile =  './unknown essentiality genes/newly known lethal genes/prediction_of_lethal_genes.chimerge_significance_level0.95_boolean_transformed_linear_svm_ROC_area=0.949'
        infile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known lethal genes\\newly_known_lethal_genes.csv'
        infile2 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\known essentiality genes\\training_all_genesinfo_corrected_classes.csv'
        infile3 = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\predictions\\linear_svm_normalized_predictions_with_gene_names'
        outfile = 'C:\\Users\\David\\Dropbox\\datasets\\essential genes prediction\\unknown essentiality genes\\newly known lethal genes\\prediction_of_lethal_genes.linear_svm_normalized'

        infileL = [line.strip() for line in open(infile)]

        if len(infileL)==0:
                print(infile+" is empty.")
                sys.exit(-1)

        infile2L = [line.strip() for line in open(infile2)]

        if len(infile2L)==0:
                print(infile2+" is empty.")
                sys.exit(-1)

        infile3L = [line.strip() for line in open(infile3)]

        if len(infile3L)==0:
                print(infile3+" is empty.")
                sys.exit(-1)

        lethal_genes = set()
        for gene in infileL[1:len(infileL)]:
                lethal_genes.add(gene)

        training_set_genes = {}
        
        for line in infile2L[1:len(infile2L)]:
                vals = line.split(',')		
                gene_name = vals[0]
                gene_class = vals[len(vals)-1] 
                training_set_genes[gene_name] = gene_class

        hash_predictions = {}#key=gene name, value=prediction
        for line in infile3L:
        #Plekhg2  1        1:?   2:Viable       0.617
                m = re.match('^.+\d+:(Viable).+$',line)
                m2 = re.match('^.+\d+:(Lethal).+$',line)
                if m:
                        valsL = line.split(' ')
                        gene_name = valsL[0]
                        hash_predictions[gene_name] = m.group(1)
                elif m2:
                        valsL = line.split(' ')
                        gene_name = valsL[0]
                        hash_predictions[gene_name] = m2.group(1)
        write_to_file(lethal_genes,training_set_genes,hash_predictions,outfile)

def write_to_file(lethal_genes,training_set_genes,hash_predictions,outfile):
	lethal = 0
	genes_not_in_training_set = 0 # count lethal genes not in training set
	n=0
	wrong_viable_genes_in_training_set=''
	fw = open(outfile,"w")
	fw.write("gene names,in training set,prediction\n")
	ks = list(training_set_genes.keys())
	for gene_name in list(lethal_genes):
		if gene_name in ks:
			in_training_set = 'y'
			if training_set_genes[gene_name] == 'Viable':
				print(gene_name+" is Viable in training set. It is Lethal in newly known lethal gene set.\n")
				n += 1
				wrong_viable_genes_in_training_set += gene_name+","
		else:
			in_training_set = 'n'
			genes_not_in_training_set += 1
		if hash_predictions.get(gene_name)!= None:
			prediction = hash_predictions[gene_name]
			if prediction == 'Lethal':
				lethal += 1
		else:
			#print("hash_predictions has no "+gene_name+" in_training_set: "+in_training_set+"\n")
			prediction = '?'
			#lethal += 1#count the missing prediction as Lethal class prediction
		fw.write(gene_name+','+in_training_set+","+prediction+"\n")
	
	ratio = float(lethal)/float(genes_not_in_training_set)
	percent = str(ratio*100)
	print("ratio: "+str(ratio)+", percentage: "+percent+"\n")
	print(str(n)+" genes are lethal in newly known gene set, but they are viable in training set.\n")
	print(wrong_viable_genes_in_training_set)
	print("no. of genes in training set: "+str(len(lethal_genes)-genes_not_in_training_set)+"\n")
	fw.write("Percentage of correct prediction of lethal genes: "+percent+"%\n")
	fw.close()

if __name__ == "__main__":
	main()


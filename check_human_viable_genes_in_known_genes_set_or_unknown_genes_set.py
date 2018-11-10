#Check whether the human viable genes are in the known essentiality mouse gene set or in the unknown essentiality mouse gene set.

import sys

def main():	
	infile='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'
	infile2 ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/merged_gene_protein_features_with_ids_missing_vals_non-mouse_genes_removed_and_duplicates_removed.csv'
	infile3='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes.csv'
	infile4='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes2.csv'
	outfile='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/all_human_viable_genes.csv'
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
	
	infile4L = [line.strip() for line in open(infile4)]
	if len(infile4L)==0:
		print(infile4+" is empty.")
		sys.exit(-1)
	
	known_mouse_genes = set()
	unknown_mouse_genes = set()
	human_genes = set()
	original_gene_names = {}
	human_genes_in_known_mouse_genes = set()
	human_genes_in_unknown_mouse_genes = set()
	human_genes_not_in_known_and_unknown_mouse_genes = set()
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		known_mouse_genes.add(genename)
	for line in infile2L[1:len(infile2L)]:
		valsL = line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		#if genename not in known_mouse_genes:
		unknown_mouse_genes.add(genename)
	for line in infile3L:
		original_genename = line
		genename = original_genename.lower()
		human_genes.add(genename)
		original_gene_names[genename] = original_genename
		if genename in known_mouse_genes:
 			human_genes_in_known_mouse_genes.add(genename)
		elif genename in unknown_mouse_genes:
			human_genes_in_unknown_mouse_genes.add(genename)
		else:
			human_genes_not_in_known_and_unknown_mouse_genes.add(genename)
	for line in infile4L:
		original_genename = line
		genename = original_genename.lower()
		human_genes.add(genename)
		original_gene_names[genename] = original_genename
		if genename in known_mouse_genes:
 			human_genes_in_known_mouse_genes.add(genename)
		elif genename in unknown_mouse_genes:
			human_genes_in_unknown_mouse_genes.add(genename)
		else:
			human_genes_not_in_known_and_unknown_mouse_genes.add(genename)
	'''
	orthologsL = [('ZNF721','Zfp721'),('ZNF140','Zfp140'),('SLC5A4','Slc5a4a'),('OR51G1','Olfr578'),('C4orf21','Zgrf1'),('ABCB1','Abcb1a'),('C19orf45','1700019B03Rik'),('KIAA1462','9430020K01Rik'),('ZNF793','Zfp793'),('BOD1L1','Bod1l'),('CES3','Ces3b'),('OR6B3','Olfr1414')]
	hash_orthologs = {}#key= human gene, value = mouse ortholog of human gene
	human_genes_whose_orthologs_in_known_mouse_genes = set()
	human_genes_whose_orthologs_in_unknown_mouse_genes = set()
	for pair in orthologsL:
		human_gene = pair[0]
		mouse_gene = pair[1]
		hash_orthologs[human_gene] = mouse_gene
	for human_gene in list(human_genes_not_in_known_and_unknown_mouse_genes):
		original_human_gene = original_gene_names[human_gene]
		if hash_orthologs.get(original_human_gene)!= None:
			mouse_ortholog = hash_orthologs[original_human_gene]
			if mouse_ortholog.lower() in known_mouse_genes:
				human_genes_whose_orthologs_in_known_mouse_genes.add(original_human_gene)
				human_genes_not_in_known_and_unknown_mouse_genes.remove(human_gene)#if the ortholog of a human gene is in the known mouse gene set or unknown mouse gene set, the human gene is in the known mouse gene set or unknown mouse gene set.
			elif  mouse_ortholog.lower() in unknown_mouse_genes:
				human_genes_whose_orthologs_in_unknown_mouse_genes.add(original_human_gene)
				human_genes_not_in_known_and_unknown_mouse_genes.remove(human_gene)
	'''
	fw = open(outfile,'w')
	fw.write('GeneName,in Known Mouse Genes Set,in Unknown Mouse Genes Set\n')
	for human_gene in list(human_genes_in_known_mouse_genes):
		original_human_gene = original_gene_names[human_gene]
		fw.write(original_human_gene+','+'Yes,'+'No\n')
	for human_gene in list(human_genes_in_unknown_mouse_genes):
		original_human_gene = original_gene_names[human_gene]
		fw.write(original_human_gene+','+'No,'+'Yes\n')
	'''
	#print orthologs to file
	for human_gene in list(human_genes_whose_orthologs_in_known_mouse_genes):
		original_human_gene = original_gene_names[human_gene.lower()]
		ortholog = hash_orthologs[original_human_gene]
		fw.write(original_human_gene+','+'Yes(ortholog: '+ortholog+'),'+'No\n')
	for human_gene in list(human_genes_whose_orthologs_in_unknown_mouse_genes):
		original_human_gene = original_gene_names[human_gene.lower()]
		ortholog = hash_orthologs[original_human_gene]
		fw.write(original_human_gene+','+'No,'+'Yes(ortholog: '+ortholog+')\n')
	'''
	for human_gene in list(human_genes_not_in_known_and_unknown_mouse_genes):
		original_human_gene = original_gene_names[human_gene.lower()]
		fw.write(original_human_gene+',No,No\n')
	fw.close()
	
	print(str(len(known_mouse_genes))+' known mouse genes')
	print(str(len(unknown_mouse_genes))+' unknown mouse genes')
	print(str(len(human_genes))+' human genes')
	print(str(len(human_genes_in_known_mouse_genes))+' human genes in known mouse genes.')
	print(str(len(human_genes_in_unknown_mouse_genes))+' human genes in unknown mouse genes.')
	print(str(len(human_genes_not_in_known_and_unknown_mouse_genes))+' human genes are not in the known and the unknown mouse genes sets')
	#print(str(len(orthologsL))+' of the '+str(len(human_genes_not_in_known_and_unknown_mouse_genes))+' human genes have mouse orthologs.')
	#print('The set of the known and the unknown mouse genes contains the ortholog(s) of '+str(len(human_genes_whose_orthologs_in_known_mouse_genes)+len(human_genes_whose_orthologs_in_unknown_mouse_genes))+' human genes.')
	#print('The known mouse genes set contains the ortholog(s) of '+str(len(human_genes_whose_orthologs_in_known_mouse_genes))+' human genes.')
	#print('The unknown mouse genes set contains the ortholog(s) of '+str(len(human_genes_whose_orthologs_in_unknown_mouse_genes))+' human genes.')
	
if __name__ == "__main__":
        main()


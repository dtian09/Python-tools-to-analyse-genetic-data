#program to: count the number of genes which encode 1 protein only
#            count the number of genes which encode more than 1 proteins
# 	     count the number of proteins of all the genes
#input: an uniprot-mapped fasta file
#output: number of genes which have encode 1 protein only,
#        number of genes which have encode more than 1 proteins
# fasta sequence header format:
#>tr|A2AHY8|A2AHY8_MOUSE Pet2 protein OS=Mus musculus GN=Pet2 PE=2 SV=1
#
#a fasta file may have noise so that the gene name of a protein is not one of the genes to be collected. Another noise is that the gene name of a protein is missing. In this case, the gene name of the previous protein in the fasta file is the gene name of that protein.
#
#For example, Supt4a is to be collected, but Supt4h1a is not to be collected. In this case, Supt4h1a is Supt4a and the protein of Supt4h1a is a protein of Supt4a.
#
#>tr|G3UVU8|G3UVU8_MOUSE MCG7669, isoform CRA_a OS=Mus musculus GN=Supt4a PE=4 SV=1
#MALETVPKDLRHLRACLLCSLVKTIDQFEYDGCDNCDAYLQMKGNREMVYDCTSSSFDGC
#LARAVALGSCCGPLQQGRLGLTLCLPL
#>sp|P63271|SPT4A_MOUSE Transcription elongation factor SPT4-A OS=Mus musculus GN=Supt4h1a PE=2 SV=1
#MALETVPKDLRHLRACLLCSLVKTIDQFEYDGCDNCDAYLQMKGNREMVYDCTSSSFDGI
#IAMMSPEDSWVSKWQRVSNFKPGVYAVSVTGRLPQGIVRELKSRGVAYKSRDTAIKT
#
#
import sys
import re

def main():
	one_protein_genes=0 #number of genes encoding one protein
	numerous_proteins_genes=0    #number of genes encoding numerous proteins
	total_proteins=0    #number of proteins of all the genes	
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_ids.csv'#id of genes to collect from the fasta file
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes_ids.csv'
	#fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes.fasta2'
	#fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/uniprot-1700067K01Rik.fasta'
	#fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/examples.fasta'
	fasta_file = '/home/david/Dropbox/datasets/essential genes prediction/new viable genes/new_viable_genes.fasta'
	genesid_fileL = [line.strip() for line in open(genesid_file)]
	fasta_fileL = [line.strip() for line in open(fasta_file)]
	genenames = set()
	for line in genesid_fileL[1:len(genesid_fileL)]:
		valsL=line.split(',')
		original_genename = valsL[0]
		genename = original_genename.lower()
		genenames.add(genename)
	(proteins,uniprotids_entrynames,proteins_of_genes_not_to_collect,proteins_with_no_gene_names,genes_not_collected) = get_proteins(list(genenames),fasta_fileL)
	if len(genes_not_collected)>0:
		print('in genes.py: fasta_file does not contain protein information of these genes in genesid_file:')
		for gene in list(genes_not_collected):
			print(gene)
	genes = list(proteins.keys())
	for gene in genes:
		num_proteins = len(proteins[gene])#key=gene name, value= set of tuples (protein name,uniprot id, uniprot entry name)
		total_proteins += num_proteins
		if num_proteins==1:
			one_protein_genes+=1
		else:
			numerous_proteins_genes+=1
	print ("number of genes which encode 1 protein only: %d" %one_protein_genes)
	print ("number of genes which encode more than one proteins: %d" %numerous_proteins_genes)
	print ("number of genes whose proteins are collected: %d " %len(genes))
	print ("total number of different proteins of all the genes: %d " %total_proteins)
	print(str(len(proteins_of_genes_not_to_collect))+' proteins in the fasta file are encoded by genes not to be collected.')
	for protein in list(proteins_of_genes_not_to_collect):
		print(protein)
	print(str(len(proteins_with_no_gene_names))+' proteins have no gene names in the fasta file.')
	for protein in list(proteins_with_no_gene_names):
		print(protein)
		
def get_proteins(genes_to_collectL,fasta_fileL):
	#input: names of genes to collect from fasta file
	#	fasta file
	proteins = {} #a hashtable which only stores the proteins (in fasta_file) of those genes (in genes_id file) to collect. 
		      #key=gene name, value= set of tuples (protein, uniprot id, uniprot entry name)		      
	proteins_of_genes_not_to_collect = set()
	proteins_with_no_gene_names = set()
	genes_collected = set()
	uniprotids_entrynames = {}
	for line in fasta_fileL:
		(gene_protein_uniprotid_entry_name,proteins_of_genes_not_to_collect,proteins_with_no_gene_names) = get_protein_header(genes_to_collectL,genes_collected,proteins_of_genes_not_to_collect,proteins_with_no_gene_names,line)
		if gene_protein_uniprotid_entry_name !=None:
			gene_name = gene_protein_uniprotid_entry_name[0]
			gene_name_lower_case = gene_name.lower()
			protein_name = gene_protein_uniprotid_entry_name[1]
			uniprotid = gene_protein_uniprotid_entry_name[2]
			entry_name = gene_protein_uniprotid_entry_name[3]
			uniprotids_entrynames[uniprotid] = entry_name
			if proteins.get(gene_name_lower_case) != None:
				proteins_of_gene = proteins[gene_name_lower_case]
				proteins_of_gene.add((protein_name,uniprotid,entry_name))
				proteins[gene_name_lower_case] = proteins_of_gene
			else:
				proteins_of_gene = set()
				proteins_of_gene.add((protein_name,uniprotid,entry_name))
				proteins[gene_name_lower_case] = proteins_of_gene
	genes_not_collected = set(genes_to_collectL).difference(genes_collected)
	return (proteins,uniprotids_entrynames,proteins_of_genes_not_to_collect,proteins_with_no_gene_names,genes_not_collected)

def get_protein_header(genes_to_collectL,genes_collected,proteins_of_genes_not_to_collect,proteins_with_no_gene_names,line):
# fasta sequence header format:
#>tr|A2AHY8|A2AHY8_MOUSE Pet2 protein OS=Mus musculus GN=Pet2 PE=2 SV=1
	m = re.match('>.+\\|(.+)(\\|.+)', line)
	if m:
		uniprot_id = m.group(1)
		subStr = m.group(2)
		m2 = re.match('\\|([\w]+_[\w]+)\\s(.+)\\sOS=.+\\s+GN=(.+)\\s+PE.+',subStr)
		if m2 != None:
			entry_name = m2.group(1)
			protein_name = m2.group(2)
			gene_name = m2.group(3)
			gene_name = gene_name.lower()
			if gene_name not in genes_to_collectL:
				#print('Protein '+uniprot_id+' is encoded by gene '+gene_name+' which is not to be collected')
				proteins_of_genes_not_to_collect.add(uniprot_id)
			else:
				genes_collected.add(gene_name)
			return ([gene_name,protein_name,uniprot_id,entry_name],proteins_of_genes_not_to_collect,proteins_with_no_gene_names)
		else:
			#This protein has no gene name
			m3 = re.match('\\|([\w]+_[\w]+)\\s(.+)\\sOS=.+\\s+PE.+',subStr)#e.g. >sp|Q8C669|PELI1_MOUSE E3 ubiquitin-protein ligase pellino homolog 1 OS=Mus musculus PE=1 SV=2
			if m3:
				entry_name = m3.group(1)
				protein_name = m3.group(2)
				#print('Protein '+uniprot_id+' has no gene name')
				proteins_with_no_gene_names.add(uniprot_id)
				return (None,proteins_of_genes_not_to_collect,proteins_with_no_gene_names)
			else:
				print(line+" matches 1st pattern, but does not match 2nd and 3rd patterns in get_protein_header")
				return (None,proteins_of_genes_not_to_collect,proteins_with_no_gene_names)
	else:
		return (None,proteins_of_genes_not_to_collect,proteins_with_no_gene_names)

if __name__ == "__main__":
        main()


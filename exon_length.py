#Compute the feature Exon_length: the sum of length of all Exons of the longest transcript of each gene
#
#The Exon_length (sum of length of the Exons of the longest transcript of a gene) is equal to the length of the longest transcript of the gene.
#
#The length of the longest transcript of a gene can be retrieved from Ensemble using Biomart.
#
#To compute total exons length, for each gene find the longest transcript, then sum up the lengths of all the exons of the longest transcript to get Exon length.
#
#input: transcript_length_exon_rank_exon_start_end.csv
#output: a csv file containing Exon length of genes
#
#format of csv file
#
#Ensembl Gene ID,Ensembl Transcript ID,Exon Rank in Transcript,Transcript length,Exon Chr Start (bp),Exon Chr End (bp)
#ENSMUSG00000017390,ENSMUST00000017534,1,2786,78322968,78323222
#ENSMUSG00000017390,ENSMUST00000017534,2,2786,78324390,78324513
#ENSMUSG00000017390,ENSMUST00000017534,3,2786,78324604,78324815
#ENSMUSG00000017390,ENSMUST00000017534,4,2786,78324909,78324963
#ENSMUSG00000017390,ENSMUST00000017534,5,2786,78325049,78325209
#
#This program uses exon_count.py

import sys
import exon_count

def main():
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/test set/G6b_transcript_length_exon_rank_exon_start_end.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/G6b_exon_length.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/transcript_length_exon_rank_exon_start_end.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length.csv'	
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/transcript_length_exon_rank_exon_start_end1013.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length1013.csv'	
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/transcript_length_exon_rank_exon_start_end10.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length10.csv'	
	genesL = [line.strip() for line in open(infile)]
	#hash_longest_transcripts: key = gene id, value = transcript id of the longest transcript of the gene	
	hash_longest_transcripts = exon_count.get_longest_transcripts(genesL[1:len(genesL)])
	get_total_exon_length(hash_longest_transcripts,genesL[1:len(genesL)],outfile)
	#compute_total_exon_length(hash_longest_transcripts,genesL[1:len(genesL)],outfile)

def get_total_exon_length(hash_longest_transcripts,genesL,outfile):
	#Retrieve the length of the longest transcript of a gene from transcript_length_exon_rank_exon_start_end.csv
	#total exons length is the transcript length of the longest transcript of a gene
	hash_exon_length = {} #key = transcribe id of a longest transcript of a gene, value = total exon length of the longest transcript of the gene (exon length)
	longest_transcripts_ids = list(hash_longest_transcripts.values())
	for gene in genesL:
		geneL = gene.split(',')
		transcript_id = geneL[1]
		if transcript_id in longest_transcripts_ids:
			hash_exon_length[transcript_id] = geneL[3]
	fw=open(outfile,'w')
	fw.write("Ensemble_ID,ExonLength\n")
	gene_ids = list(hash_longest_transcripts.keys())
	for gene_id in gene_ids:
		transcript_id = hash_longest_transcripts[gene_id]
		fw.write(gene_id+','+str(hash_exon_length[transcript_id])+'\n')
	fw.close()
		
def compute_total_exon_length(hash_longest_transcripts,genesL,outfile):
	hash_exon_start_end_bps_of_longest_transcripts = {} #key = transcript id of a longest transcript, value = set of tuples (exon start bp, exon end bp) of the transcript
	longest_transcripts_ids = list(hash_longest_transcripts.values())
	for gene in genesL:
		geneL = gene.split(',')
		transcript_id = geneL[1]
		if transcript_id in longest_transcripts_ids:
			if hash_exon_start_end_bps_of_longest_transcripts.get(transcript_id)!= None:
				exon_start_end_bps = hash_exon_start_end_bps_of_longest_transcripts[transcript_id]
				exon_start_end_bps.add((int(geneL[4]),int(geneL[5])))
				hash_exon_start_end_bps_of_longest_transcripts[transcript_id] = exon_start_end_bps
			else:
				exon_start_end_bps = set()
				exon_start_end_bps.add((int(geneL[4]),int(geneL[5])))
				hash_exon_start_end_bps_of_longest_transcripts[transcript_id] = exon_start_end_bps
	#print(hash_exon_start_end_bps_of_longest_transcripts)
	hash_exon_length = {} #key = transcribe id of a longest transcript of a gene, value = total exon length of the longest transcript of the gene (exon length)
	for transcript_id in longest_transcripts_ids:
		exon_start_end_bps_Set = hash_exon_start_end_bps_of_longest_transcripts[transcript_id]
		exon_start_end_bpsL = list(exon_start_end_bps_Set)
		total_exons_length = 0
		print(transcript_id)
		for exon_start_end_bps in exon_start_end_bpsL:
			start_bp = exon_start_end_bps[0]
			end_bp = exon_start_end_bps[1]
			exon_length = end_bp - start_bp + 1
			total_exons_length += exon_length
			print(str(end_bp) +' - '+str(start_bp) + '+1')
		print('total exon length: '+str(total_exons_length))
		hash_exon_length[transcript_id] = total_exons_length
	fw=open(outfile,'w')
	fw.write("Ensemble_ID,ExonLength\n")
	gene_ids = list(hash_longest_transcripts.keys())
	for gene_id in gene_ids:
		transcript_id = hash_longest_transcripts[gene_id]
		fw.write(gene_id+','+str(hash_exon_length[transcript_id])+'\n')
	fw.close()

if __name__ == "__main__":
	main()	

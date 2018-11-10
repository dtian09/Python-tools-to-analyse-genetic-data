#compute Exon count of genes from 'transcript_length_exon_rank_exon_start_end.csv'.
#For each gene find the longest transcript, then find the largest Exon rank (Exon count) of the longest transcript. 
#
#input: transcript_length_exon_rank_exon_start_end.csv
#output: a csv file containing Exon count of genes
#
#format of csv file
#
#Ensembl Gene ID,Ensembl Transcript ID,Exon Rank in Transcript,Transcript length,Exon Chr Start (bp),Exon Chr End (bp)
#ENSMUSG00000017390,ENSMUST00000017534,1,2786,78322968,78323222
#ENSMUSG00000017390,ENSMUST00000017534,2,2786,78324390,78324513
#ENSMUSG00000017390,ENSMUST00000017534,3,2786,78324604,78324815
#ENSMUSG00000017390,ENSMUST00000017534,4,2786,78324909,78324963
#ENSMUSG00000017390,ENSMUST00000017534,5,2786,78325049,78325209

import sys

def main():
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/test set/G6b_transcript_length_exon_rank_exon_start_end.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/G6b_exon_count.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/transcript_length_exon_rank_exon_start_end.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_count.csv'
	#infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/transcript_length_exon_rank_exon_start_end1013.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_count1013.csv'
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/transcript_length_exon_rank_exon_start_end10.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_count10.csv'
	genesL = [line.strip() for line in open(infile)]
	#hash_longest_transcripts: key = gene id, value = transcript id of the longest transcript of the gene
	hash_longest_transcripts = get_longest_transcripts(genesL[1:len(genesL)])
	#print(hash_longest_transcripts)
	exon_count(hash_longest_transcripts,genesL[1:len(genesL)],outfile)
	
def get_longest_transcripts(genesL):
	#input: a list with each element a line in a csv file
	#output: hash_longest_transcripts: key = gene id, value = transcript id of the longest transcript of the gene
	#
	hash_transcripts_length = {} #key = gene id, value = a set of tuples (transcript id of a transcript,length of the transcript)
	for gene in genesL:
		geneL = gene.split(',')
		gene_id = geneL[0]
		transcript_id = geneL[1]
		if gene_id == '?' or gene_id == '' or gene_id == 'none':
			print('missing gene_id')
		if transcript_id == '?' or transcript_id == '' or transcript_id == 'none':
			print('missing transcript_id')
		transcript_length = geneL[3]
		if hash_transcripts_length.get(gene_id) != None:
			transcripts_length = hash_transcripts_length[gene_id]
			transcripts_length.add((transcript_id,int(transcript_length)))
			hash_transcripts_length[gene_id] = transcripts_length
		else:
			transcripts_length = set()
			transcripts_length.add((transcript_id,int(transcript_length)))
			hash_transcripts_length[gene_id] = transcripts_length
	#print(hash_transcripts_length)			
	hash_longest_transcripts = {}
	gene_ids = list(hash_transcripts_length.keys())
	for gene_id in gene_ids:
		transcripts_lengthSet = hash_transcripts_length[gene_id]
		transcripts_lengthL = list(transcripts_lengthSet)
		longest_transcript_tuple = transcripts_lengthL[0]
		longest_transcript_list = [longest_transcript_tuple[0],longest_transcript_tuple[1]]
		for next_transcript in transcripts_lengthL[1:len(transcripts_lengthL)]:
			if longest_transcript_list[1] < next_transcript[1]:
	 			longest_transcript_list[0] = next_transcript[0]
				longest_transcript_list[1] = next_transcript[1]
		hash_longest_transcripts[gene_id] = longest_transcript_list[0] 
	return hash_longest_transcripts

def exon_count(hash_longest_transcripts,genesL,outfile):
	hash_exon_ranks_of_longest_transcripts = {} #key = transcript id of a longest transcript, value = set of exon ranks of the transcript
	longest_transcripts_ids = list(hash_longest_transcripts.values())
	for gene in genesL:
		geneL = gene.split(',')
		transcript_id = geneL[1]
		if transcript_id == '?' or transcript_id == '' or transcript_id == 'none':
			print('missing transcript_id')
		if transcript_id in longest_transcripts_ids:
			if hash_exon_ranks_of_longest_transcripts.get(transcript_id)!= None:
				exon_ranks = hash_exon_ranks_of_longest_transcripts[transcript_id]
				exon_ranks.add(int(geneL[2]))
				hash_exon_ranks_of_longest_transcripts[transcript_id] = exon_ranks
			else:
				exon_ranks = set()
				exon_ranks.add(int(geneL[2]))
				hash_exon_ranks_of_longest_transcripts[transcript_id] = exon_ranks
	#print(hash_exon_ranks_of_longest_transcripts)
	hash_exon_counts = {} #key = transcribe id of a longest transcript of a gene, value = largest exon rank of the longest transcript of the gene (exon count)
	for transcript_id in longest_transcripts_ids:
		exon_ranksSet = hash_exon_ranks_of_longest_transcripts[transcript_id]
		exon_ranksL = list(exon_ranksSet)
		largest_exon_rank = exon_ranksL[0]
		for next_exon_rank in exon_ranksL[1:len(exon_ranksL)]:
			if next_exon_rank > largest_exon_rank:
				largest_exon_rank = next_exon_rank
		hash_exon_counts[transcript_id] = largest_exon_rank
	#print(hash_exon_counts)
	fw=open(outfile,'w')
	fw.write("Ensemble_ID,Exon_Count\n")
	gene_ids = list(hash_longest_transcripts.keys())
	for gene_id in gene_ids:
		transcript_id = hash_longest_transcripts[gene_id]
		fw.write(gene_id+','+str(hash_exon_counts[transcript_id])+'\n')
	fw.close()
	
if __name__ == "__main__":
	main()	

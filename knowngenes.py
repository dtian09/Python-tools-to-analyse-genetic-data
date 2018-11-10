#program to get the known essentiality genes which are either essential or non-essential
#input: a fasta file which may contain known essentiality genes and unknown essentiality genes
#       a file containing gene names of either essential genes or non-essential genes
#output: a fasta file containing only the gene names of either essential genes or non-essential genes and their protein sequences
#
# fasta sequence header format:
#>tr|A2AHY8|A2AHY8_MOUSE Pet2 protein OS=Mus musculus GN=Pet2 PE=2 SV=1
#author: David Tian
import sys
import re
import genes

def main():
	fasta_fileIn = sys.argv[1]
	gene_names_file = sys.argv[2]
	fasta_fileOut = sys.argv[3]
	
	fastaFile = [line.strip() for line in open(fasta_fileIn)]
	genesNamesToKeep = [line2.strip() for line2 in open(gene_names_file)]
	genesNamesToKeep = set(genesNamesToKeep)
	sequence_part = False
	fw = open(fasta_fileOut,"w")
	if len(genesNamesToKeep) == 0:
		print("No genes to keep and all genes are removed from "+fasta_fileIn)
		fw.close()
		sys.exit()
	i=0
	while i<len(fastaFile):
		gene_name = get_gene_name(fastaFile[i])
		#print(fastaFile[i])
		if gene_name != None:#line is a protein sequence header line and contains a gene name
			if gene_name in genesNamesToKeep:#gene name should be kept
				fw.write(fastaFile[i]+"\n")
				i += 1				
				sequence_part = a_part_of_sequence(fastaFile[i])
				while sequence_part == True and i < len(fastaFile):
					if i < (len(fastaFile)-1):
						fw.write(fastaFile[i]+"\n")	
					else:
						fw.write(fastaFile[i])
					i += 1
					if  i < len(fastaFile):
						sequence_part = a_part_of_sequence(fastaFile[i])
			else:
				i += 1	
		else:
			i += 1
		#	print("This line "+fastaFile[i]+" is not a protein sequence header line")
		#	sys.exit()
	fw.close()

def a_part_of_sequence(line):
#input: index of a line
#	a list with each element a line of a fasta file
#output: True (the line is a part of a sequence) or False (otherwise)
	matchObj = re.match('[\w]+', line)
	if matchObj!= None:
		return True
	else:
		return False

def get_gene_name(line):
#fasta sequence header format:
#>tr|A2AHY8|A2AHY8_MOUSE Pet2 protein OS=Mus musculus GN=Pet2 PE=2 SV=1
#input: a line of a fasta file
#output: gene name or None
	matchObj = re.match('>.+\\|.+(\\|.+)', line)
	if matchObj!= None:
		subStr = matchObj.group(1)
		matchObj = re.match('\\|[\w]+_[\w]+\\s.+\\sOS=.+\\s+GN=(.+)\\s+PE.+',subStr)
		if matchObj != None:
			gene_name = matchObj.group(1)
			return gene_name
		else:
			print("does not match: "+subStr)
			return None
	else:
		return None	

if __name__ == "__main__":
        main()				


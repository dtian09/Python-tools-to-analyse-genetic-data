#Update a data set with charge feature in a pepstats output file
#input:  a data set containing protein sequences features including the charge feature
#	 a pepstats output file containing properties of protein sequences 
#output: the data set with the charge feature
#
#pepstats file format: 
#start of a protein sequence: PEPSTATS of A2AHY8_MOUSE from 1 to 747
#
#some genes encode more than one proteins, the properties of the longest protein are used
#author: David Tian
import genes
import re
import sys

def main():
	data_file = sys.argv[1]
	fasta_file = sys.argv[2]
	pepstats_file = sys.argv[3]
	data_file2 = sys.argv[4]

	data = [line.strip() for line in open(data_file)]
	pepstats = [line2.strip() for line2 in open(pepstats_file)]
	fr = open(fasta_file,"r")
	fw = open(data_file2,"w")
	fw.write("Gene Name,UniProtID,Protein Name,MW,residues,charge,IP,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Aliphatic,Aromatic,NonPolar,Polar,Charged,basic,Acidic\n");
	proteins = genes.get_proteins(fr)#returns a dictionary: key=gene_name, value=set of tuples (protein name, UniProt id, UniProt entry name)
	proteins_positions = {}#a dictionary: key=uniprot entry name, value = position of header line of a sequence
	proteins_positions = get_positions_of_sequences(pepstats,proteins_positions)
	i=1 #ith gene that has been read
	gene2=''
	while i<len(data): #the header line is at index 0
		line = data[i]
		gene = line.split(',')
		gene_name = gene[0]	
		protein_name = gene[2]
		proteins_of_gene = proteins[gene_name]
		found_protein=False
		
		for protein_of_gene in proteins_of_gene:
			gene_protein = re.sub(',',":",protein_of_gene[0])
			if gene_protein == protein_name:
				found_protein = True
				uniprot_entry_name = protein_of_gene[2]
				position = proteins_positions[uniprot_entry_name]
				#m=re.match("Q3U983_MOUSE",uniprot_entry_name)
				#if uniprot_entry_name == "Q3U983_MOUSE": 
				#	print(position)
				#if protein_name == "Rab GTPase-binding effector protein 1":
				#	print(uniprot_entry_name)
				#	print(position)	
				charge = get_charge(position,pepstats)
				gene2=''
				a=0
				while a<5:
					gene2 += gene[a]+","
					a+=1
				gene2 += charge+"," #feature 'charge' has index 5
				b=6
				while b<(len(gene)-1):
					gene2 += gene[b]+","
					b+=1
				gene2 += gene[len(gene)-1]	
				break
		if found_protein == False:
			print('gene: '+gene_name+" does not have protein "+protein_name+" in dictionary")
		if i<(len(data)-1):
			fw.write(gene2+'\n')
		else:
			fw.write(gene2)
		i+=1		
	fw.close()

def get_positions_of_sequences(pepstats,proteins_positions):
#input: pepstats (a list with each element a line in a pepstats output file)
#	a dictionary: key=uniprot entry name, value = position of header line of a sequence in the pepstats list	
#output: the dictionary with positions of sequences added
	i=0
	while i<2526865: #position of header line of last protein sequence
		matchObj = re.match("PEPSTATS\\sof\\s(.+)\\sfrom.+",pepstats[i])
		entry_name = matchObj.group(1)
		#if entry_name == 'Q3U983_MOUSE':
		#	print(i)
		if proteins_positions.get(entry_name) == None:
			proteins_positions[entry_name] = i
		#else:
			#print("key: "+entry_name+" exist")
		i+=48
	return proteins_positions
	
def get_charge(position,pepstats):
#input: position of header line of a sequence in the pepstats list
#       pepstats (a list with each element a line in a pepstats output file)
#output: charge of the sequence
	line = pepstats[position+3]
	#Average Residue Weight  = 115.214 	Charge   = -0.5  
	matchObj = re.match("Average\\s+Residue\\s+Weight\\s+=\\s+[0-9\\.]+\\s+Charge\\s+=\\s+([\\-0-9\\.]+)\\s*",line)
	return matchObj.group(1)
	
if __name__ == "__main__":
	main()
				


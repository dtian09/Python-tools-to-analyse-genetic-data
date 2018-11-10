#program to merge all gene features with all protein features of unknown genes
#input: merged_gene_features.csv (MGI ID, Ensemble ID) (output by merge_gene_features.py)
#	merged_protein_features.csv (Uniprot ID, gene names) (output by merge_protein_features.py)
#	genenames16999_mgiids.txt
#	Training_All_GenesInfo_balanced1.csv
#output: merged_gene_protein_features.csv

import re
import sys

def main():
	infile_gene_features = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/merged_gene_features.csv'
	infile_protein_features = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/merged_protein_features.csv'
	infile_genenames_mgiids = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/unknown genes16999/genenames16999_mgiids.txt'
	infile_training_data ='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/merged_gene_protein_features_with_ids.csv'	

	gene_featuresL = [line.strip() for line in open(infile_gene_features)]

	if len(gene_featuresL)==0:
		print(infile_gene_features+" is empty.")
		sys.exit(-1)
	
	protein_featuresL = [line.strip() for line in open(infile_protein_features)]

	if len(protein_featuresL)==0:
		print(infile_protein_features+" is empty.")
		sys.exit(-1)
	
	genenames_mgiidsL = [line.strip() for line in open(infile_genenames_mgiids)]

	if len(genenames_mgiidsL)==0:
		print(infile_genesnames_mgiids+" is empty.")
		sys.exit(-1)

	training_dataL = [line.strip() for line in open(infile_training_data)]

	if len(training_dataL)==0:
		print(infile_training_data+" is empty.")
		sys.exit(-1)


	gene_features_hash = get_features(gene_featuresL,genenames_mgiidsL)#key=gene name, value= gene features including MGI ID (no missing values)
	protein_features_hash = get_protein_features(protein_featuresL)#key=gene name, value= protein features including Uniprot ID
	#print(protein_features_hash)
	merge_features(gene_features_hash,protein_features_hash,gene_featuresL[0],protein_featuresL[0],training_dataL[0],outfile)
	
def get_features(gene_featuresL,genenames_mgiidsL):
	#put features of genes in merged_gene_features.csv in a hashtable
	mgiids_hash = {}#key = mgi id, value = gene features
	mgiids_genenames_hash = {}#key = mgi id, value = gene name
	gene_features_hash = {}#key = gene name, value = gene features including MGI ID (no missing values)
	for line in gene_featuresL[1:len(gene_featuresL)]:
		m = re.match('^([:\w]+),([?\w]+),(.+)$',line)#e.g. MGI:2447322,ENSMUSG00000103770,189778,39.81,...
		if m:
			mgi_id = m.group(1)		
			mgiids_hash[mgi_id] = line
		else:
			print(line+' does not match pattern1 in get_features.')
	for line in genenames_mgiidsL[1:len(genenames_mgiidsL)]:
		m = re.match('^([\-\.\w]+)\s+([:\w]+)$',line)
		if m:
			gene_name = m.group(1)
			mgiid = m.group(2)
			mgiids_genenames_hash[mgiid] = gene_name
		else:
			print(line+' does not match pattern2 in get_features.')
	#mgiids = list(mgiids_genenames_hash.keys())#mgi ids in genenames16999_mgiids.txt
	mgiids = list(mgiids_hash.keys())  	    #mgi ids in merged_gene_features.csv 
	for mgiid in mgiids:
		#get the gene name of the mgi id
		if mgiids_genenames_hash.get(mgiid) != None:
			genename = mgiids_genenames_hash[mgiid]
		else:
			print('mgiids_genenames_hash has no gene name for '+mgiid+'\n')
			genename = '?'
		#get the gene features of the mgi id
		if mgiids_hash.get(mgiid)!= None:
			gene_features_hash[genename] = mgiids_hash[mgiid]
		else:
			print('mgiids_hash has no information for '+mgiid+'\n')
			line = gene_featuresL[0]
			vals = line.split(',')
			i = 1
			missing_vals = '?'
			while i <= len(vals)-1:
				missing_vals += ',?'
				i += 1
			gene_features_hash[genename] = missing_vals
	return gene_features_hash

def get_protein_features(protein_featuresL):
	hash_protein_features = {} #key=gene name, value= protein features including Uniprot ID
	for line in protein_featuresL[1:len(protein_featuresL)]:
		m = re.match('^[-\.\w]+,([-\.\s\(\)\w]+),.+$',line)#e.g. Q8BH24 (Uniprot id),Tm9sf4 (gene name),74693.22,643,...
		if m:
			genename = m.group(1)
			hash_protein_features[genename] = line
		else:
			print(line+' does not match pattern in get_protein_features.')

	return hash_protein_features

def merge_features(gene_features_hash,protein_features_hash,gene_features_header,protein_features_header,training_data_header,outfile):
	gene_features = gene_features_header.split(',')
	protein_features = protein_features_header.split(',')
	hash_gene_features_indx = {}
	hash_protein_features_indx = {}
	#get the indices of gene features in merged_gene_features.csv
	i=0
	for gene_f in gene_features:
		hash_gene_features_indx[gene_f]=i
		i += 1
	#get the indices of protein features in merged_protein_features.csv
	j=0
	for protein_f in protein_features:
		hash_protein_features_indx[protein_f] = j
		j += 1
	fw = open(outfile,'w')
	fw.write(training_data_header+'\n')
	#get the gene features at the indices in merged_gene_features.csv and the protein features at the indices merged_protein_features.csv
	genenames = list(protein_features_hash.keys())
	for genename in genenames:
		if gene_features_hash.get(genename) != None:
			gene_features = gene_features_hash[genename]
			gene_featuresL = gene_features.split(',')		
		else:
			gene_featuresL = []
			for gene_f in gene_features:
				gene_featuresL.append('?')	
		protein_features = protein_features_hash[genename]
		protein_featuresL = protein_features.split(',')
		
		l = protein_featuresL[hash_protein_features_indx['GeneName']]+','+gene_featuresL[hash_gene_features_indx['MGI_ID']]+','+protein_featuresL[hash_protein_features_indx['UniProt_ID']]+','+gene_featuresL[hash_gene_features_indx['Ensemble_ID']]
		features = 'GeneLength,GC_content,Transcript_count,Exon_Count,ExonLength,IntronLength'
		features_g = features.split(',')
		for f in features_g:
			l += ','+gene_featuresL[hash_gene_features_indx[f]]
		featuresStr = 'MW,ProteinLength,Aliphatic,Aromatic,NonPolar,Polar,Charged,Basic,Acidic,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,GlycoProtein,Phosphoprotein,Acetylation,Transcription,SignalPeptide,1.-.-.-,2.-.-.-,3.-.-.-,4.-.-.-,5.-.-.-,6.-.-.-,Nuclues_Score,Cytoplasm_Score,Plasma_Score,Extracellular_Score,Golgi_Score,ER_Score,Mitochondria_Score,Peroxisome_Score,Lysosome_Score,Nucleus_UniProt,Cytoplasm_UniProt,Plasma_UniProt,Membrane_UniProt,Extracellular_UniProt,Mitochondrion_UniProt,ER_UniProt,Golgi_UniProt,Lysosome_UniProt,Peroxisome_UniProt,CellJunction_Uniprot,CellProjection_Uniprot,Transmembrane_Count,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,BN_Score_Known,EPC_Score_Known,MNC_Score_Known,DMNC_Score_Known,ASP_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,Degree_KP,TC_KP,BN_Score_KP,EPC_Score_KP,MNC_Score_KP,DMNC_Score_KP'
		features_p = featuresStr.split(',')
		for f in features_p:
			l += ','+protein_featuresL[hash_protein_features_indx[f]]
		featuresStr2 = 'Oocyte(Transcript/Million),Unfertilized_Ovum(Transcript/Million),Zygote(Transcript/Million),Cleavage(Transcript/Million),Morula(Transcript/Million),Blastocyst(Transcript/Million),Egg_Cylinder(Transcript/Million),Gastrula(Transcript/Million),Organogenesis(Transcript/Million),Fetus(Transcript/Million),Neonate(Transcript/Million),Juveline(Transcript/Million),Adult(Transcript/Million),Age'
		features_g = featuresStr2.split(',')
		for f in features_g:
			l += ','+gene_featuresL[hash_gene_features_indx[f]]

		fw.write(l+',?\n')
	fw.close()
		
if __name__ == "__main__":
	main()

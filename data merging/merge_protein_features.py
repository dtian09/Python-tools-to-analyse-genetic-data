#program to merge all the protein features of unknown genes using Uniprot ID
#
#input: uniprot-mapped_unknown_essentiality_genes.pepstats.csv
#	unknown-essentiality genes_EC_keywords.csv
#	unknown_genes_subcell_loc_same_locs_as_known_genes.csv
#	signalp.csv
#	transmembrane_count.csv
#	merged_ppi_features.csv
#
#ouput: merged_protein_features.csv

import sys
import re

def main():
	infile_pepstats = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/uniprot-mapped_unknown_essentiality_genes.pepstats.csv'
	infile_EC_keywords = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/EC and keywords/unknown-essentiality genes_EC_keywords.csv'
	infile_subcell_locs = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/subcell loc/unknown_genes_subcell_loc_same_locs_as_known_genes.csv'
	infile_signalP = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/signal peptide/signalp.csv'
	infile_transmembrane_count ='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/transmembrane count/transmembrane_count.csv'
	infile_ppi = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/PPI network features/merged_ppi_features.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/merged_protein_features.csv'
	
	infile_pepstatsL = [line.strip() for line in open(infile_pepstats)]

	if len(infile_pepstatsL)==0:
		print(infile_pepstats+" is empty.")
		sys.exit(-1)
	
	infile_EC_keywordsL = [line.strip() for line in open(infile_EC_keywords)]

	if len(infile_EC_keywordsL)==0:
		print(infile_EC_keywords+" is empty.")
		sys.exit(-1)

	infile_subcell_locsL = [line.strip() for line in open(infile_subcell_locs)]

	if len(infile_subcell_locsL)==0:
		print(infile_subcell_locs+" is empty.")
		sys.exit(-1)

	infile_signalPL = [line.strip() for line in open(infile_signalP)]

	if len(infile_signalPL)==0:
		print(infile_signalP+" is empty.")
		sys.exit(-1)

	infile_transmembrane_countL = [line.strip() for line in open(infile_transmembrane_count)]

	if len(infile_transmembrane_countL)==0:
		print(infile_transmembrane_count+" is empty.")
		sys.exit(-1)
	
	infile_ppiL = [line.strip() for line in open(infile_ppi)]

	if len(infile_ppiL)==0:
		print(infile_ppi+" is empty.")
		sys.exit(-1)
	
	hash_pepstats = {}
	hash_EC_keywords = {}
	hash_subcell_locs = {}
	hash_signalP = {}
	hash_transmembrane_count = {}
	hash_ppi = {}
	
	hash_pepstats = get_pepstats_features(infile_pepstatsL,hash_pepstats)
	hash_EC_keywords = get_features(infile_EC_keywordsL,hash_EC_keywords)
	hash_subcell_locs = get_features(infile_subcell_locsL,hash_subcell_locs)
	hash_signalP = get_features(infile_signalPL,hash_signalP)
	hash_transmembrane_count = get_features(infile_transmembrane_countL,hash_transmembrane_count)
	hash_ppi = get_features(infile_ppiL,hash_ppi)
	
	fw = open(outfile,'w')
	features_names='UniProt_ID,GeneName,MW,ProteinLength,Charge,IP,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Aliphatic,Aromatic,NonPolar,Polar,Charged,Basic,Acidic'
	features_names = features_names+','+'1.-.-.-,2.-.-.-,3.-.-.-,4.-.-.-,5.-.-.-,6.-.-.-,Glycoprotein,Phosphoprotein,Acetylation,Transcription'
	features_names = features_names+','+'Nucleus_Uniprot,Cytoplasm_Uniprot,Plasma_Uniprot,Membrane_Uniprot,Extracellular_Uniprot,Mitochondrion_Uniprot,ER_Uniprot,Golgi_Uniprot,Lysosome_Uniprot,Peroxisome_Uniprot,CellJunction_Uniprot,CellProjection_Uniprot,Extracellular_Score,Cytoplasm_Score,Plasma_Score,Lysosome_Score,Peroxisome_Score,Mitochondria_Score,Nuclues_Score,ER_Score,Golgi_Score'
	features_names = features_names+','+'SignalPeptide'
	features_names = features_names+','+'Transmembrane_Count'
	features_names = features_names+','+'ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,ASP_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,TC_KP,BN_Score_Known,BN_Score_KP,Degree_KP,DMNC_Score_Known,DMNC_Score_KP,EPC_Score_Known,EPC_Score_KP,MNC_Score_Known,MNC_Score_KP'
	#debug
	#features_names = 'Uniprot ID,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,ASK_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,TC_KP,BN_known,BN_kp,Degree_kp,DMNC_known,DMNC_kp,EPC_known,EPC_kp,MNC_known,MNC_kp'

	EC_keywords=''
	subcell_locs =''
	signalP=''
	transmembrane_count=''
	ppi =''

	fw.write(features_names+'\n')
	ks = list(hash_pepstats.keys())
	for k in ks:
		if hash_EC_keywords.get(k) == None:
			names = 'EC_1,EC_2,EC_3,EC_4,EC_5,EC_6,Glycoprotein,Phosphoprotein,Acetylation,Transcription'
			namesL = names.split(',')
			for name in namesL[0:len(namesL)-1]:
				EC_keywords = EC_keywords+'?,'
			EC_keywords = EC_keywords+'?'
		else:
			EC_keywords = hash_EC_keywords[k]
		if hash_subcell_locs.get(k) == None:
			names = 'Nucleus_Uniprot,Cytoplasm_Uniprot,Plasma_Uniprot,Membrane_Uniprot,Extracellular_Uniprot,Mitochondrion_Uniprot,ER_Uniprot,Golgi_Uniprot,Lysosome_Uniprot,Peroxisome_Uniprot,CellJunction_Uniprot,CellProjection_Uniprot,Extracellular_Score,Cytoplasm_Score,Plasma_Score,Lysosome_Score,Peroxisome_Score,Mitochondria_Score,Nuclues_Score,ER_Score,Golgi_Score'
			namesL = names.split(',')
			for name in namesL[0:len(namesL)-1]:
				subcell_locs = subcell_locs+'?,'
			subcell_locs = subcell_locs+'?'
		else:
			subcell_locs = hash_subcell_locs[k]
		if hash_signalP.get(k) == None:
			signalP='?'
		else:
			signalP = hash_signalP[k]
		if hash_transmembrane_count.get(k) == None:
			transmembrane_count = '?'
		else:
			transmembrane_count = hash_transmembrane_count[k]
		if hash_ppi.get(k) == None:
			names = 'ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,ASK_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,TC_KP,BN_known,BN_kp,Degree_kp,DMNC_known,DMNC_kp,EPC_known,EPC_kp,MNC_known,MNC_kp'
			namesL = names.split(',')
			for name in namesL[0:len(namesL)-1]:
				ppi = ppi+'?,'
			ppi = ppi+'?'
		else:
			ppi = hash_ppi[k]
		fw.write(k+','+hash_pepstats[k]+','+EC_keywords+','+subcell_locs+','+signalP+','+transmembrane_count+','+ppi+'\n')
		EC_keywords = ''
		subcell_locs = ''
		ppi=''
	fw.close()

def get_pepstats_features(infileL,hash_features):
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		gene_name = vals[0]
		uniprot_id = vals[1]
		l = gene_name
		for val in vals[3:len(vals)]:#skip protein name field
			l = l+','+val
		hash_features[uniprot_id] = l	
	return hash_features

def get_features(infileL,hash_features):
	for line in infileL[1:len(infileL)]:
		m = re.match('^([\w]+),(.+)$',line)
		if m:
			uniprot_id = m.group(1)		
			hash_features[uniprot_id] = m.group(2)
		else:
			print(line+' does not match pattern in get_features.')
	return hash_features
	
if __name__ == "__main__":
	main()		
		

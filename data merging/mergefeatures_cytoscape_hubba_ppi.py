#program to merge cytoscape and Hubba ppi features
#input: merged_cytoscape_features.csv
#	merged_hubba_features.csv
#output: merged_ppi_features.csv
#
#format of merged_cytoscape_features.csv
#
#uniprot_id,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,ASP_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,TC_KP
#Q8BH24,-1,-1,-1,-1,-1,-1,4.05111168,0.00003333,0.24684582,0,0.25
#Q9JKY0,2.66666667,0,0.375,0,1,0,3.10375671,0.00027371,0.3221902,0.1952381,0.10029263
#
#format of merged_hubba_features.csv
#
#Uniprot_id,BN_known,BN_kp,Degree_kp,DMNC_known,DMNC_kp,EPC_known,EPC_kp,MNC_known,MNC_kp
#Q8C3Y4,0.00000,0.00000,1.00000,0.00000,0.00000,5.19400,1.02200,1.00000,1.00000
#Q3U9G9,0.00000,0.00000,19.00000,0.00000,0.27027,4.69200,12.17500,1.00000,14.00000
#Q9JKY0,0.00000,0.00000,21.00000,0.00000,0.33191,1.26100,15.18600,1.00000,17.00000

import sys
import re
import mergefeatures_cytoscape_ppi

def main():
	infile1 = sys.argv[1]#merged cytoscape ppi features
	infile2 = sys.argv[2]#merged hubba ppi features
	outfile = sys.argv[3]

	infile1L = [line.strip() for line in open(infile1)]

	if len(infile1L)==0:
		print(infile1+" is empty.")
		sys.exit(-1)
	
	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	
	hash_cytoscape_ppi = {}
	hash_hubba_ppi = {}

	hash_cytoscape_ppi = mergefeatures_cytoscape_ppi.get_ppi_features(infile1L,hash_cytoscape_ppi)
	hash_hubba_ppi = mergefeatures_cytoscape_ppi.get_ppi_features(infile2L,hash_hubba_ppi)
	fw = open(outfile,'w')
	fw.write('uniprot_id,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,ASK_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,TC_KP,BN_known,BN_kp,Degree_kp,DMNC_known,DMNC_kp,EPC_known,EPC_kp,MNC_known,MNC_kp\n')
	ks = list(hash_cytoscape_ppi.keys())
	for k in ks:
		if hash_hubba_ppi.get(k)!= None:
			fw.write(k+','+hash_cytoscape_ppi[k]+','+hash_hubba_ppi[k]+'\n')
		else:#protein has no interaction with other proteins and is not in the ppi network (not in the Hubba network features file infile2), its ppi features take missing values
			fw.write(k+','+hash_cytoscape_ppi[k]+',-?,-?,-?,-?,-?,-?,-?,-?,-?\n')
	fw.close()
	
if __name__ == "__main__":
	main()		
		

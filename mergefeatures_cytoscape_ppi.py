#program to merge known ppi features with kp features together
#input: known ppi features file
#	kp ppi features file
#output: merged_cytoscape_features.csv
#
#format of input files
#name,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known
#P29037,3.59259259,0.0507687,0.27835052,0.07142857,11,0.12890625
#P70670,3.36684303,0.02281426,0.29701414,0.33333333,4,0.34177215
#P84228,3.37213404,0.07069945,0.29654812,0.03846154,13,0.0951417
#
 
import collectfeatures_ppi

def main():
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15.csv'
	#infile1 = '/home/david/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15_cytoscape_known_ppi.csv'
	#infile2 = '/home/david/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15_cytoscape_kp_ppi.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15_cytoscape_ppi.csv'
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	infile1 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/known_ppi_cytoscape_features.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/kp_ppi_cytoscape_features.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/ppi_cytoscape_features.csv'	
	genesidL = [line.strip() for line in open(genesid_file)]
	infile1L = [line.strip() for line in open(infile1)]
	infile2L = [line.strip() for line in open(infile2)]
	hash_known_ppi = get_ppi_features(infile1L)
	#print(str(len(hash_known_ppi)))
	hash_kp_ppi = get_ppi_features(infile2L)
	#print(str(len(hash_kp_ppi)))
	uniprotids = collectfeatures_ppi.get_retrieved_longest_reviewed_proteins(genesidL)
	genenames_hash = get_genenames(genesidL)
	fw = open(outfile,'w')
	fw.write('GeneName,uniprot_id,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,ASP_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,Degree_KP,TC_KP\n')
	for uniprotid in list(uniprotids):
		if hash_known_ppi.get(uniprotid)!= None:
			known_ppis = hash_known_ppi[uniprotid]
		else:
			known_ppis = '?,?,?,?,?,?'#protein has no interaction with other proteins and is not in the ppi network (not in infile1), its ppi features take missing values
			#known_ppis = '-1,-1,-1,-1,-1,-1'
		if hash_kp_ppi.get(uniprotid)!= None:
			kp_ppis = hash_kp_ppi[uniprotid]
		else:
			kp_ppis = '?,?,?,?,?,?'
			#kp_ppis = '-1,-1,-1,-1,-1,-1'	
		fw.write(genenames_hash[uniprotid]+','+uniprotid+','+known_ppis+','+kp_ppis+'\n')
	fw.close()
	
def get_genenames(genesidL):
	genenames_hash = {}#key=uniprot id, value=gene name
	for line in genesidL[1:len(genesidL)]:
		vals = line.split(',')
		genename = vals[0]
		uniprotid = vals[2]
		genenames_hash[uniprotid]=genename
	return genenames_hash
		
def get_ppi_features(infileL):
	hash_ppi = {}
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		#SUID	AverageShortestPathLength	BetweennessCentrality	ClosenessCentrality
		#0	1				2			3
		#ClusteringCoefficient	Degree	Eccentricity	IsSingleNode	name
		#4			5	6		7		8
		#NeighborhoodConnectivity	NumberOfDirectedEdges	NumberOfUndirectedEdges
		#9				10			11
		#PartnerOfMultiEdgedNodePairs	Radiality	selected	SelfLoops	
		#12				13		14		15	
		#shared name	Stress	TopologicalCoefficient
		#16		17	18
		#features to collect: 		  	  
            #uniprot_id,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known
		uniprotid = vals[8]
		uniprotid = uniprotid.strip('"')
		asp = vals[1]
		bc = vals[2]
		closeness = vals[3]
		clustering = vals[4]
		degree = vals[5]
		tc = vals[18]
		hash_ppi[uniprotid] = asp+','+bc+','+closeness+','+clustering+','+degree+','+tc
	return hash_ppi
			
if __name__ == "__main__":
	main()		
		

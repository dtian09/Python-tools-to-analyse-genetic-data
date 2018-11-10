#compare features of 2 data files

def main():
	file1='/home/david/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15_cytoscape_ppi.csv'
	file2='/home/david/Dropbox/datasets/essential genes prediction/test set/BlindTestSet1_19Aug15_cytoscape_ppi_mitra.csv'
	outfile='/home/david/Dropbox/datasets/essential genes prediction/test set/compare_results.csv'
	file1L = [line.strip() for line in open(file1)]
	file2L = [line.strip() for line in open(file2)]
	(file1_features,hash_uniprotids1) = get_features(file1L)
	(file2_features,hash_uniprotids2) = get_features(file2L)
	names = file1L[0]
	namesL = names.split(',')
	gene_id = namesL[0]
	compare_features(names,gene_id,file1_features,file2_features,hash_uniprotids1,hash_uniprotids2,outfile)

def get_features(fileL):
	hash_features = {}
	hash_uniprotids = {}
	for line in fileL[1:len(fileL)]:
		vals = line.split(',')
		gene_id = vals[0]
		hash_features[gene_id] = vals[2:len(vals)]
		hash_uniprotids[gene_id] = vals[1]
	return (hash_features,hash_uniprotids)

def compare_features(names,gene_id,file1_features,file2_features,hash_uniprotids1,hash_uniprotids2,outfile):
	fw = open(outfile,'w')
	fw.write(names+','+names+','+names+'\n')
	gene_ids1 = set(file1_features.keys())
	gene_ids2 = set(file2_features.keys())
	if gene_ids1 == gene_ids2:
		print('file1 and file2 have identical '+gene_id)
	else:
		gene_ids_in_file1_not_in_file2 = gene_ids1.difference(gene_ids2)
		gene_ids_in_file2_not_in_file1 = gene_ids2.difference(gene_ids1)
		print('file2 does not contain '+str(len(gene_ids_in_file1_not_in_file2))+' ids of file1')
		print('file1 does not contain '+str(len(gene_ids_in_file2_not_in_file1))+' ids of file2')
	for gene_id in list(gene_ids1):
		features1 = file1_features[gene_id]
		features2 = file2_features[gene_id]
		uniprotid1 = hash_uniprotids1[gene_id]
		uniprotid2 = hash_uniprotids2[gene_id]
		#write features of file1
		fw.write(gene_id)
		fw.write(','+uniprotid1)
		i=0
		while i < len(features1):
			fw.write(','+features1[i])
			i += 1		
		#write features of file2
		fw.write(','+gene_id)
		fw.write(','+uniprotid2)
		i2=0
		while i2 < len(features2):
			fw.write(','+features2[i2])
			i2 += 1
		#compare features of file1 and file2	
		fw.write(','+gene_id)
		fw.write(','+uniprotid1)
		i3=0
		while i3 < len(features1):
			diff = float(features1[i3]) - float(features2[i3])
			diff = abs(diff)
			if features1[i3] == features2[i3]:
				fw.write(',0')
			elif diff < 0.01:#the difference is very small, the feature values are very similar. The difference may be due that different version of Cytoscape is used to compute the ppi features
				fw.write(',0')				
			else:
				fw.write(',1')
			i3+=1
		fw.write('\n')
	fw.close()

if __name__ == "__main__":
	main()						
			

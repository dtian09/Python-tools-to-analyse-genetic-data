#get the different and common annotations between the functional annotation clusters of the 539 human essential genes and the functional annotation clusters of the 539 mouse viable genes

import re
import sys

def main():
	clusters_annotations_human_genes='/home/david/Dropbox/datasets/essential genes prediction/human essential genes/DAVID functional annotation clustering/genes_of_different_essentiality_DAVIDfunctional annotation clustering_human_genes.csv'
	clusters_annotations_mouse_genes='/home/david/Dropbox/datasets/essential genes prediction/human essential genes/DAVID functional annotation clustering/genes_of_different_essentiality_DAVIDfunctional annotation clustering_mouse_genes.csv'	
	humanL = [line.strip() for line in open(clusters_annotations_human_genes)]
	mouseL = [line.strip() for line in open(clusters_annotations_mouse_genes)]
	(human_cluster_annots,human_cluster_annots_counts) = get_annotations(humanL)	
	(mouse_cluster_annots,mouse_cluster_annots_counts) = get_annotations(mouseL)
	common_annots = human_cluster_annots.intersection(mouse_cluster_annots)	
	human_cluster_annots2 = human_cluster_annots.difference(common_annots)
	mouse_cluster_annots2 = mouse_cluster_annots.difference(common_annots)
	#print(human_cluster_annots_counts)
	#print(mouse_cluster_annots_counts)
	print('functional annotation clustering of the 539 human genes and the 539 mouse genes with different essentiality\n')
	print('There are '+str(len(common_annots))+' common annotations of 539 human genes and 539 mouse genes:')
	print('annotation\tcount')	
	for annot in list(common_annots):
		count = human_cluster_annots_counts[annot] + mouse_cluster_annots_counts[annot]
		print(annot+'\t'+str(count))
	print('\nThe human genes clusters have '+str(len(human_cluster_annots2))+' annotations which the mouse genes clusters do not have:')
	print('annotation\tcount')		
	for annot in list(human_cluster_annots2):
		count = human_cluster_annots_counts[annot]
		print(annot+'\t'+str(count))
	print('\nThe mouse genes clusters have '+str(len(mouse_cluster_annots2))+' annotations which the human genes clusters do not have:')
	print('annotation\tcount')		
	for annot in list(mouse_cluster_annots2):
		count = mouse_cluster_annots_counts[annot]
		print(annot+'\t'+str(count))
	

def get_annotations(clusterannotL):
	cluster_annots = set()
	annots_counts = {}
	i=0
	j=1
	k=0
	while i < len(clusterannotL):
		line = clusterannotL[i]
		#print(str(i)+': '+line)
		if line == 'Category,Term,Count,%,PValue,Genes,List Total,Pop Hits,Pop Total,Fold Enrichment,Bonferroni,Benjamini,FDR':
			k=0
			i += 1
			line = clusterannotL[i]
			while line !=',,,,,,,,,,,,':
				vals = line.split(',')
				annot = vals[1]
				count = vals[2]
				m = re.match('^[0-9]+$',count)
				if m:
					cluster_annots.add(annot)		
					annots_counts[annot] = int(count)
				else:#get the next field until is is the count
					annot = annot+','+vals[2]
					for val in vals[3:len(vals)]:
						count = val
						m = re.match('^[0-9]+$',count)
						if m:
							cluster_annots.add(annot)
							annots_counts[annot] = int(count)
							break;
						else:
							annot = annot+','+val
				k += 1
				i += 1
				if i < len(clusterannotL):
					line = clusterannotL[i]
				else:
					break
			'''
			if line == ',,,,,,,,,,,,' or i == len(clusterannotL):
				print('cluster '+str(j)+' has '+str(k)+' annotations')
				j += 1
			'''
			i += 1
		else:
			i += 1
	return (cluster_annots,annots_counts)
	
if __name__ == "__main__":
	main()


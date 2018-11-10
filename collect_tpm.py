#Program to collect the Transcript Per Million at the 13 developmental stages of the unknown genes from Unigene 
#
#transcript per million = (gene EST * 1000000/total EST in pool)
# e.g. Gene EST/Total EST in pool = 1/19348
#      Transcript per Million = TRUNC(1*1000000/19348) 
#                             = TRUNC(51.685)
#                             = 51, which matches the value present in UniGene
#input:	genenames_mgiids_uniprotids_ensembleids_unigeneids.csv
#	Mm.profiles ((gene EST/total EST in pool) of 13 developmental stage)
#output: a data file with the ids and the transcript per million features
#
#format of Mm.profiles
#
#> 1|Developmental Stage
#oocyte	0 / 19348	
#unfertilized ovum	0 / 20312	
#zygote	0 / 28807	
#cleavage	1 / 27537	
#morula	0 / 36903	
#blastocyst	9 / 68210	
#egg cylinder	0 / 12123	
#gastrula	4 / 29408	
#organogenesis	8 / 130865	
#fetus	67 / 673862	
#neonate	20 / 108168	
#juvenile	32 / 286633	
#adult	76 / 1035996	

import sys
import re

def main():
	idsfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids_unigeneids.csv'
	Mmfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/transcript_per_million/Mm.profiles' #EST profiles
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/transcript_per_million/unknown_genes_transcript_per_million.csv'

	idsfileL = [line.strip() for line in open(idsfile)]

	if len(idsfileL)==0:
		print(idsfile+" is empty.")
		sys.exit(-1)

	profsL = [line.strip() for line in open(Mmfile)]

	if len(profsL)==0:
		print(Mmfile+" is empty.")
		sys.exit(-1)
	hash_unigeneids_positions = get_positions_of_unigeneids_in_Mmfile(profsL)
	#print(hash_unigeneids_positions)
	(n,k) = get_transcript_per_million_features(idsfileL,hash_unigeneids_positions,profsL,outfile)
	print('total no. of non-missing Unigene ids in idsfile: '+str(n))
	print(str(k)+' non-missing Unigene ids of idsfile have no EST profiles in Unigene and are not in Mm.profiles.')

def get_positions_of_unigeneids_in_Mmfile(profsL):
	hash_unigeneids_positions = {}
	i=0
	for line in profsL:
		m = re.match('^>\s+([\d]+)\s*\|\s*Developmental\s*Stage\s*$',line)
		if m:
			unigeneid = m.group(1)
			unigeneid = 'Mm.'+unigeneid
			hash_unigeneids_positions[unigeneid] = i	
		i += 1
	return hash_unigeneids_positions

def get_transcript_per_million_features(idsfileL,hash_unigeneids_positions,profsL,outfile):
	n = 0 #total no. of unigene ids in idsfile
	k = 0 #no. of unigene ids of idsfile not in the Mm.profiles
	fw=  open(outfile,'w')
	gene_ids = idsfileL[0]
	fw.write(gene_ids+',Oocyte(Transcript/Million),Unfertilized_Ovum(Transcript/Million),Zygote(Transcript/Million),Cleavage(Transcript/Million),Morula(Transcript/Million),Blastocyst(Transcript/Million),Egg_Cylinder(Transcript/Million),Gastrula(Transcript/Million),Organogenesis(Transcript/Million),Fetus(Transcript/Million),Neonate(Transcript/Million),Juveline(Transcript/Million),Adult(Transcript/Million)\n')
	for line in idsfileL[1:len(idsfileL)]:
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		uniprotid = vals[2]
		ensembleid = vals[3]
		unigeneid = vals[4]
		fw.write(genename+','+mgiid+','+uniprotid+','+ensembleid+','+unigeneid+',')
		if unigeneid != 'none':
			n += 1
			if hash_unigeneids_positions.get(unigeneid) != None:
				i = hash_unigeneids_positions[unigeneid]
				i += 1
				get_est_features_of_gene(i,profsL,fw)
			else:
				k += 1
				print(unigeneid+' has no EST profile in Unigene and is not in Mm.profiles.')
				j = 0
				while j < 12:
					fw.write('?,')#unigene Mm.profiles has no EST information about this Unigene id
					j += 1
				fw.write('?\n')
		else:#The gene name does not have a Unigene id
			j = 0
			while j < 12:
				fw.write('?,')
				j += 1
			fw.write('?\n')
	fw.close()
	return (n,k)
	
def get_est_features_of_gene(i,profsL,fw):
	#input: index of line containing oocyte feature e.g. "oocyte	0 / 19348"
	#	a list with each element a line in Mm.profiles
	#	output file
	#output: the index of the next line after the 13th (Gene EST/Total EST in pool) feature	
	line = profsL[i]
	m = re.match('^oocyte\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	#oocyte	0 / 19348	
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of oocyte is not matched at line: '+str(i))
	#unfertilized ovum	0 / 20312	
	i+=1
	line = profsL[i]
	m = re.match('^unfertilized\s+ovum\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of unfertilized is not matched at line: '+str(i))
	#zygote	0 / 28807	
	i+=1
	line = profsL[i]
	m = re.match('^zygote\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of zygote is not matched at line: '+str(i))
	#cleavage	1 / 27537	
	i+=1
	line = profsL[i]
	m = re.match('^cleavage\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of cleavage is not matched at line: '+str(i))
	#morula	0 / 36903
	i+=1
	line = profsL[i]
	m = re.match('^morula\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of morula is not matched at line: '+str(i))
	#blastocyst	9 / 68210	
	i+=1
	line = profsL[i]
	m = re.match('^blastocyst\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of blastocyst is not matched at line: '+str(i))	
	#egg cylinder	0 / 12123
	i+=1
	line = profsL[i]
	m = re.match('^egg\s+cylinder\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of egg cylinder is not matched at line: '+str(i))
	#gastrula	4 / 29408
	i+=1
	line = profsL[i]
	m = re.match('^gastrula\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of gastrula is not matched at line: '+str(i))
	#organogenesis	8 / 130865	
	i+=1
	line = profsL[i]
	m = re.match('^organogenesis\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of organogenesis is not matched at line: '+str(i))
	#fetus	67 / 673862
	i+=1
	line = profsL[i]
	m = re.match('^fetus\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of fetus is not matched at line: '+str(i))
	#neonate	20 / 108168	
	i+=1
	line = profsL[i]
	m = re.match('^neonate\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of neonate is not matched at line: '+str(i))
	#juvenile	32 / 286633	
	i+=1
	line = profsL[i]
	m = re.match('^juvenile\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+',')
	else:
		print('pattern of juvenile is not matched at line: '+str(i))
	#adult	76 / 1035996	
	i+=1
	line = profsL[i]
	m = re.match('^adult\s+([\d]+)\s*/\s+([\d]+)\s*$',line)
	if m:
		tpm = transcript_per_million(m.group(1),m.group(2))
		fw.write(str(tpm)+'\n')
	else:
		print('pattern of adault is not matched at line: '+str(i))
	i+=1#increment index to the index of next line
	return i

def transcript_per_million(gene_est,total_est_in_pool):
	return int(gene_est)*1000000/int(total_est_in_pool)

'''	
def get_ids(idsfileL):
	hash_unigeneids = {}
	for line in idsfileL[1:len(idsfileL)]:
		vals = line.split(',')				
		genename = vals[0]
		mgiid = vals[1]
		uniprotid = vals[2]
		ensembleid = vals[3]
		unigeneid = vals[4]
		if unigeneid != 'none' or unigeneid != '?':
			hash_unigeneids[genename] = unigeneid
'''		
		
'''
def get_genenames_unigene_ids(genename_indx,unigene_indx,idsfileL):
	hash_ids = {}
	all_genenames = set()
	genes_withno_unigene_ids = set()
	for line in idsfileL[1:len(idsfileL)]:
		valsL = line.split(',')
		genename = valsL[genename_indx]
		unigene_id = valsL[unigene_indx]
		if unigene_id != 'none' or unigene_id != '?':
			hash_ids[unigene_id] = genename
		else:
			genes_withno_unigene_ids.add(genename)
		all_genenames.add(genename)
	return (hash_ids,all_genenames,genes_withno_unigene_ids)
'''
'''
def get_est_features_of_genes(hash_unigeneids,profsL,fw):
	#input: hash_unigeneids: key=unigene id, value=gene name
	#	a list with each element a line in Mm.profiles
	#output: transcript per million features
	#	
	i=0 #index of line containing oocyte feature e.g. "oocyte	0 / 19348"
	unigene_ids_of_collected_genes = set()
	unigene_ids_of_genes_to_collect = set(list(hash_ids.keys()))
	while unigene_ids_of_genes_to_collect != set() and i < len(profsL):
		line = profsL[i]
		m=re.match('^>\s+([\d]+)\s*|\s*Developmental\s*stage\s*$',line)
		if m:
			unigene_id = m.group(1)
			unigene_id = 'Mm.'+unigene_id
			if unigene_id in unigene_ids_of_genes_to_collect:
				unigene_ids_of_genes_to_collect.remove(unigene_id)
				if hash_ids.get(unigene_id) != None:
					gene_name = hash_unigeneids[unigene_id]
					fw.write(gene_name+','+unigene_id+',')
					i+=1
					i = get_est_features_of_gene(i,profsL,fw)
					unigene_ids_of_collected_genes.add(unigene_id)
					fw.write('\n')
				else:
					print(unigene_id+' is not in Mm.profiles.\n')
			else:
				i+=1
		else:
			i+=1
	if len(unigene_ids_of_collected_genes) < len(unigene_ids_of_genes_to_collect):
		print('no. of genes to collect: '+str(len(unigene_ids_of_genes_to_collect))+', no. of genes collected: '+str(len(unigene_ids_of_collected_genes)))
		print('genes not collected: '+str(unigene_ids_of_genes_to_collect.difference(unigene_ids_of_collected_genes)))
	return fw
'''

if __name__ == "__main__":
	main()


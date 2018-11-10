#fill in missing ids of genes
import sys
import re

def main():
	
	t = 0 #total no. of missing ids in infile
	k = 0 #no. of missing ids have been filled in
	'''
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/more_genenames_mgiids_from_mgi.txt'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids2.csv'
	'''
	
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	infile2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_all_uniprotids2.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids2.csv'
		
	infileL = [line.strip() for line in open(infile)]
	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	infile2L = [line.strip() for line in open(infile2)]
	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)
	fw = open(outfile,'w')
	hash_ids = {}
	'''
	#Get MGI ids
	id_indx_infile = 'mgiid'
	for line in infile2L[1:len(infile2L)]:
		vals = line.split('\t')
		if len(vals)==1:
			print(line+' has no id.')
		genename = vals[0]
		if genename == '':
			print('gene name is empty.')
		genename = genename.lower()
		input_type = vals[1]
		if input_type == 'current symbol' or input_type == 'old symbol' or input_type == 'synonym' or input_type == 'related synonym':
			mgiid = vals[2]
			if mgiid == '':
				print('gene id is empty.')
			else:
				hash_ids[genename] = mgiid
	'''
	#get Uniprot ids
	id_indx_infile = 'uniprotid'
	for line in infile2L[1:len(infile2L)]:
		vals = line.split(',')
		genename = vals[0]
		uniprotid = vals[1]
		if genename == '':
			print('gene name is empty.')
		genename = genename.lower()
		if uniprotid == '':
			print('Uniprot id is empty.')
		if genename !='' and uniprotid != '':
			hash_ids[genename] = uniprotid
	
	#fill in missing gene ids
	fw.write(infileL[0]+'\n')
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		original_genename = vals[0]
		genename_lower = original_genename.lower()
		if id_indx_infile == 'mgiid':
			id_indx = 1
		elif id_indx_infile == 'uniprotid':
			id_indx = 2
		elif id_indx_infile == 'ensembleid':
			id_indx = 3
		else:
			print('invalid id_indx_infile: '+id_indx_infile)
		gene_id = vals[id_indx]
		if gene_id == 'none':
			t += 1
			if hash_ids.get(genename_lower) != None:
				k += 1
				gene_id = hash_ids[genename_lower]
				if id_indx_infile == 'mgiid':
					mgiid = gene_id
					new_line = original_genename+','+mgiid
					for val in vals[2:len(vals)]:
						new_line += ','+val
				elif id_indx_infile == 'uniprotid':
					mgiid = vals[1]
					uniprotid = gene_id
					new_line = original_genename+','+mgiid+','+uniprotid
					for val in vals[3:len(vals)]:
						new_line += ','+val
				elif id_indx_infile == 'ensembleid':
					mgiid = vals[1]
					uniprotid = vals[2]
					ensembleid = gene_id
					new_line = original_genename+','+mgiid+','+uniprotid+','+ensembleid
					for val in vals[4:len(vals)]:
						new_line += ','+val
				else:
					print('invalid id_indx_infile: '+id_indx_infile)
					new_line = 'invalid id_indx_infile'
				fw.write(new_line+'\n')
			else:
				fw.write(line+'\n')
		else:
			fw.write(line+'\n')
	fw.close()
	print('total no. of missing ids: '+str(t))
	print(str(k)+' missing ids have been filled in.')

if __name__ == "__main__":
	main()
	 	


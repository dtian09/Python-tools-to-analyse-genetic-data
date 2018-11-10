#merge gene features of unknown genes
#input: exon_count.csv (ensemble id)
#	exon_length.csv (ensemble id)
#	genelength.csv (mgi id)
#	intronlength.csv (mgi id)
#	GC_transcript_count.csv (mgi id)
#	transcript_length.csv (ensemble id)
#	est_data.csv (mgi id, unigene id)
#	age.csv (mgi id)
#	mgi_ids_ensemble_ids.csv
#output: merged_gene_features.csv
import sys
import re

def main():
	infile_exon_count = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_count.csv'
	infile_exon_length = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length.csv'
	infile_gene_length = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/genelength.csv'
	infile_intronlength = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/intronlength.csv'
	infile_GC_transcript_count = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/GC_transcript_count.csv'
	infile_transcript_length = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/transcript_length.csv'
	infile_est_data = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/EST features/est_data.csv'
	infile_age = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/age/age.csv'
	infile_mgiids_ensembleids = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/mgi_ids_ensemble_ids.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/merged_gene_features.csv'

	exon_countL = [line.strip() for line in open(infile_exon_count)]

	if len(exon_countL)==0:
		print(infile_exon_count+" is empty.")
		sys.exit(-1)
	
	exon_lengthL = [line.strip() for line in open(infile_exon_length)]

	if len(exon_lengthL)==0:
		print(infile_exon_length+" is empty.")
		sys.exit(-1)

	gene_lengthL = [line.strip() for line in open(infile_gene_length)]

	if len(gene_lengthL)==0:
		print(infile_gene_length+" is empty.")
		sys.exit(-1)

	intronlengthL = [line.strip() for line in open(infile_intronlength)]

	if len(intronlengthL)==0:
		print(infile_intronlength+" is empty.")
		sys.exit(-1)

	GC_transcript_countL = [line.strip() for line in open(infile_GC_transcript_count)]

	if len(GC_transcript_countL)==0:
		print(infile_GC_transcript_count+" is empty.")
		sys.exit(-1)
	
	transcript_lengthL = [line.strip() for line in open(infile_transcript_length)]

	if len(transcript_lengthL)==0:
		print(infile_transcript_length+" is empty.")
		sys.exit(-1)

	est_dataL = [line.strip() for line in open(infile_est_data)]

	if len(est_dataL)==0:
		print(infile_est_data+" is empty.")
		sys.exit(-1)

	ageL = [line.strip() for line in open(infile_age)]

	if len(ageL)==0:
		print(infile_age+" is empty.")
		sys.exit(-1)

	mgiids_ensembleidsL = [line.strip() for line in open(infile_mgiids_ensembleids)]

	if len(mgiids_ensembleidsL)==0:
		print(infile_mgiids_ensembleids+" is empty.")
		sys.exit(-1)

	exon_count_hash = get_features(exon_countL)#key = ensemble_id, value = exon count
	exon_length_hash = get_features(exon_lengthL)#key = ensemble_id, value = exon length
	gene_length_hash = get_features(gene_lengthL)#key = mgi id, value = gene length
	intron_length_hash = get_features(intronlengthL)#key = mgi id, value = intron length
	GC_transcript_count_hash = get_features(GC_transcript_countL)#key = mgi_id, value = GC content, transcript count
	transcript_length_hash = get_features(transcript_lengthL)#key=ensemble id, value = transcript length
	est_hash = get_features(est_dataL)#key=mgi id, value = est 
	age_hash = get_features(ageL)#key=mgi id, value = age
	mgiids_ensembleids_hash = get_features(mgiids_ensembleidsL)#key=mgi id, value = ensemble id
	
	mgiids = list(gene_length_hash.keys())
	fw = open(outfile,'w')
	fw.write('MGI_ID,Ensemble_ID,GeneLength,GC_content,Transcript_count,Exon_Count,ExonLength,IntronLength,Oocyte(Transcript/Million),Unfertilized_Ovum(Transcript/Million),Zygote(Transcript/Million),Cleavage(Transcript/Million),Morula(Transcript/Million),Blastocyst(Transcript/Million),Egg_Cylinder(Transcript/Million),Gastrula(Transcript/Million),Organogenesis(Transcript/Million),Fetus(Transcript/Million),Neonate(Transcript/Million),Juveline(Transcript/Million),Adult(Transcript/Million),Age\n')
	for mgiid in mgiids:
		fw.write(mgiid+',')
		if mgiids_ensembleids_hash.get(mgiid)!= None:
			ensemble_id = mgiids_ensembleids_hash[mgiid]	
		else:
			ensemble_id = '?' #missing id
		fw.write(ensemble_id+',')
		if gene_length_hash.get(mgiid)!= None:
			gene_length = gene_length_hash[mgiid]
		else:
			gene_length = '?' #missing value
		fw.write(gene_length+',')
		if GC_transcript_count_hash.get(mgiid)!= None:
			GC_transcript_count = GC_transcript_count_hash[mgiid]
			vals = GC_transcript_count.split(',')
			GC = vals[0]
			transcript_count = vals[1]	
		else:
			GC = '?'
			transcript_count = '?'
		fw.write(GC+','+transcript_count+',')
		if ensemble_id == '?':
			exon_count = '?'
			exon_length = '?'
		else:
			if exon_count_hash.get(ensemble_id)!= None:
				exon_count = exon_count_hash[ensemble_id]
			else:
				exon_count = '?'
			if exon_length_hash.get(ensemble_id)!= None:
				exon_length = exon_length_hash[ensemble_id]
			else:
				exon_length = '?'
		fw.write(exon_count+','+exon_length+',')
		if intron_length_hash.get(mgiid)!= None:
			intron_length = intron_length_hash[mgiid]
		else:
			intron_length = '?'
		fw.write(intron_length+',')
		if est_hash.get(mgiid)!= None:
			est = est_hash[mgiid]
			m = re.match('^([\.\w]+),(.+)$',est)
			if m:
				est_features = m.group(2)#get the est features
		else:
			est_features='?'
			i=1
			while i <= 12:
				est_features +=',?'
				i+=1
		fw.write(est_features+',')
		if age_hash.get(mgiid)!= None:
			age = age_hash[mgiid]
		else:
			age = '?'
		fw.write(age+'\n')
	fw.close()

def get_features(infileL):
	hash_features = {}
	for line in infileL[1:len(infileL)]:
		m = re.match('^([:\w]+),(.+)$',line)
		if m:
			gene_id = m.group(1)		
			hash_features[gene_id] = m.group(2)
		else:
			print(line+' does not match pattern in get_features.')
	return hash_features

if __name__ == "__main__":
	main()
	
	
	

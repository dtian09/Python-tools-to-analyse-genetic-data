#insert gene features into mouse_genes database
import MySQLdb

def main():
	db = MySQLdb.connect("localhost","david","david","mouse_genes")
	#update_ids(db)
	#insert_tpm_features(db)#insert transcript per million into transcript_per_million table	
	#insert_gene_length(db)
	#insert_exon_count(db)
	#insert_exon_length(db)
	#insert_GC_content_and_transcript_count(db)
	#insert_intron_length(db)
	#insert_pepstats_features(db)
	#insert_cytoscape_ppi_features(db)
	#insert_hubba_ppi_features(db)
	#insert_subcell_location_features(db)
	#insert_ec_keywords_features(db)
	#insert_signalpeptide_features(db)
	#insert_transmembrane_count_features(db)
	insert_age(db)
	db.close()

def insert_age(db):
	print('insert age feature')
	file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/age/age.csv'
	fileL = [line.strip() for line in open(file)]
	cursor = db.cursor()
	i=0
	for line in fileL[1:len(fileL)]:
		i += 1
		vals = line.split(',')	
		ensembleid = vals[0]
		age = vals[1]
		try:
			sql = "update gene_features set Age = '%s' where ensembleid='%s'" % (age,ensembleid)	
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()
			sql = "insert into gene_features(ensembleid,Age) values('%s','%s')" % (ensembleid,age)
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))

def insert_transmembrane_count_features(db):
	print('insert transmembrane count features')
	file='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/transmembrane count/transmembrane_count.csv'
	fileL = [line.strip() for line in open(file)]
	cursor = db.cursor()
	i=0
	for line in fileL[1:len(fileL)]:
		i += 1
		vals = line.split(',')	
		uniprotid = vals[0]
		t_count = vals[1]
		sql = "insert into transmembrane(uniprotid,Transmembrane_Count) values('%s', '%s')" % (uniprotid,t_count)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
		except:
			db.rollback()
			print('rollback '+str(i))

def insert_signalpeptide_features(db):
	print('insert signalpeptide features')
	sigpepfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/signal peptide/uniprot_signalp.csv'
	sigpepL = [line.strip() for line in open(sigpepfile)]
	cursor = db.cursor()
	i=0
	for line in sigpepL[1:len(sigpepL)]:
		i += 1
		vals = line.split(',')	
		uniprotid = vals[0]
		sigpep = vals[1]
		sql = "insert into signalpeptides(uniprotid,signalpeptide) values('%s', '%s')" % (uniprotid,sigpep)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
		except:
			db.rollback()
			print('rollback '+str(i))

def insert_ec_keywords_features(db):
	print('insert ec keywords features')
	ec_keywords='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/EC and keywords/EC_keywords.csv'
	ec_keywordsL = [line.strip() for line in open(ec_keywords)]
	cursor = db.cursor()
	i=0
	for line in ec_keywordsL[1:len(ec_keywordsL)]:
		i += 1
		vals = line.split(',')	
		uniprotid = vals[0]
		ec1 = vals[1]
		ec2 = vals[2]
		ec3 = vals[3]
		ec4 = vals[4]
		ec5 = vals[5]
		ec6 = vals[6]
		gly = vals[7]
		pho = vals[8]
		ace = vals[9]
		tra = vals[10]
		sql = "insert into EC_keywords_features(uniprotid,ec1,ec2,ec3,ec4,ec5,ec6,GlycoProtein,Phosphoprotein,Acetylation,Transcription) values('%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s')" % (uniprotid,ec1,ec2,ec3,ec4,ec5,ec6,gly,pho,ace,tra)
		#print(sql)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
		except:
			db.rollback()
			print('rollback '+str(i))

def insert_subcell_location_features(db):
	print('insert subcell location features')
	locfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/subcell loc/uniprot_wolfpsort_loc2.csv'
	locfileL = [line.strip() for line in open(locfile)]
	cursor = db.cursor()
	i=0
	for line in locfileL[1:len(locfileL)]:
		i += 1
		vals = line.split(',')	
		uniprotid = vals[0]
		Extracellular_Score = vals[1]
		Cytoplasm_Score = vals[2]
		Plasma_Score = vals[3]
		Lysosome_Score = vals[4]
		Peroxisome_Score = vals[5]
		Mitochondria_Score= vals[6]
		Nucleus_Score = vals[7]
		ER_Score = vals[8]
		Golgi_Score = vals[9]
		Nucleus_Uniprot	= vals[10]
		Cytoplasm_Uniprot = vals[11]
		Membrane_Uniprot = vals[12]
		Extracellular_Uniprot = vals[13]
		Mitochondrion_Uniprot = vals[14]
		ER_Uniprot = vals[15]
		Golgi_Uniprot = vals[16]
		Lysosome_Uniprot = vals[17]
		Peroxisome_Uniprot = vals[18]
		CellJunction_Uniprot = vals[19]
		CellProjection_Uniprot = vals[20]	
		sql = "insert into subcellular_localizations(uniprotid,Extracellular_Score,Cytoplasm_Score,Plasma_Score,Lysosome_Score,Peroxisome_Score,Mitochondria_Score,Nucleus_Score,ER_Score,Golgi_Score,Nucleus_Uniprot,Cytoplasm_Uniprot,Membrane_Uniprot,Extracellular_Uniprot,Mitochondrion_Uniprot,ER_Uniprot,Golgi_Uniprot,Lysosome_Uniprot,Peroxisome_Uniprot,CellJunction_Uniprot,CellProjection_Uniprot) \
values('%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s')" % \
(uniprotid,Extracellular_Score,Cytoplasm_Score,Plasma_Score,Lysosome_Score,Peroxisome_Score,Mitochondria_Score,Nucleus_Score,ER_Score,Golgi_Score,Nucleus_Uniprot,Cytoplasm_Uniprot,Membrane_Uniprot,Extracellular_Uniprot,Mitochondrion_Uniprot,ER_Uniprot,Golgi_Uniprot,Lysosome_Uniprot,Peroxisome_Uniprot,CellJunction_Uniprot,CellProjection_Uniprot)
		#print(sql)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
		except:
			db.rollback()
			print('rollback '+str(i))

def insert_cytoscape_ppi_features(db):
	print('insert cytoscape ppi features')
	ppifile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/ppi_cytoscape_features.csv'
	ppifileL = [line.strip() for line in open(ppifile)]
	cursor = db.cursor()
	i=0
	for line in ppifileL[1:len(ppifileL)]:
		i += 1
		vals = line.split(',')
		genename = vals[0]
		uniprotid = vals[1]
		asp_known = vals[2]
		bc_known = vals[3]
		closenesscentrality_known = vals[4]
		clusteringcoefficient_known = vals[5]
		degree_known = vals[6]
		tc_known = vals[7]
		asp_kp = vals[8]
		bc_kp = vals[9]
		closenesscentrality_kp = vals[10]
		clusteringcoefficient_kp = vals[11]
		degree_kp = vals[12]
		tc_kp = vals[13]
		sql = "insert into cytoscape_ppi_features(genename,uniprotid,ASP_Known,BC_Known,ClosenessCentrality_Known,ClusteringCoefficient_Known,Degree_Known,TC_Known,ASP_KP,BC_KP,ClosenessCentrality_KP,ClusteringCoefficient_KP,Degree_KP,TC_KP) \
values('%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s')" % \
(genename,uniprotid,asp_known,bc_known,closenesscentrality_known,clusteringcoefficient_known,degree_known,tc_known,asp_kp,bc_kp,closenesscentrality_kp,clusteringcoefficient_kp,degree_kp,tc_kp)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
		except:
			db.rollback()
			print('rollback '+str(i))
	
def insert_hubba_ppi_features(db):
	print('insert hubba ppi features')
	ppifile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/ppi_hubba_features.csv'
	ppifileL = [line.strip() for line in open(ppifile)]
	cursor = db.cursor()
	i=0
	for line in ppifileL[1:len(ppifileL)]:
		i += 1
		vals = line.split(',')
		genename = vals[0]
		uniprotid = vals[1]
		bn_known = vals[2]
		epc_known = vals[3]
		mnc_known = vals[4]
		dmnc_known = vals[5]
		bn_kp = vals[6]
		epc_kp = vals[7]
		mnc_kp = vals[8]
		dmnc_kp = vals[9]
		sql = "insert into hubba_ppi_features(genename,uniprotid,BN_Score_Known,EPC_Score_Known,MNC_Score_Known,DMNC_Score_Known,BN_Score_KP,EPC_Score_KP,MNC_score_KP,DMNC_Score_KP) \
values('%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s')" % \
(genename,uniprotid,bn_known,epc_known,mnc_known,dmnc_known,bn_kp,epc_kp,mnc_kp,dmnc_kp)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
		except:
			db.rollback()
			print('rollback '+str(i))

def insert_tpm_features(db):
	print('insert tpm features')
	tpmfile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/transcript_per_million/unknown_genes_transcript_per_million.csv'
	tpmL = [line.strip() for line in open(tpmfile)]
	cursor = db.cursor()
	i=0
	for line in tpmL[1:len(tpmL)]:
		i += 1
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		uniprotid = vals[2]
		ensembleid = vals[3]
		unigeneid = vals[4]
		Oocyte = vals[5]
		Unfertilized_Ovum = vals[6]
		Zygote = vals[7]
		Cleavage = vals[8]
		Morula = vals[9]
		Blastocyst = vals[10]
		Egg_Cylinder = vals[11]
		Gastrula = vals[12]
		Organogenesis = vals[13]
		Fetus = vals[14]
		Neonate = vals[15]
		Juveline = vals[16]
		Adult = vals[17]
		try:#try starts a database transaction
			sql = "UPDATE transcript_per_million SET mgiid = '%s', uniprotid = '%s', ensembleid = '%s', unigeneid = '%s', Oocyte = '%s', Unfertilized_Ovum = '%s', Zygote = '%s', Cleavage = '%s', Morula = '%s', Blastocyst = '%s', Egg_Cylinder = '%s', Gastrula = '%s', Organogenesis = '%s', Fetus = '%s', Neonate = '%s', Juveline = '%s', Adult = '%s' where genename='%s'" % (mgiid,uniprotid,ensembleid,unigeneid,Oocyte,Unfertilized_Ovum,Zygote,Cleavage,Morula,Blastocyst,Egg_Cylinder,Gastrula,Organogenesis,Fetus,Neonate,Juveline,Adult,genename)
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()#end the current transaction before start a new transaction
			sql = "INSERT INTO transcript_per_million(genename,mgiid,uniprotid,ensembleid,unigeneid,Oocyte,Unfertilized_Ovum,Zygote,Cleavage,Morula,Blastocyst,Egg_Cylinder,Gastrula,Organogenesis,Fetus,Neonate,Juveline,Adult)\
     				VALUES ('%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s')" % \
       				(genename,mgiid,uniprotid,ensembleid,unigeneid,Oocyte,Unfertilized_Ovum,Zygote,Cleavage,Morula,Blastocyst,Egg_Cylinder,Gastrula,Organogenesis,Fetus,Neonate,Juveline,Adult)
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
	print('insert tpm finished')

def insert_gene_length(db):
	print('insert gene length')
	#file1 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length.csv'
	#file1 = '/home/david/Dropbox/datasets/essential genes prediction/test set/gene_length.csv'
	file1 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/gene_length10.csv'	
	cursor = db.cursor()
	i = 0
	genelengthL = [line.strip() for line in open(file1)]
	for line in genelengthL[1:len(genelengthL)]:
		i += 1
		vals = line.split(',')
		mgiid = vals[0]
		ensembleid = vals[1]
		genelength = vals[2]		
		try:
			sql = "update gene_features set ensembleid = '%s', genelength = '%s' where mgiid='%s'" % (ensembleid,genelength,mgiid)			
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()
			sql = "insert into gene_features(mgiid,ensembleid,genelength) values('%s','%s','%s')" % (mgiid,ensembleid,genelength)			
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))		
	print('insert gene length finished')	

def insert_exon_count(db):
	print('insert exon count')
	#file1 = '/home/david/Dropbox/datasets/essential genes prediction/test set/exon_count.csv'
	#file1 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_count.csv'
	#file1 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_count1013.csv'		
	file1 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_count10.csv'			
	cursor = db.cursor()
	i = 0
	exoncountL = [line.strip() for line in open(file1)]
	for line in exoncountL[1:len(exoncountL)]:
		i += 1
		vals = line.split(',')
		ensembleid = vals[0]
		exoncount = vals[1]
		try:
			sql = "update gene_features set exoncount = '%s' where ensembleid='%s'" % (exoncount,ensembleid)
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()
			sql = "insert into gene_features(ensembleid,exoncount) values('%s','%s')" % (ensembleid,exoncount)
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
	print('insert exon count finished')	

def insert_exon_length(db):
	print('insert exon length')
	#file = '/home/david/Dropbox/datasets/essential genes prediction/test set/exon_length.csv'
	#file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length.csv'
	#file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length1013.csv'	
	file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/exon_length10.csv'		
	cursor = db.cursor()
	i = 0
	exonlengthL = [line.strip() for line in open(file)]
	for line in exonlengthL[1:len(exonlengthL)]:
		i += 1
		vals = line.split(',')
		ensembleid = vals[0]
		exonlength = vals[1]
		try:
			sql = "update gene_features set exonlength = '%s' where ensembleid='%s'" % (exonlength,ensembleid)
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()
			sql = "insert into gene_features(ensembleid,exonlength) values('%s','%s')" % (ensembleid,exonlength)
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
	print('insert exon length finished')	

def insert_GC_content_and_transcript_count(db):
	print('insert GC content and transcript count')
	#file = '/home/david/Dropbox/datasets/essential genes prediction/test set/GC_transcript_count.csv'
	#file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/GC_transcript_count.csv'	
	file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/GC_transcript_count10.csv'	
	cursor = db.cursor()
	i = 0
	geneL = [line.strip() for line in open(file)]
	for line in geneL[1:len(geneL)]:
		i += 1
		vals = line.split(',')
		mgiid = vals[0]
		GCcontent = vals[1]
		transcriptcount = vals[2]
		try:
			sql = "update gene_features set GCcontent = '%s', transcriptcount = '%s' where mgiid = '%s'" % (GCcontent,transcriptcount,mgiid)
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()
			sql = "insert into gene_features(mgiid,GCcontent,transcriptcount) values('%s','%s','%s')" % (mgiid,GCcontent,transcriptcount)
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
	print('insert % GC content and transcript count finished')	

def insert_intron_length(db):
	print('insert intron length')
	#file = '/home/david/Dropbox/datasets/essential genes prediction/test set/intron_length.csv'
	#file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/intron_length.csv'	
	#file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/intron_length1013.csv'		
	file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/gene features/intron_length10.csv'		
	cursor = db.cursor()
	i = 0
	intronlengthL = [line.strip() for line in open(file)]
	for line in intronlengthL[1:len(intronlengthL)]:
		i += 1
		vals = line.split(',')
		mgiid = vals[0]
		ensembleid = vals[1]
		intronlength = vals[2]
		try:
			sql = "update gene_features set intronlength = '%s', ensembleid='%s' where mgiid = '%s'" % (intronlength,ensembleid,mgiid)
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()
			sql = "insert into gene_features(mgiid,ensembleid,intronlength) values('%s','%s','%s')" % (mgiid,ensembleid,intronlength)
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
	print('insert intron length finished')

def insert_pepstats_features(db):
	print('insert pepstats features')
	db.autocommit(False)
	#file = '/home/david/Dropbox/datasets/essential genes prediction/test set/pepstats.csv'
	#file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/unknowngenes.csv'	
	file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/pepstats/33proteins_not_in_fasta_file.csv'	
	cursor = db.cursor()
	i = 0
	pepstatsL = [line.strip() for line in open(file)]
	for line in pepstatsL[1:len(pepstatsL)]:
	#GeneName,UniProt_ID,MW,ProteinLength,IP,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Aliphatic,Aromatic,NonPolar,Polar,Charged,basic,Acidic
		i += 1
		vals = line.split(',')
		genename = vals[0]
		uniprotid = vals[1]
		MW = vals[2]
		proteinlength = vals[3]
		IP = vals[4]
		A = vals[5]
		C = vals[6]
		D = vals[7]
		E = vals[8]
		F = vals[9]
		G = vals[10]
		H = vals[11]
		I = vals[12]
		K = vals[13]
		L = vals[14]
		M = vals[15]
		N = vals[16]
		P = vals[17]
		Q = vals[18]
		R = vals[19]
		S = vals[20]
		T = vals[21]
		V = vals[22]
		W = vals[23]
		Y = vals[24]
		Aliphatic = vals[25]
		Aromatic = vals[26]
		NonPolar = vals[27]
		Polar = vals[28]
		Charged = vals[29]
		Basic = vals[30]
		Acidic = vals[31]
		try:
			sql = "INSERT INTO pepstats_features(genename,uniprotid,MW,ProteinLength,IP,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Aliphatic,Aromatic,NonPolar,Polar,Charged,Basic,Acidic) VALUES('%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s','%s','%s','%s', '%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s')" % (genename,uniprotid,MW,proteinlength,IP,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Aliphatic,Aromatic,NonPolar,Polar,Charged,Basic,Acidic)
			cursor.execute(sql)
			db.commit()
			print('insert '+str(i))
		except:
			db.rollback()
			sql = "update pepstats_features set MW = '%s', uniprotid='%s', ProteinLength = '%s', IP = '%s', A = '%s', C = '%s', D = '%s', E = '%s', F = '%s', G = '%s', H = '%s', I = '%s', K = '%s', L = '%s', M = '%s', N = '%s', P = '%s', Q = '%s', R = '%s', S = '%s', T = '%s', V = '%s', W = '%s', Y = '%s', Aliphatic = '%s', Aromatic = '%s',NonPolar = '%s',Polar = '%s',Charged = '%s',Basic = '%s', Acidic = '%s' where genename='%s'" % (MW,uniprotid,proteinlength,IP,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,Aliphatic,Aromatic,NonPolar,Polar,Charged,Basic,Acidic,genename)
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
	print('insert pepstats features finished')

def update_ids(db):
	print('update ids')
	file='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	db.autocommit(False)
	cursor = db.cursor()
	i = 0
	idsL = [line.strip() for line in open(file)]
	for line in idsL[1:len(idsL)]:
		i += 1
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		try:
			sql = "update gene_features set genename = '%s' where mgiid='%s'" % (genename,mgiid)
			cursor.execute(sql)
			db.commit()
			print('update '+str(i))
		except:
			db.rollback()
			print('rollback '+str(i))	
	print('update ids finished')

if __name__ == "__main__":
	main()
	
	


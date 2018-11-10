#retrieve data from MySQL database
import MySQLdb

def main():
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/notes/human genes/2096_human_lethal_genes_which_are_not_in_known_mouse_gene_set.csv'
	#mouse_viable_genes_file ='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/mouse_viable_genes.csv'
	#mouse_lethal_genes_file ='/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/mouse_lethal_genes.csv'
	#all_human_viable_genes_file = '/home/david/Dropbox/datasets/essential genes prediction/notes/human genes/all_human_viable_genes.csv'
	#human_lethal_genes_lek_paper = '/home/david/Dropbox/datasets/essential genes prediction/notes/human genes/human_lethal_genes_lek_paper.csv'
	#human_lethal_genes_blomen_paper = '/home/david/Dropbox/datasets/essential genes prediction/notes/human genes/human_lethal_genes_blomen_paper.csv'
	#human_lethal_genes_wang_paper = '/home/david/Dropbox/datasets/essential genes prediction/notes/human genes/human_lethal_genes_wang_paper.csv'
	#unknown_genes_ids = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv' 
	#all_genenames = '/home/david/Dropbox/datasets/essential genes prediction/all_lethal_viable_genenames.csv'
	#genes = '/home/david/Dropbox/datasets/essential genes prediction/genes_of_null_chromosome_in_db.csv'
	#unknown_genes = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes.csv' 
	db = MySQLdb.connect("localhost","david","david","mouse_genes")
	#select_ids_of_new_lethal_and_new_viable_genes_from_mouse_gene_table(db)
	#select_ids_of_transcript_per_million_table(db)
	#select_ids_of_gene_features_table(db)
	fw=open(outfile,'w')
	#fw=open(mouse_viable_genes_file,'w')
	#fw=open(mouse_lethal_genes_file,'w')
	#fw=open(all_human_viable_genes_file,'w')
	#fw=open(human_lethal_genes_lek_paper,'w')
	#fw=open(human_lethal_genes_blomen_paper,'w')
	#fw=open(human_lethal_genes_wang_paper,'w')
	#fw=open(unknown_genes_ids,'w')
	#fw=open(all_genenames,'w')
	#fw=open(unknown_genes,'w')
	#fw=open(genes,'w')
	#select_new_lethal_new_viable_genes(db,fw)#select features of new lethal new viable genes from db
	#select_human_lethal_genes_not_in_mouse_genes(db,fw)#select genenames of human lethal genes which are not in known mouse genes set
	#select_mouse_lethal_genes(db,fw)
	#select_mouse_viable_genes(db,fw)
	#select_human_viable_genes_from_all_papers(db,fw)
	#select_human_lethal_genes_from_Lek_paper(db,fw)
	#select_human_lethal_genes_from_Blomen_paper(db,fw)
	#select_human_lethal_genes_from_Wang_paper(db,fw)
	#select_ids_of_unknown_genes(db,fw)
	#select_all_genenames(db,fw)
	#select_genes_with_null_chromosome(db,fw)
	#select_unknown_genes(db,fw)
	db.close()

def select_genes_with_null_chromosome(db,fw):
	sql = "select genename,mgiid,ensembleid from mouse_gene where isnull(chromosome)"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			mgiid = row[1]
			ensembleid = row[2]
			fw.write(genename+','+mgiid+','+ensembleid+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_genes_with_null_chromosome finished')	

def select_all_genenames(db,fw):
	sql = "select genename from mouse_gene"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			fw.write(genename+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_all_genenames finished')
			
def select_ids_of_unknown_genes(db,fw):
	sql = "select genename,mgiid,uniprotid,ensembleid from mouse_gene where essentialitytype=\'predicted\'"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			uniprotid = row[1]
			mgiid = row[2]
			ensembleid = row[3]
			fw.write(genename+','+mgiid+','+uniprotid+','+ensembleid+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_ids_of_unknown_genes finished')		

def select_human_lethal_genes_from_Wang_paper(db,fw):
	sql = "select distinct genename from aac7041_SM_Table_S3 where essentiality=\'Lethal\'"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			fw.write(genename+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_human_lethal_genes_from_Wang_paper done')

def select_human_lethal_genes_from_Blomen_paper(db,fw):
	sql = "select distinct genename from human_gene_kbm7_table_s1 where essentiality=\'Lethal\' \
		union \
		select distinct genename from human_gene_hap1_table_s2 where essentiality=\'Lethal\' \
		order by genename"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			fw.write(genename+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_human_lethal_genes_from_Blomen_paper done')

def select_human_lethal_genes_from_Lek_paper(db,fw):
	sql = "select genename from human_gene where essentiality=\'Lethal\'"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			fw.write(genename+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_human_lethal_genes_from_Lek_paper done')

def select_human_viable_genes_from_all_papers(db,fw):
	sql="select distinct genename from human_gene where essentiality = \'Viable\' \
		union \
		select distinct genename from human_gene_kbm7_table_s1 where essentiality=\'Viable\' \
		union \
		select distinct genename from human_gene_hap1_table_s2 where essentiality=\'Viable\' \
		union \
		select distinct genename from aac7041_SM_Table_S3 where essentiality=\'Viable\' \
		order by genename"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			fw.write(genename+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_human_viable_genes_from_all_papers finished')

def select_mouse_lethal_genes(db,fw):
	sql = "select genename,mgiid,ensembleid,essentiality from mouse_gene where essentiality=\'Lethal\' and essentialitytype=\'known essentiality\'"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			mgiid = row[1]
			ensembleid = row[2]
			essentiality = row[3]
			fw.write(genename+','+mgiid+','+ensembleid+','+essentiality+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_mouse_lethal_genes finished')

def select_mouse_viable_genes(db,fw):
	sql = "select genename,mgiid,ensembleid,essentiality from mouse_gene where essentiality=\'Viable\' and essentialitytype=\'known essentiality\'"
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename = row[0]
			mgiid = row[1]
			ensembleid = row[2]
			essentiality = row[3]
			fw.write(genename+','+mgiid+','+ensembleid+','+essentiality+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_mouse_viable_genes finished')

def select_human_lethal_genes_not_in_mouse_genes(db,fw):
	sql = "select h.genename from human_gene as h where h.essentiality='Lethal' and h.genename not in (select genename from mouse_gene where essentialitytype = \'known essentiality\')"	
	cursor = db.cursor()
	try:	
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			genename=row[0]
			fw.write(genename+'\n')
		fw.close()
	except:
		print('sql query failed')
	print('select_human_genes finished')

def select_unknown_genes(db,fw):
	#GeneName	MGI_ID	UniProt_ID	Ensemble_ID	GeneLength	GC_content	Transcript_count	Exon_Count	ExonLength	IntronLength	MW	ProteinLength	Aliphatic	Aromatic	NonPolar	Polar	Charged	Basic	Acidic	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	Oocyte(Transcript/Million)	Unfertilized_Ovum(Transcript/Million)	Zygote(Transcript/Million)	Cleavage(Transcript/Million)	Morula(Transcript/Million)	Blastocyst(Transcript/Million)	Egg_Cylinder(Transcript/Million)	Gastrula(Transcript/Million)	Organogenesis(Transcript/Million)	Fetus(Transcript/Million)	Neonate(Transcript/Million)	Juveline(Transcript/Million)	Adult(Transcript/Million)	Age	Class
	sql = "select distinct m.genename as GeneName,m.mgiid as MGI_ID,m.uniprotid as Uniprot_ID,m.ensembleid as Ensemble_ID,gf.genelength as GeneLength,gf.GCcontent as GC_content,gf.transcriptcount as Transcript_count,gf.exoncount as Exon_Count,gf.exonlength as ExonLength,gf.intronlength as IntronLength,p.MW,p.ProteinLength,p.Aliphatic,p.Aromatic,p.NonPolar,p.Polar,p.Charged,p.Basic,p.Acidic,p.A,p.C,p.D,p.E,p.F,p.G,p.H,p.I,p.K,p.L,p.M,p.N,p.P,p.Q,p.R,p.S,p.T,p.V,p.W,p.Y,ec_kw.GlycoProtein,ec_kw.Phosphoprotein,ec_kw.Acetylation,ec_kw.Transcription,sig.SignalPeptide,ec_kw.ec1,ec_kw.ec2,ec_kw.ec3,ec_kw.ec4,ec_kw.ec5,ec_kw.ec6,sub.Nucleus_Score,sub.Cytoplasm_Score,sub.Plasma_Score,sub.Extracellular_Score,sub.Golgi_Score,sub.ER_Score,sub.Mitochondria_Score,sub.Peroxisome_Score,sub.Lysosome_Score,sub.Nucleus_UniProt,sub.Cytoplasm_UniProt,sub.Plasma_UniProt,sub.Membrane_UniProt,sub.Extracellular_Uniprot,sub.Mitochondrion_UniProt,sub.ER_UniProt,sub.Golgi_UniProt,sub.Lysosome_UniProt,sub.Peroxisome_UniProt,sub.CellJunction_UniProt,sub.CellProjection_UniProt,t.Transmembrane_Count,cyt.ASP_Known,cyt.BC_Known,cyt.ClosenessCentrality_Known,cyt.ClusteringCoefficient_Known,cyt.Degree_Known,cyt.TC_Known,h.BN_Score_Known,h.EPC_Score_Known,h.MNC_Score_Known,h.DMNC_Score_Known,cyt.ASP_KP,cyt.BC_KP,cyt.ClosenessCentrality_KP,cyt.ClusteringCoefficient_KP,cyt.Degree_KP,cyt.TC_KP,h.BN_Score_KP,h.EPC_Score_KP,h.MNC_Score_KP,h.DMNC_Score_KP,tpm.Oocyte as 'Oocyte(Transcript/Million)',tpm.Unfertilized_Ovum as 'Unfertilized_Ovum(Transcript/Million)',tpm.Zygote as 'Zygote(Transcript/Million)',tpm.Cleavage as 'Cleavage(Transcript/Million)',tpm.Morula as 'Morula(Transcript/Million)',tpm.Blastocyst as 'Blastocyst(Transcript/Million)',tpm.Egg_Cylinder as 'Egg_Cylinder(Transcript/Million)',tpm.Gastrula as 'Gastrula(Transcript/Million)',tpm.Organogenesis as 'Organogenesis(Transcript/Million)',tpm.Fetus as 'Fetus(Transcript/Million)',tpm.Neonate as 'Neonate(Transcript/Million)',tpm.Juveline as 'Juveline(Transcript/Million)',tpm.Adult as 'Adult(Transcript/Million)',m.essentiality as Class,m.essentialitytype \
from mouse_gene as m \
left join transcript_per_million as tpm \
on m.mgiid = tpm.mgiid \
left join gene_features as gf \
on m.mgiid = gf.mgiid \
left join pepstats_features as p \
on m.uniprotid = p.uniprotid \
left join EC_keywords_features as ec_kw \
on m.uniprotid = ec_kw.uniprotid \
left join cytoscape_ppi_features as cyt \
on m.uniprotid = cyt.uniprotid \
left join hubba_ppi_features as h \
on m.uniprotid = h.uniprotid \
left join signalpeptides as sig \
on m.uniprotid = sig.uniprotid \
left join subcellular_localizations as sub \
on m.uniprotid = sub.uniprotid \
left join transmembrane as t \
on m.uniprotid = t.uniprotid \
where m.essentialitytype=\'predicted\'"

	cursor = db.cursor()
	try:
		cursor.execute(sql)
		names=''
		for d in cursor.description:#column names
			names += ','+d[0]
		names = names.strip(',')
		fw.write(names+'\n')	
		rows = cursor.fetchall()
		n = 0
		for row in rows:
			n += 1
			gene=''
			i=0
			while i < len(row):
				gene += ','+row[i]
				i += 1
			gene = gene.strip(',')
			#print gene
			fw.write(gene+'\n')
		fw.close()
		print('no. of rows retrieved: '+str(n))
	except:
		print('sql query failed')
	print('select_new_lethal_new_viable_genes finished')

def select_new_lethal_new_viable_genes(db,fw):
#GeneName	MGI_ID	UniProt_ID	Ensemble_ID	GeneLength	GC_content	Transcript_count	Exon_Count	ExonLength	IntronLength	MW	ProteinLength	Aliphatic	Aromatic	NonPolar	Polar	Charged	Basic	Acidic	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	Oocyte(Transcript/Million)	Unfertilized_Ovum(Transcript/Million)	Zygote(Transcript/Million)	Cleavage(Transcript/Million)	Morula(Transcript/Million)	Blastocyst(Transcript/Million)	Egg_Cylinder(Transcript/Million)	Gastrula(Transcript/Million)	Organogenesis(Transcript/Million)	Fetus(Transcript/Million)	Neonate(Transcript/Million)	Juveline(Transcript/Million)	Adult(Transcript/Million)	Age	Class
	sql = 'select m.genename as GeneName,m.mgiid as MGI_ID,m.uniprotid as UniProt_ID,m.ensembleid as Ensemble_ID,gf.genelength as GeneLength,gf.GCcontent as GC_content,gf.transcriptcount as Transcript_count,gf.exoncount as Exon_Count,gf.exonlength as ExonLength,gf.intronlength as IntronLength,p.MW,p.ProteinLength,p.Aliphatic,p.Aromatic,p.NonPolar,p.Polar,p.Charged,p.Basic,p.Acidic,p.A,p.C,p.D,p.E,p.F,p.G,p.H,p.I,p.K,p.L,p.M,p.N,p.P,p.Q,p.R,p.S,p.T,p.V,p.W,p.Y,tpm.Oocyte as \'Oocyte(Transcript/Million)\',tpm.Unfertilized_Ovum as \'Unfertilized_Ovum(Transcript/Million)\',tpm.Zygote as \'Zygote(Transcript/Million)\',tpm.Cleavage as \'Cleavage(Transcript/Million)\',tpm.Morula as \'Morula(Transcript/Million)\',tpm.Blastocyst as \'Blastocyst(Transcript/Million)\',tpm.Egg_Cylinder as \'Egg_Cylinder(Transcript/Million)\',tpm.Gastrula as \'Gastrula(Transcript/Million)\',tpm.Organogenesis as \'Organogenesis(Transcript/Million)\',tpm.Fetus as \'Fetus(Transcript/Million)\',tpm.Neonate as \'Neonate(Transcript/Million)\',tpm.Juveline as \'Juveline(Transcript/Million)\',tpm.Adult as \'Adult(Transcript/Million)\',m.essentiality as Class from mouse_gene as m inner join transcript_per_million as tpm on m.genename = tpm.genename inner join gene_features as gf on tpm.mgiid = gf.mgiid inner join pepstats_features as p on m.uniprotid = p.uniprotid where m.essentialitytype=\'known essentiality since May 2015\''
	cursor = db.cursor()
	try:
		cursor.execute(sql)
		names=''
		for d in cursor.description:
			names += ','+d[0]
		names = names.strip(',')
		fw.write(names+'\n')	
		rows = cursor.fetchall()
		n = 0
		for row in rows:
			n += 1
			gene=''
			i=0
			while i < len(row):
				gene += ','+row[i]
				i += 1
			gene = gene.strip(',')
			#print gene
			fw.write(gene+'\n')
		fw.close()
		print('no. of rows retrieved: '+str(n))
	except:
		print('sql query failed')
	print('select_new_lethal_new_viable_genes finished')

def select_ids_of_new_lethal_and_new_viable_genes_from_mouse_gene_table(db):
	sql = 'select genename, mgiid, uniprotid, ensembleid, unigeneid, essentiality from mouse_gene where essentialitytype = \'known essentiality since May 2015\''
	cursor = db.cursor()
	try:
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			#print row
			print row[0]+','+row[1]+','+row[2]+','+row[3]+','+row[4]
	except:
		print('sql query failed')
	print('select_ids_of_new_lethal_and_new_viable_genes_from_mouse_gene_table finished')

def select_ids_of_transcript_per_million_table(db,idfile):	
	#select the ids of transcript_per_million table
	sql = 'select genename, mgiid, uniprotid, ensembleid, unigeneid from transcript_per_million'
	cursor = db.cursor()
	try:
		cursor.execute(sql)
		rows = cursor.fetchall()
		for row in rows:
			#print row
			print row[0]+','+row[1]+','+row[2]+','+row[3]+','+row[4]
	except:
		print('sql query failed')
	print('\nselect_ids_of_transcript_per_million_table finished')

if __name__ == "__main__":
	main()
	

#insert human gene essentiality data into MySQL database
import MySQLdb
import re

def main():
	#genesfile='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes_ensembleids.csv' #viable genes from MacArthur paper and Sulem et al paper
	#genesfile2='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes2_ensembleids.csv' #viable genes from "Homozygous loss-of-function variants in European cosmopolitan and isolate populations" paper by Kaiser
	#genesfile3='/home/david/Dropbox/datasets/essential genes prediction/human essential genes/human_essential_genes_ensembleids.csv'# lethal genes from paper "Analysis of protein-coding genetic variation in 60,706 humans"
	#genesfile4='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes.csv' #viable genes from MacArthur paper and Sulem et al paper
	#genesfile5='/home/david/Dropbox/datasets/essential genes prediction/human viable genes/human_viable_genes2.csv' #viable genes from "Homozygous loss-of-function variants in European cosmopolitan and isolate populations" paper
	#genesfile6='/home/david/Dropbox/datasets/essential genes prediction/human essential genes/human_essential_genes.csv'# lethal genes from paper "Analysis of protein-coding genetic variation in 60,706 humans"
	#genestable1='/home/david/Dropbox/papers/gene prediction/human genetics paper/gene essentiality and synthetic lethality in haploid human cells/aac7557_SM_Table_S1.csv'#lethal and viable genes from paper "gene essentiality and synthetic lethality in haploid human cells"
	#genestable2='/home/david/Dropbox/papers/gene prediction/human genetics paper/gene essentiality and synthetic lethality in haploid human cells/aac7557_SM_Table_S2.csv'#lethal and viable genes from paper "gene essentiality and synthetic lethality in haploid human cells"
	genestableS3='/home/david/Dropbox/papers/gene prediction/human genetics paper/Identification and characterization of essential genes in the human genome/aac7041_SM_Table_S3.csv'#lethal and viable genes from paper "Identification and characterization of essential genes in the human genome" by Wang et al
	#genesL = [line.strip() for line in open(genesfile)]
	#genes2L = [line.strip() for line in open(genesfile2)]
	#genes3L = [line.strip() for line in open(genesfile3)]
	#genes4L = [line.strip() for line in open(genesfile4)]
	#genes5L = [line.strip() for line in open(genesfile5)]
	#genes6L = [line.strip() for line in open(genesfile6)]
	#genestable1L = [line.strip() for line in open(genestable1)]
	#genestable2L = [line.strip() for line in open(genestable2)]
	genestableS3L = [line.strip() for line in open(genestableS3)]
	db = MySQLdb.connect("localhost","david","david","mouse_genes")
	#insert_genes(genesL,'Viable',db)
	#insert_genes(genes2L,'Viable',db)
	#insert_genes(genes3L,'Lethal',db)
	#human_viable_genes_ensembleids.csv, human_viable_genes2_ensembleids.csv, human_essential_genes_ensembleids.csv contain ensemble ids of some genenames only. So insert the genenames with no ensemble ids into tables
	#insert_genes2(genes4L,'Viable',db)
	#insert_genes2(genes5L,'Viable',db)
	#insert_genes2(genes6L,'Lethal',db)
	#table1 = 'human_gene_kbm7_table_s1'
	#insert_genes3(genestable1L,table1,'Lethal',db)
	#table2 = 'human_gene_hap1_table_s2'
	#insert_genes3(genestable2L,table2,'Lethal',db)
	tableS3 = 'aac7041_SM_Table_S3'
	insert_genes4(genestableS3L,tableS3,db)
	db.close()

def insert_genes(genesL,essentiality,db):
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:	
		vals = line.split(',')
		genename = vals[0]
		ensembleid = vals[1]
		sql = "INSERT INTO human_gene(genename,ensembleid,essentiality)\
	 		VALUES('%s', '%s', '%s')" % \
       			(genename,ensembleid,essentiality)
		try:
			cursor.execute(sql)
			db.commit()
			print('commit insert the gene at line'+str(i))
		except:
   			#Rollback in case there is any error e.g. the table contains the genename (primary key) to insert already
			db.rollback()
			print('rollback at line '+str(i))
		i += 1
	print('insert_genes done')

def insert_genes2(genesL,essentiality,db):
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL:	
		genename = line
		sql = "insert into human_gene(genename,essentiality)\
			values('%s','%s')" % \
			(genename,essentiality)
		try:
			cursor.execute(sql)
			db.commit()
			print('commit insert the gene at line'+str(i))
		except:#if the gene to insert is already in the table and has essentiality Viable, update its essentiality to Lethal (viable genes can become lethal genes later after experimentation)
			sql = "UPDATE human_gene SET essentiality='Lethal' where genename='%s' and essentiality='Viable'" % (genename)
			cursor.execute(sql)
			db.commit()
			print('update essentiality to Lethal at line '+str(i))
			'''
			db.rollback()
			#print('rollback at line '+str(i))
			'''
		i += 1
	print('insert_genes2 done')

def insert_genes3(genesL,table,essentiality,db):
	cursor = db.cursor()
	i=1
	for line in genesL:	
		vals = line.split(',')
		genename = vals[0]
		ensembleid = vals[1]
		selected = vals[len(vals)-1]
		if selected == 'YES':
			essentiality='Lethal'
		else:
			essentiality='Viable'
		sql = "insert into "+table+"(genename,ensembleid,essentiality) values('%s','%s','%s')" % (genename,ensembleid,essentiality)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert the gene at line'+str(i))
		except:
   			db.rollback()#end the current transaction if sql statement fails before start a new transaction
			sql = "update "+table+" set essentiality=\'"+essentiality+"\', ensembleid='%s' where genename='%s'" % (ensembleid,genename)
			cursor.execute(sql)
			db.commit()
			print('update line '+str(i))
		i+=1
	print('insert_genes3 finished')

def insert_genes4(genesL,table,db):
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:
		vals = line.split(',')
		genename = vals[0]
		KBM7_CS = vals[2]
		KBM7_adjusted_pvalue = vals[3]
		if float(KBM7_CS) < -0.1 and float(KBM7_adjusted_pvalue) < 0.05:#"Genes with a CS < -0.1 and corrected p< 0.05 in the KBM7 cell line were defined as cell-essential in downstream analyses." p8 of "Supplementary Material for Identification and characterization of essential genes in the human genome" by Wang
			essentiality='Lethal'
		else:
			essentiality='Viable'
		sql = "insert into "+table+"(genename,essentiality) values('%s','%s')" % (genename,essentiality)
		try:
			cursor.execute(sql)
			db.commit()
			print('insert the gene at line'+str(i))
		except:
   			db.rollback()#end the current transaction if sql statement fails before start a new transaction
			print('rollback at line: '+str(i)+' genename: '+genename)
		i+=1
	print('insert_genes4 done')

if __name__ == "__main__":
	main()
	
	


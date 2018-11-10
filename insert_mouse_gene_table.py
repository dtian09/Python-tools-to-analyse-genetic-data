#insert data into mouse_gene table in mouse_genes db
import MySQLdb
import re
#1. store the known genes into database
#2. store the unknown genes predictions into database
#3. store the new lethal and new viable genes into database (new lethal and new viable genes became known essentiality genes after unknown genes were collected so the new lethal and new viable genes are also in the unknown genes file)
#4. store the ids of the unknown genes into database
#
#database transactions are used to insert the data into database so that if a gene record exists in the database and it is inserted again, the record is skipped without stop running the program with duplicate primary key errors. 
 
def main():
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/known essentiality genes/training_all_genesinfo.csv'
	genesfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/predictions/prediction_reduced_48_features_data_with_gene_names.random_forest.space'
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/new_lethal_new_viable_genes_not_in_train_set_missing_vals.csv'
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids_unigeneids.csv'
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/all_lethal_viable_chromosome_names_startloc_endloc.csv'
	#genesfile ='/home/david/Dropbox/datasets/essential genes prediction/chromosome_names_startloc_endloc6_2.csv'
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/mgiid_ensembleid12.csv'
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/ids2.txt'
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ids10.csv'
	#genesfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/unigeneids_none.csv'
	genesL = [line.strip() for line in open(genesfile)]	
	db = MySQLdb.connect("localhost","david","david","mouse_genes")
	#update_unigeneiids(genesL,db)
	#update_ensembleids(genesL,db)
	#insert_gene_startloc_endloc_on_chromosome(genesL,db)
	#update_gene_startloc_endloc_on_chromosome(genesL,db)
	#update_mgiids(genesL,db)
	#insert_known_genes(genesL,db)
	insert_predictions(genesL,db)
	#insert_new_lethal_new_viable(genesL,db)
	#insert_ids_of_unknown_genes(genesL,db)
	db.close()

def update_unigeneiids(genesL,db):
	print('update unigeneiids')
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for genename in genesL[1:len(genesL)]:
		sql = "update mouse_gene set unigeneid = 'none' where genename = '%s'" % (genename)
		try:#start a transaction
			cursor.execute(sql)
			db.commit()#commit the transaction and end this transaction
			print('update unigene id at line '+str(i))
		except:
			print('update unigene id fails at line '+str(i)+' '+genename)
			db.rollback()#end this transaction
		i+=1

def update_mgiids(genesL,db):
	print('update mgiids')
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL:
		vals = line.split('\t')
		genename = vals[0]
		mgiid = vals[2]
		sql = "update mouse_gene set mgiid = '%s' where genename = '%s'" % (mgiid,genename)
		try:#start a transaction
			cursor.execute(sql)
			db.commit()#commit the transaction and end this transaction
			print('update gene at line '+str(i))
		except:
			print('update gene fails at line '+str(i)+' '+genename)
			db.rollback()#end this transaction
		i+=1

def update_ensembleids(genesL,db):
	print('update ensembleids')
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:
		vals = line.split(',')
		mgiid = vals[0]
		ensembleid = vals[1]
		sql = "update mouse_gene set ensembleid = '%s' where mgiid = '%s'" % (ensembleid,mgiid)
		try:#start a transaction
			cursor.execute(sql)
			db.commit()#commit the transaction and end this transaction
			print('update gene at line '+str(i))
		except:
			print('update gene fails at line '+str(i)+' '+genename)
			db.rollback()#end this transaction
		i+=1

def update_gene_startloc_endloc_on_chromosome(genesL,db):
	print('update_gene_startloc_endloc_on_chromosome')
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		ensembleid = vals[2]
		chromosome = vals[3]
		gene_start = vals[4]
		gene_end = vals[5]
		sql = "update mouse_gene set chromosome = '%s', startloc = '%s', endloc = '%s' where ensembleid= '%s'" % (chromosome,gene_start,gene_end,ensembleid)
		try:#start a transaction
			cursor.execute(sql)
			db.commit()#commit the transaction and end this transaction
			print('update gene at line '+str(i))
		except:
			print('update gene fails at line '+str(i)+' '+genename)
			db.rollback()#end this transaction
			try:#start a transaction
				sql = "INSERT INTO mouse_gene(genename,mgiid,ensembleid,chromosome,startloc,endloc) VALUES ('%s','%s','%s','%s', '%s', '%s')" % \
					(genename,mgiid,ensembleid,chromosome,gene_start,gene_end)
				cursor.execute(sql)
				db.commit()
				print('insert gene succeeds')
   			except:
				#Rollback in case there is any error e.g. the table contains the genename (primary key) to insert already
				db.rollback()#end this transaction
				print('rollback at line '+str(i))
		i += 1

def insert_gene_startloc_endloc_on_chromosome(genesL,db):
	print('insert_gene_startloc_endloc_on_chromosome')
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:	
		vals = line.split(',')
		genename = vals[0]
		chromosome = vals[1]
		gene_start = vals[2]
		gene_end = vals[3]
		sql = "update mouse_gene set chromosome = '%s', startloc = '%s', endloc = '%s' where genename = '%s'" % (chromosome,gene_start,gene_end,genename)
		try:#start a transaction
			cursor.execute(sql)
			db.commit()#commit the transaction and end this transaction
			print('commit at line '+str(i))
		except:
			print('update genename, chromosome, gene start, gene end fails at line '+str(i))
			db.rollback()#end this transaction
			try:#start a transaction
				sql = "INSERT INTO mouse_gene(genename,chromosome,startloc,endloc) VALUES ('%s', '%s', '%s', '%s')" % \
					(genename,chromosome,gene_start,gene_end)
				cursor.execute(sql)
				db.commit()
				print('insert genename, chromosome, gene start, gene end succeeds')
   			except:
				#Rollback in case there is any error e.g. the table contains the genename (primary key) to insert already
				db.rollback()#end this transaction
				print('rollback at line '+str(i))
		i += 1

def insert_known_genes(genesL,db):
	print('insert_known_genes')
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:	
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		uniprotid = vals[2]
		ensembleid = vals[3]
		essentiality = vals[len(vals)-1]
		sql = "INSERT INTO mouse_gene(genename,mgiid,ensembleid,uniprotid,essentiality,essentialitytype)\
     			VALUES ('%s', '%s', '%s', '%s', '%s', '%s')" % \
       			(genename,mgiid,ensembleid,uniprotid,essentiality,'known essentiality')
		try:
			cursor.execute(sql)
			db.commit()
			print('commit at line '+str(i))
		except:
   			#Rollback in case there is any error e.g. the table contains the genename (primary key) to insert already
			db.rollback()
			print('rollback at line '+str(i))
		i += 1

def insert_predictions(genesL,db):
	print('insert predictions')
	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:	
		#vals = line.split('\t')
		vals = line.split(' ')
		genename = vals[0]
		pred = vals[1]
		conf = vals[2]
		#1:Lethal
		m = re.match('[12]+\:([LethalVib]+)',pred)
		if m:
			pred = m.group(1)
			sql = "UPDATE mouse_gene set essentiality = '%s', essentialitytype='predicted', confidence = '%s' WHERE genename = '%s'" % (pred,conf,genename)
			try:
				cursor.execute(sql)
				db.commit()
				print('Update predictions of unknown genes at line'+str(i))
			except:
				db.rollback()
				sql = "INSERT INTO mouse_gene(genename,essentiality,essentialitytype,confidence)\
				VALUES ('%s', '%s', '%s', '%f')" % \
				(genename,pred,'predicted',conf)
				print('Insert predictions of unknown genes at line '+str(i))
		i += 1
	print('insert predictions finished')

def insert_new_lethal_new_viable(genesL,db):
	print('insert_new_lethal_new_viable')
	#update the essentialities and essentialitytypes of new lethal and new viable genes in the unknown gene set to their essentialities obtained in May 2015 
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		uniprotid = vals[2]
		ensembleid = vals[3]
		essentiality = vals[len(vals)-1]		
		try:
			sql = "UPDATE mouse_gene SET mgiid='%s', ensembleid='%s', uniprotid='%s', essentiality='%s',essentialitytype='known essentiality since May 2015', confidence=Null where genename ='%s'" % (mgiid,ensembleid,uniprotid,essentiality,genename)
			cursor.execute(sql)
		except:
			#Rollback in case there is any error e.g. the table does not contain the record to update
			db.rollback()#end this transaction
			print('update ids and essentialities fails at line '+str(i))
			try:
				db.autocommit(False)
				sql = "INSERT INTO mouse_gene(genename,mgiid,ensembleid,uniprotid,essentiality,essentialitytype)\
     				VALUES ('%s', '%s', '%s', '%s', '%s', '%s')" % \
       				(genename,mgiid,ensembleid,uniprotid,essentiality,'known essentiality since May 2015')
				cursor.execute(sql)
				db.commit()
   			except:
				#Rollback in case there is any error e.g. the table contains the genename (primary key) to insert already
				db.rollback()
				print('rollback insert ids and essentialities of new lethal and new viable genes at line '+str(i))
		i += 1

def insert_ids_of_unknown_genes(genesL,db):
	print('insert_ids_of_unknown_genes')
	#Update the ids of unknown genes in the gene table with the ids of the unknown genes in genenames_mgiids_uniprotids_ensembleids_unigeneids.csv
 	db.autocommit(False)
	cursor = db.cursor()
	i=1
	for line in genesL[1:len(genesL)]:
		vals = line.split(',')
		genename = vals[0]
		mgiid = vals[1]
		uniprotid = vals[2]
		ensembleid = vals[3]
		unigeneid = vals[4]
		try:
			#sql = "UPDATE gene SET mgiid='%s', ensembleid='%s', uniprotid='%s' where genename='%s' and (isnull(mgiid) or isnull(ensembleid) or isnull(uniprotid))" % (mgiid,ensembleid,uniprotid,genename)
			sql = "UPDATE mouse_gene SET mgiid='%s', ensembleid='%s', uniprotid='%s', unigeneid='%s' where genename='%s'" % (mgiid,ensembleid,uniprotid,unigeneid,genename)
			cursor.execute(sql)
			db.commit()
			print('commit update mgi ids, ensemble ids, uniprot ids of unknown genes at line '+str(i))
		except:
			#Rollback in case there is any error e.g. the table does not contain the record to update
			db.rollback()
			#print('rollback update mgi ids, ensemble ids, uniprot ids of unknown genes at line '+str(i))				
		i += 1

if __name__ == "__main__":
	main()


#create a sql query condition

def main():
	#infile='/home/david/Dropbox/datasets/essential genes prediction/mouse_genes_and_non_mouse_genes'
	infile='/home/david/Dropbox/datasets/essential genes prediction/MGIalleleQuery_20160316_125157.txt'
	infileL = [line.strip() for line in open(infile)]
	sql_cond=''
	n = 0
	'''
	for line in infileL[1:len(infileL)]:
		valsL = line.split('\t')
		genename = valsL[0]
		sql_cond += "'"+genename+"',"
		n += 1 
	'''
	for line in infileL[1:len(infileL)]:
		valsL = line.split('\t')
		mgiid = valsL[0]
		sql_cond += "'"+mgiid+"',"
		n += 1
	print("no. of genenames: "+str(n))
	sql_cond = sql_cond.strip(',')
	sql_cond = '('+sql_cond+')'
	print(sql_cond)
	
if __name__ == "__main__":
        main()


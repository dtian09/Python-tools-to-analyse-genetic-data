#create a Unigene search query for Unigene ids of genes
import sys

def main():
	infile='/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_withno_unigeneids'
	infileL = [line.strip() for line in open(infile)]
	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	query=''
	#no_of_unigeneids = 0
	for line in infileL[2:len(infileL)]:
		#valsL = line.split('\t')
		#unigeneid = valsL[1]
		#query += ' OR '+unigeneid
		#no_of_unigeneids += 1
		query +=' OR '+line
	print(query)
	#print('no. of unigene ids: '+str(no_of_unigeneids))

if __name__ == "__main__":
        main()


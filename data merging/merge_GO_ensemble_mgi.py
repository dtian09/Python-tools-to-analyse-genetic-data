#Merge GO annotations collected from ensemble and MGI databases
#
#input: output file of collectfeatures_GO_ensemble.py (a data file containing GO annotations in database e.g. ensemble)
#	output file of collectfeatures_GO_mgi.py (a data file containing GO annotations in another database e.g. mgi)
#output: the merged data file containing all GO annotations of ensemble and mgi
#
#format of input data files:
#MGI_id,term1,term2,...,termN,class
#mgi:1,0,1,0,...,1,?
#mgi:2,1,0,1,...,0,?
#
#the 1st line contains "MGI_id", GO terms names and "class"
#
#author: David Tian	tiand@cs.man.ac.uk

import sys
import re
import collectfeatures_GO

def main():
	in_file = sys.argv[1]
	in_file2 = sys.argv[2]
	outfile = sys.argv[3]

	genesL = [line.strip() for line in open(in_file)]
	
	if len(genesL)==0:
		print(in_file+" is empty.")
		sys.exit(-1)

	genesL2 = [line.strip() for line in open(in_file2)]
	
	if len(genesL2)==0:
		print(in_file2+" is empty.")
		sys.exit(-1)

	#get all GO annotations of 1st data file
	h1 = genesL[0]#get the 1st line of file
	m = re.match("[\w:]+,(.+),class",h1)
	if m:
		GOannots1 = m.group(1)
		GOannotsAll1 = GOannots1.split(',')
	else:
		print(h1+" does not match pattern 1")
		sys.exit(-1)
	#get all GO annotations of 2nd data file
	h2 = genesL2[0]#get the 1st line of file
	m = re.match("[\w:]+,(.+),class",h2)
	if m:
		GOannots2 = m.group(1)
		GOannotsAll2 = GOannots2.split(',')
	else:
		print(h2+" does not match pattern 2")
		sys.exit(-1)
		
	ht = create_hashtable(genesL[1:len(genesL)],GOannotsAll1)#hashtable: key=mgi id, value = set of GO annots of a gene
	ht2 = create_hashtable(genesL2[1:len(genesL2)],GOannotsAll2)#hashtable: key=mgi id, value = set of GO annots of a gene
	GOannotsAll1set = set(GOannotsAll1)
	GOannotsAll2set = set(GOannotsAll2)
	GOannotsAllset = GOannotsAll1set.union(GOannotsAll2set)
	GOannotsAllL = list(GOannotsAllset)
	
	fw = open(outfile,"w")
	fw.write("MGI_Gene_ID,")
	for GOannot in GOannotsAllL:
		fw.write(GOannot+",")
	fw.write("class\n")
	collectfeatures_GO.write_GOannotations3(ht,ht2,GOannotsAllL,fw,'?')
	fw.close()
	outfileL = [line.strip() for line in open(outfile)]
	print("no. of genes in the file "+outfile+": "+str(len(outfileL)-1))
	print("no. of GO annotations of genes in "+in_file+": "+str(len(GOannotsAll1)))
	print("no. of GO annotations of genes in "+in_file2+": "+str(len(GOannotsAll2)))
	print("no. of GO annotations of genes in "+outfile+": "+str(len(GOannotsAllL)))
	
def create_hashtable(genesL,GOannotsAllL):
	#input: a list with each element a line in a data file
	#	a list of all GO annotations in the data file
	#output: a hashtable: key=mgi id, value = set of GO annots of a gene
	#	 a set of all GO annotations of the genes
	ht={}#hashtable: key=mgi id, value = set of GO annots of a gene
	for line in genesL:
		GOannots = set()
		m = re.match("([\w:]+),([01,]+),?",line)
		if m:
			key=m.group(1)
			vals = m.group(2)
			vals = vals.rstrip(',')#remove any ending ","
			valsL = vals.split(',')
			i=0
			for val in valsL:
				if val == '1':
					GOannots.add(GOannotsAllL[i])
				i+=1
			if ht.get(key)== None:
				ht[key] = GOannots
			else:
				GOannots2 = ht[key]
				ht[key] = GOannots2.union(GOannots)				
		else:
			print(line+" does not match pattern 3")
	return ht

if __name__ == "__main__":
	main()


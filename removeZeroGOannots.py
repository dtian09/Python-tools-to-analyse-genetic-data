#Remove the GO annotations which have 0 frequencies
#input: a data file with GO annotations and biological process, cellular component and molecular function, class as columns
#output: the data file with zero frequency GO annotations removed
#foramt of data file:
#Ensemble_Gene_ID,GOterm1,GOterm2,...,GOtermN,biological_process,cellular_component,molecular_function,class
#ENSMUSG00000020154,1,1,0,0,0,1,1,0,1,1,...,1

import sys

def main():
	infile = sys.argv[1]
	outfile = sys.argv[2]
	fw = open(outfile,'w')
 	
	data = [line.strip() for line in open(infile)]
	if len(data)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	line = data[0]
	cols = line.split(',')
	del cols[0] #delete the ensemble id
	features = len(cols)-4
	bp = len(cols)-4
	cc = len(cols)-3
	mf = len(cols)-2
	c = len(cols)-1
	
	frequencies = []
	i=0
	while i < features:
		frequencies.append(0)
		i += 1
	
	for line in data[1:len(data)]:
		vals = line.split(',')
		del vals[0]#remove ensemble id
		i=0
		while i < features:
			frequencies[i] += int(vals[i])
			i += 1
	line = data[0]
	cols = line.split(',')
	ensemble_id = cols[0]
	del cols[0]
	fw.write('Ensemble_Gene_ID,')
	i=0
	while i < features:
		if frequencies[i] > 0:
			fw.write(cols[i]+',')
		i += 1
	fw.write("biological_process,cellular_component,molecular_function,class\n")
		
	for line in data[1:len(data)]:
		vals = line.split(',')
		fw.write(vals[0]+',')
		del vals[0] #remove ensemble id
		i = 0
		while i < features:
			if frequencies[i] > 0:
				fw.write(vals[i]+',')
			i += 1
		fw.write(vals[bp]+','+vals[cc]+','+vals[mf]+','+vals[c]+"\n")
	fw.close()

if __name__ == "__main__":
	main()

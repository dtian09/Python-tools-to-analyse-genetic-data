#replace values of genes with other values

import sys

def main():
	infile='/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_ids.csv'
	outfile='/home/david/Dropbox/datasets/essential genes prediction/new lethal genes/new_lethal_genes_ids2.csv'

	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)
	fw = open(outfile,'w')

	hash_new_values = {}
	hash_new_values['Rpa1']='Mm.180734'
	hash_new_values['Ddost']='Mm.7236'
	hash_new_values['Eprs']='Mm.154511'
	hash_new_values['Gtf2a1']='Mm.275728'
	hash_new_values['Zfp207']='Mm.401332'
	hash_new_values['Vps13d']='Mm.240310'
	hash_new_values['Mapkap1']='Mm.270866'
	hash_new_values['Cdc26']='Mm.109930'
	hash_new_values['Patz1']='Mm.275563'
	hash_new_values['Ints2']='Mm.440906'
	hash_new_values['L3mbtl2']='Mm.491152'
	hash_new_values['Sgpl1']='Mm.412319'
	hash_new_values['Tulp3']='Mm.291783'
	hash_new_values['1700067K01Rik']='Mm.492236'
	hash_new_values['Ssbp1']='Mm.276356'
	hash_new_values['Gpr107']='Mm.236074'
	n = 0
	f_indx = 4
	fw.write(infileL[0]+'\n')
	for line in infileL[1:len(infileL)]:
		vals = line.split(',')
		genename = vals[0]
		f_val = vals[f_indx]
		if hash_new_values.get(genename)!= None:
			fw.write(vals[0]+','+vals[1]+','+vals[2]+','+vals[3]+','+hash_new_values[genename]+'\n')
			print('gene name: '+vals[0]+', old value: '+f_val+', new value: '+hash_new_values[genename])
			n += 1
		else:
			fw.write(line+'\n')
	fw.close()
	print('no. of values changed: '+str(n))

if __name__ == "__main__":
	main()

	

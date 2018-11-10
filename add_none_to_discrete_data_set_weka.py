#program to add 'none' to PPI network features of balanced1307_discretized_sig_level=0_95.arff
import re

def main():
	#discrete_data_file = '/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_discretized_sig_level=0_95.arff'
	#numeric_data_file = '/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_missing_values.arff'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_discretized_sig_level=0_95_none_values2.arff'
	#discrete_data_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_discretized_by_cuts_sig_level=0_95.arff'
	#numeric_data_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_missing_values.arff'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_discretized_by_cuts_sig_level=0_95_none_values2.arff'
	#discrete_data_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes_discretized_by_cuts_sig_level=0_95.arff'
	#numeric_data_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes_missing_values.arff'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/103 features data/unknowngenes_discretized_by_cuts_sig_level=0_95_none_values2.arff'
	discrete_data_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_viable_genes_test_set2_discretized_by_cuts_sig_level=0_95.arff'
	numeric_data_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_viable_genes_test_set2_missing_vals.arff'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_viable_genes_test_set2_discretized_by_cuts_sig_level=0_95_none_values2.arff'

	discrete_dataL = [line.strip() for line in open(discrete_data_file)]
	numeric_dataL = [line.strip() for line in open(numeric_data_file)]
	add_none_values(discrete_dataL,numeric_dataL,outfile)

def add_none_values(discrete_dataL,numeric_dataL,outfile):
	#Add 'none' values to the PPI network features definitions from ASP_Known to DMNC_Score_KP
	#Replace the missing values of the PPI network features with 'none'
	fw = open(outfile,"w")
	i=0
	attr_indx=1
	while (i < len(discrete_dataL)):
		line = discrete_dataL[i]
		line2 = numeric_dataL[i]
		m = re.match("^@attribute\s+Class\s+\{[LethalVib\,]+\}$",line)#match class attribute
		m2 = re.match("^(@attribute\s+[^\{]+\s+\{[^\{]+)\}$",line)  #@attribute GC_content {'\'(-inf-32.07)\'','\'[32.07-35.02)\''}
		m3 = re.match("^[?\-\d\.\,LethalVib]+$",line2)#match an instance in continuous data file containing missing values ?
		if m:#match class attribute, write it to file
			fw.write(discrete_dataL[i]+"\n")
		elif m2:
			#print(str(attr_indx)+': '+line)	
			if attr_indx >= 69 and attr_indx <= 88:#match an attribute definition, add 'none' to the definition
				attr_defn = m2.group(1)+",none}"
			else:
				attr_defn = line
			attr_indx += 1
			fw.write(attr_defn+'\n')
		elif m3:
			l = line.split(',')
			l2 = line2.split(',')
			new_instance = l[0]
			col=1
			while col <= 67:
				new_instance = new_instance+','+l[col]
				col += 1
			col=68
			while col <= 87:
				if l2[col] == '?':
					new_instance = new_instance+',none'
				else:
					new_instance = new_instance+','+l[col]
				col+=1
			col=88
			while col <=101:
				new_instance = new_instance+','+l[col]
				col+=1
	 		class_label = l[102]
			fw.write(new_instance+','+class_label+'\n')
		else:
			fw.write(line+'\n')
		i+=1
	fw.close()
'''
		elif m3:#match an instance in continuous data file, replace '?' with 'none' in discrete data file
			l = line.split(',')
			l2 = line2.split(',')
			if l2[0] == '?':
				new_instance = 'none'
			else:
				new_instance = l[0]
			j=1
			for val in l2[1:len(l2)-1]:
				if val == '?':
					new_instance = new_instance+',none'
				else:
					new_instance = new_instance+','+l[j]
				j+=1
			class_label = l[len(l)-1]
			fw.write(new_instance+','+class_label+'\n')	
		else:
			fw.write(discrete_dataL[i]+"\n")
		i+=1
'''
	
if __name__ == "__main__":
	main()
			

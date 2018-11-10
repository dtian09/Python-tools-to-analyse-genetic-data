#program to add 'none' to PPI network features of balanced1307_discretized_sig_level=0_95.arff
def main():
	#discrete_data_file = '/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_discretized_sig_level=0_95.csv'
	#numeric_data_file = '/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_missing_values.csv'
	#outfile = '/home/david/Dropbox/datasets/essential genes prediction/train set/103 features data/balanced1307_discretized_sig_level=0_95_none_values2.csv'
	discrete_data_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_discretized_by_cuts_sig_level=0_95.csv'
	numeric_data_file = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_missing_values.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/test set/103 features data/BlindTestSet1_19Aug15_discretized_by_cuts_sig_level=0_95_none_values2.csv'
	discrete_dataL = [line.strip() for line in open(discrete_data_file)]
	numeric_dataL = [line.strip() for line in open(numeric_data_file)]
	add_none_values(discrete_dataL,numeric_dataL,outfile)

def add_none_values(discrete_dataL,numeric_dataL,outfile):
	#Add 'none' values to the PPI network features definitions from ASP_Known to DMNC_Score_KP
	#Replace the missing values of the PPI network features with 'none'
	fw = open(outfile,"w")
	fw.write(discrete_dataL[0]+'\n')
	i=1
	while (i < len(discrete_dataL)):
		line = discrete_dataL[i]
		line2 = numeric_dataL[i]
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
		i+=1
	fw.close()
	
if __name__ == "__main__":
	main()
			

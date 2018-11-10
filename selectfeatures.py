#select features which have frequencies >= a cut off threshold
#
#input: a boolean data file
#	has_id_col (0,1)
#	has_class_attr (0,1)
#	a cut off threshold e.g. 100
#output: a boolean data file containing only the features with frequencies equal to or above the cut off threshold
#
#data file format:
#id,f1,f2,f3,...,class
#MGI:1,0,1,0,...,?
#MGI:2,1,1,0,...,?
#MGI:3,0,1,0,...,?
#MGI:4,0,0,0,...,?
#MGI:5,0,1,0,...,?
#MGI:6,1,1,1,...,?
#
#or
#
#id,f1,f2,f3,...fn
#MGI:1,0,1,0,...,1
#MGI:2,1,1,0,...,0
#MGI:3,0,1,0,...,0
#MGI:4,0,0,0,...,1
#MGI:5,0,1,0,...,1
#MGI:6,1,1,1,...,0
#
import sys

def main():
	data_file = sys.argv[1]
	has_id_col = sys.argv[2]#first column is ids (0 or 1) 
	has_class_attr= sys.argv[3]#data file has class attribute (0 or 1)
	threshold = sys.argv[4]
	out_file = sys.argv[5]

	if has_class_attr != '1' and has_class_attr != '0':
		print("has_class_attr should be 0 or 1")
		sys.exit(-1)
		
	dataL = [line.strip() for line in open(data_file)]
	
	if len(dataL)==0:
		print(data_file+" is empty.")
		sys.exit(-1)

	if has_class_attr == '1':
		no_of_features = len(dataL)-2
	else:
		no_of_features = len(dataL)-1
			
	if int(threshold)>= no_of_features:
		print("cut off threshold >= total number of features of data file!")
		print("total no. of features: "+no_of_features+", cut off threshold: "+threshold)
		sys.exit(-1)
		
	(freq,featuresL) = count_freq_of_features(dataL,has_id_col,has_class_attr)
	print("list of all features: "+str(featuresL))
	no_of_selected_features = select_features(featuresL,freq,threshold,dataL,has_id_col,has_class_attr,out_file)
	
	print("no. of features selected: "+no_of_selected_features)
	
def count_freq_of_features(dataL,has_id_col,has_class_attr):
	#input: a list with each element a line of the data file
	#	has id col (0 or 1)
	#	has_class_attr (0 or 1)
	#output: hash table (key=feature, value=frequency of feature)
	#initialize a hash table
	freq = {}
	features = dataL[0]
	featuresL = features.split(",")
	if has_class_attr == '0':#no class attribute
		if has_id_col == '0':#no id field
			featuresL2 = featuresL
		else:
			featuresL2 = featuresL[1:len(featuresL)]#skip the id field
		for f in featuresL2:
			freq[f]=0
	else:#has class attribute
		if has_id_col == '0':#no id field
			featuresL2 = featuresL
		else:
			featuresL2 = featuresL[1:len(featuresL)-1]#skip the id field and the class attribute
		for f in featuresL2:
			freq[f]=0
	#insert frequencies into hash table
	if has_class_attr == '0':#no class attribute
		for line in dataL[1:len(dataL)]:#skip the first line
			valsL = line.split(",")
			i=0
			if has_id_col == '0':#no id field
				valsL2 = valsL
			else:
				valsL2 = valsL[1:len(valsL)]#skip the id 
			for f in featuresL2:
				freq[f] += int(valsL2[i])
				i+=1
	else:#has class attribute
		for line in dataL[1:len(dataL)]:#skip the first line
			line = line.rstrip(",") #remove the last ','
			valsL = line.split(",")
			i=0
			if has_id_col == '0':#no id field
				valsL2 = valsL[0:len(valsL)-1]#skip the class attribute
			else:
				valsL2 = valsL[1:len(valsL)-1]#skip the id and the class attribute
			for f in featuresL2:
				freq[f] += int(valsL2[i])
				i+=1
	return (freq,featuresL2)		

def select_features(featuresL,freq,threshold,dataL,has_id_col,has_class_attr,out_file):
	#input: a list of all features of the data file 
	#	hash table (key=feature, value = frequency of feature)
	#	cut off threshold
	#	a list with each element a line in the data file
	# 	has_id_col (0 or 1)
	#	has_class_attr (0 or 1)
	#output: output data file
	no_of_selected_features=0
	
	fw=open(out_file,"w")
	#write the first line of the reduced data file
	line1 = dataL[0]
	line1L = line1.split(",")
	if has_id_col == '1':
		fw.write(line1L[0]+",")#write the id field
		if has_class_attr == '0':#no class attribute
			for f in featuresL[0:len(featuresL)-1]:#skip the last feature
				if freq[f]>=int(threshold):	
					fw.write(f+",")
					no_of_selected_features +=1
			last_f = featuresL[len(featuresL)-1]
			if freq[last_f] >= int(threshold):#write last feature
				fw.write(last_f)
				no_of_selected_features +=1
			fw.write("\n")
		else:#has class attribute
			for f in featuresL:
				if freq[f]>=int(threshold):	
					fw.write(f+",")
					no_of_selected_features +=1
			fw.write(line1L[len(line1L)-1])#write class attribute
			fw.write("\n")
	else:#no id col
		if has_class_attr == '0':#no class attribute
			for f in featuresL[0:len(featuresL)-1]:#skip the last feature
				if freq[f]>=int(threshold):	
					fw.write(f+",")
					no_of_selected_features +=1
			last_f = featuresL[len(featuresL)-1]
			if freq[last_f] >= int(threshold):#write last feature
				fw.write(last_f)
				no_of_selected_features +=1
			fw.write("\n")
		else:#has class attribute
			for f in featuresL:
				if freq[f]>=int(threshold):	
					fw.write(f+",")
					no_of_selected_features +=1
			fw.write(line1L[len(line1L)-1])#write class attribute
			fw.write("\n")
	#write each line of the reduced data file
	if has_id_col == '1':
		if has_class_attr == '0':#no class attribute
			for line in dataL[1:len(dataL)]:#skip first line containing the features names
				valsL = line.split(",")
				fw.write(valsL[0]+",")#write id field
				i=0 #feature index
				for val in valsL[1:len(valsL)-1]:#skip the id and the last feature
					f = featuresL[i]
					if freq[f] >= int(threshold):
						fw.write(val+",")
					i+=1
				last_f = featuresL[len(featuresL)-1]
				if freq[last_f] >= int(threshold):#write last feature
					fw.write(valsL[len(valsL)-1])
				fw.write("\n")
		else:#has class attribute
			for line in dataL[1:len(dataL)]:#skip first line containing the features names
				valsL = line.split(",")
				fw.write(valsL[0]+",")#write id field
				i=0 #feature index
				for val in valsL[1:len(valsL)-1]:#skip the id and the class (last column)
					f = featuresL[i]
					if freq[f] >= int(threshold):
						fw.write(val+",")
					i+=1
				fw.write(valsL[len(valsL)-1])#write class
				fw.write("\n")
	else:#no id col
		if has_class_attr == '0':#no class attribute
			for line in dataL[1:len(dataL)]:#skip first line containing the features names
				valsL = line.split(",")
				i=0 #feature index
				for val in valsL[0:len(valsL)-1]:#skip the last feature
					f = featuresL[i]
					if freq[f] >= int(threshold):
						fw.write(val+",")
					i+=1
				last_f = featuresL[len(featuresL)-1]
				if freq[last_f] >= int(threshold):#write last feature
					fw.write(valsL[len(valsL)-1])
				fw.write("\n")
		else:#has class attribute
			for line in dataL[1:len(dataL)]:#skip first line containing the features names
				valsL = line.split(",")
				i=0 #feature index
				for val in valsL[0:len(valsL)-1]:#skip the class (last column)
					f = featuresL[i]
					if freq[f] >= int(threshold):
						fw.write(val+",")
					i+=1
				fw.write(valsL[len(valsL)-1])#write class
				fw.write("\n")

	fw.close()
	return str(no_of_selected_features)

if __name__ == "__main__":
	main()

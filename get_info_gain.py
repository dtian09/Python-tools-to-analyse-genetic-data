#get information gain of all features and a feature subset to plot a boxplot in Matlab
#input: information gain of all features output by weka
#	the list of selected features e.g. features in a random forest
#output: a vector of information gain of all the features
#	 a vector of the information gain of the selected features
#	 a vector of 1s of same length as no. of selected features
#e.g.
# information_gain_of_all_features = [0.0001,0.01020,0.35000,0.73000,0.0100,0.0020,0.5000,0.3000,0.420100,0.20020,0.75000,0.003000]
# information_gain_of_selected_features = [0.0100,0.0020,0.5000,0.3000]
# boxplot(information_gain_of_all_features)
# hold on
# plot([1,1,1,1],information_gain_of_selected_features,'r+')

import sys
import re

def main():
	infile='/home/david/Dropbox/datasets/essential genes prediction/train set/93 features data/info_gain_all_features'
	infile2='/home/david/Dropbox/datasets/essential genes prediction/train set/93 features data/random_forest_selected_features'
	
	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	infile2L = [line.strip() for line in open(infile2)]

	if len(infile2L)==0:
		print(infile2+" is empty.")
		sys.exit(-1)

	info_gain_hash={}
	info_gain_of_all_features=[]
	info_gain_of_selected_features=[]
	onesL0=[]
	onesL=[]
	all_features=set()
	selected_features=set()
	for line in infileL:
		m = re.match('^\s*(\d+.*\d*)\s+\d+\s+([\w\(\)/.-]+)\s*$',line)
		if m:
			info_gain = m.group(1)
			feature = m.group(2)
			all_features.add(feature)
			info_gain_hash[feature]=float(info_gain)
			info_gain_of_all_features.append(float(info_gain))
			onesL0.append(1)
		else:
			print(line+' does not match pattern')
	print('info_gain_hash: '+str(len(info_gain_hash)))
	for feature in infile2L:
		feature = feature.strip()
		selected_features.add(feature)
		if info_gain_hash.get(feature)!=None:
			info_gain = info_gain_hash[feature]
			info_gain_of_selected_features.append(float(info_gain))
			onesL.append(1)
		else:
			print(feature+' is not a key in info_gain_hash')
	
	print('info_gain_of_all_features=')
	info_gain_of_all_features.sort()
	print(info_gain_of_all_features)
	print('info_gain_of_selected_features=')
	info_gain_of_selected_features.sort()
	print(info_gain_of_selected_features)
	print('onesL0=')
	print(onesL0)
	print('onesL=')
	print(onesL)
	print('no. of features in info_gain_of_all_features: '+str(len(info_gain_of_all_features)))
	print('no. of features in info_gain_of_selected_features: '+str(len(info_gain_of_selected_features)))
	removed_features = all_features.difference(selected_features)
	print('features removed: '+str(removed_features))

if __name__ == "__main__":
	main()


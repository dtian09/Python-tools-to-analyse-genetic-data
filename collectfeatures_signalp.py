#Collect signal peptide feature from Uniprot and signal peptide prediction of signalP tool
#
#input: a tab file containing signal peptide feature of proteins from Uniprot
#	a text file containing signal peptide predictions of the proteins from signalP tool
#output: a csv file containing signal peptide information of the proteins in the Unprot tab file and the signal peptide predictions of the proteins whose signal peptide information is missing in the tab file
#
#tab file format:
#Entry	Signal peptide
#Q8CDK1	
#E9PWT2	
#Q3UJR4	
#Q8JZZ0	SIGNAL 1 22 {ECO:0000255}.
#Q9DAP7	
#
#SignalP file format:
#
# SignalP-4.1 euk predictions
# name                     Cmax  pos  Ymax  pos  Smax  pos  Smean   D     ?  Dmaxcut    Networks-used
#tr|Q8CDK1|Q8CDK1_MOUSE     0.111  57  0.118  12  0.154   1  0.115   0.117 N  0.450      SignalP-noTM
#tr|E9PWT2|E9PWT2_MOUSE     0.112  24  0.151  11  0.247   9  0.190   0.172 N  0.450      SignalP-noTM
#tr|Q3UJR4|Q3UJR4_MOUSE     0.731  31  0.784  31  0.926  14  0.838   0.813 Y  0.450      SignalP-noTM
#
import re
import sys
import collectfeatures_ppi

def main():
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'		
	tab_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/signal peptide/uniprot_signal_peptide.tab'
	signalp_file1 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/signal peptide/10000genes.signalp'
	signalp_file2 = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/signal peptide/5495genes.signalp'
	out_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/signal peptide/uniprot_signalp.csv'
	genesidL = [line.strip() for line in open(genesid_file)]
	tab_fileL = [line.strip() for line in open(tab_file)]
	signalp_file1L = [line.strip() for line in open(signalp_file1)]
	signalp_file2L = [line.strip() for line in open(signalp_file2)]
	uniprotids = collectfeatures_ppi.get_uniprotids(genesidL)
	#uniprot_hash: key=uniprot_id, value=1 or ?(missing signalp)
	uniprot_hash = get_uniprot_signalpeptide(tab_fileL[1:len(tab_fileL)])#skip 1st line
	#signalp_hash: key=uniprot_id, value=1 or 0
	signalp_hash = {}
	signalp_hash = get_signalp_prediction(signalp_hash,signalp_file1L[2:len(signalp_file1L)])#skip 1st and 2nd lines
	signalp_hash = get_signalp_prediction(signalp_hash,signalp_file2L[2:len(signalp_file2L)])#skip 1st and 2nd lines
	fw = open(out_file,"w")
	printdata(uniprotids,uniprot_hash,signalp_hash,fw)
	fw.close()
	
def get_uniprot_signalpeptide(tab_fileL):
	#signalp_uniprot_hash: key=uniprot_id, value=1 or ?(missing signalp)
	uniprot_hash = {}
	for line in tab_fileL:#get uniprot id of each longest protein
		m = re.match('^([\w]+)\s*$',line)#missing signal peptide information
		if m:
			uniprot_id = m.group(1)
			uniprot_hash[uniprot_id]='?'
		else:
			m2 = re.match('^([\w]+)\s+["SIGNAL]+.+$',line)
			if m2:
				uniprot_id = m2.group(1)
				uniprot_hash[uniprot_id]='1'
			else:
				print(line+' does not match pattern m2 in get_uniprot_signalp')
	return uniprot_hash
				
def get_signalp_prediction(signalp_hash,signalp_fileL):
	for line in signalp_fileL:
		cols = line.split('|')
		uniprot_id = cols[1]
		info = cols[2]
		m = re.match('^.+\s+([YN])\s*.+$',info)
		if m:
			pred = m.group(1)
			if pred == 'Y':
				signalp_hash[uniprot_id] = '1'
			else:
				signalp_hash[uniprot_id] = '0'
		else:
			print(info+' does not match pattern m in get_predicted_signalp')
	return signalp_hash
				
def printdata(uniprotids,uniprot_hash,signalp_hash,fw):
	fw.write('UniProt_ID,SignalPeptide\n')
	for uniprotid in list(uniprotids):
		fw.write(uniprotid+',')
		if uniprot_hash.get(uniprotid)!=None:
			if uniprot_hash[uniprotid] == '?':#uniprot does not have signal peptide information for the protein
				if signalp_hash.get(uniprotid) != None:
					fw.write(signalp_hash[uniprotid]+'\n')#write predicted value of signalP tool
				else:#the protein is not predicted by signalP tool 
					fw.write('?\n')#write a missing value
			else:
				fw.write(uniprot_hash[uniprotid]+'\n')#write the value from Uniprot
		else:#Uniprot does not have the Uniprot id
			if signalp_hash.get(uniprotid)!= None:
				fw.write(signalp_hash[uniprotid]+'\n')#write predicted value of signalP
			else:#the protein is not predicted by signalP 
				fw.write('?\n')#write a missing value
			
if __name__ == "__main__":
	main()
							

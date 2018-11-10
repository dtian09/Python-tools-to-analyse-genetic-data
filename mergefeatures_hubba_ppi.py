#program to merge Hubba PPI features
#input: BN_Score_Known.htm
#	BN_Score_KP.htm
#	Degree_KP.htm
#	DMNC_Score_Known.htm
#	DMNC_Score_KP.htm
#	EPC_Score_Known.htm
#	EPC_Score_KP.htm
#	MNC_Score_Known.htm
#	MNC_Score_KP.htm
#output: merged_hubba_features.csv
#format of .htm files
#
#<center><br><h3>Detail Results of BN Method</h3>
#<TABLE border="1">
#<TR bgcolor="lightblue"><TD><B>Name</B></TD><TD><B>Score</B></TD><TD><B>Rank</B></TD></TR>
#<TR><TD>P20263</TD><TD align="right">567.00000</TD><TD align="right">1</TD></TR>
#<TR bgcolor="lightyellow"><TD>P59240</TD><TD align="right">124.00000</TD><TD align="right">2</TD></TR>
#<TR><TD>P58252</TD><TD align="right">99.00000</TD><TD align="right">3</TD></TR>
#<TR bgcolor="lightyellow"><TD>P50516</TD><TD align="right">69.00000</TD><TD align="right">4</TD></TR>
#<TR><TD>Q62425</TD><TD align="right">63.00000</TD><TD align="right">5</TD></TR>
#<TR bgcolor="lightyellow"><TD>P84228</TD><TD align="right">51.00000</TD><TD align="right">6</TD></TR>
#<TR><TD>P70670</TD><TD align="right">46.00000</TD><TD align="right">7</TD></TR>
#</TABLE>
#<form METHOD="post"><input TYPE="button" VALUE="Go Back" OnClick="history.go( -1 );return true;"></form></center>
import collectfeatures_ppi
import mergefeatures_cytoscape_ppi
import re

def main():
	BN_score_known = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/BN_Score_Known.htm'
	BN_score_kp = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/BN_Score_KP.htm'
	DMNC_score_known = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/DMNC_Score_Known.htm'
	DMNC_score_kp = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/DMNC_Score_KP.htm'
	EPC_score_known = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/EPC_Score_Known.htm'
	EPC_score_kp = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/EPC_Score_KP.htm'
	MNC_score_known = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/MNC_Score_Known.htm'
	MNC_score_kp = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/MNC_Score_KP.htm'
	geneids_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/ppi network features/hubba/ppi_hubba_features.csv'
	BN_score_knownL = [line.strip() for line in open(BN_score_known)]
	BN_score_kpL = [line.strip() for line in open(BN_score_kp)]
	DMNC_score_knownL = [line.strip() for line in open(DMNC_score_known)]
	DMNC_score_kpL = [line.strip() for line in open(DMNC_score_kp)]
	EPC_score_knownL = [line.strip() for line in open(EPC_score_known)]
	EPC_score_kpL = [line.strip() for line in open(EPC_score_kp)]
	MNC_score_knownL = [line.strip() for line in open(MNC_score_known)]
	MNC_score_kpL = [line.strip() for line in open(MNC_score_kp)]
	geneidsL = [line.strip() for line in open(geneids_file)]
	fw = open(outfile,'w')
	BN_known = read_into_hashtable(BN_score_knownL)
	BN_kp = read_into_hashtable(BN_score_kpL)
	DMNC_known = read_into_hashtable(DMNC_score_knownL)
	DMNC_kp = read_into_hashtable(DMNC_score_kpL)
	EPC_known = read_into_hashtable(EPC_score_knownL)
	EPC_kp = read_into_hashtable(EPC_score_kpL)
	MNC_known = read_into_hashtable(MNC_score_knownL)
	MNC_kp = read_into_hashtable(MNC_score_kpL)
	genenames_hash = mergefeatures_cytoscape_ppi.get_genenames(geneidsL)
	uniprotids = collectfeatures_ppi.get_retrieved_longest_reviewed_proteins(geneidsL)
	merge_features(genenames_hash,uniprotids,BN_known,BN_kp,DMNC_known,DMNC_kp,EPC_known,EPC_kp,MNC_known,MNC_kp,fw)

def merge_features(genenames_hash,uniprotids,BN_known,BN_kp,DMNC_known,DMNC_kp,EPC_known,EPC_kp,MNC_known,MNC_kp,fw):
	fw.write('GeneName,UniProt_ID,BN_Score_Known,EPC_Score_Known,MNC_Score_Known,DMNC_Score_Known,BN_Score_KP,EPC_Score_KP,MNC_Score_KP,DMNC_Score_KP\n');
	for k in list(uniprotids):
		fw.write(genenames_hash[k]+','+k+',')
		if BN_known.get(k)!= None:
			fw.write(BN_known[k]+',')
		else:
			fw.write('?,')
		if EPC_known.get(k) != None:
			fw.write(EPC_known[k]+',')
		else:
			fw.write('?,')
		if MNC_known.get(k) != None:
			fw.write(MNC_known[k]+',')
		else:
			fw.write('?,')
		if DMNC_known.get(k) != None:
			fw.write(DMNC_known[k]+',')
		else:
			fw.write('?,')
		if BN_kp.get(k) != None:
			fw.write(BN_kp[k]+',')
		else:
			fw.write('?,')
		if EPC_kp.get(k) != None:
			fw.write(EPC_kp[k]+',')
		else:
			fw.write('?,')
		if MNC_kp.get(k) != None:
			fw.write(MNC_kp[k]+',')
		else:
			fw.write('?,')
		if DMNC_kp.get(k) != None:
			fw.write(MNC_kp[k]+'\n')
		else:
			fw.write('?\n')
	fw.close()
	
def read_into_hashtable(infileL):
	hashtable = {}
	for line in infileL:
		m = re.match("<tr><td>([\w]+)</td><td align=\"right\">([\d.]+)</td><td align=\"right\">[\d]+</td></tr>",line)
		if m:
			uniprot_id = m.group(1)
			score = m.group(2)
			hashtable[uniprot_id]=score
		else:
			m2 = re.match("<tr bgcolor=\"lightyellow\"><td>([\w]+)</td><td align=\"right\">([\d.]+)</td><td align=\"right\">[\d]+</td></tr>",line)
			if m2:
				uniprot_id = m2.group(1)
				score = m2.group(2)
				hashtable[uniprot_id]=score
	return hashtable

if __name__ == "__main__":
	main()





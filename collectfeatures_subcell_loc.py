#collect subcellular locations from uniprot and Wolf psort predicted locations
#
#input: a ids file containing uniprot ids of the proteins whose subcellular locations are to be collected 
#	an Uniprot xml file containing subcellular locations of proteins
#	Wolfpsort prediction file
#output: a data file containing subcell locations (0 or 1) and Wolf psort prediction scores in csv format
#
#notes:
#The subcellular locations of the proteins in the tab file are obtained from from the xml file.
#The xml file contains the longest proteins of the unknown genes and may also contain other poteins.
#The tab file contains the longest proteins only. 
#<comment type="subcellular location"> holds subcellular locations of a protein. If a protein does not have subcellular locations, there is no <comment type="subcellular location"> 
#
#xml file format:
#<uniprot>
#<entry>
#<accession>P897KD1</accession>
#<accession>M601WP1</accession>
#<accession>R056DL2</accession>
#</entry>
#<entry>
#<accession>Q897KD1</accession>
#<accession>D190DL21</accession>
#<comment type="function">
#<text evidence="1">Ligand for members of the frizzled family of seven transmembrane receptors.
#</text>
#</comment>
#<comment type="subcellular location">
#<subcellularLocation>
#<location evidence="1">Postsynaptic cell membrane</location>
#<location evidence="1">Cell membrane</location>
#</subcellularLocation>
#<subcellularLocation>
#<location evidence="1">Secreted</location>
#<location evidence="1">Extracellular space</location>
#<location evidence="1">Extracellular matrix</location>
#</subcellularLocation>
#<subcellularLocation>
#<location evidence="1">Midbody</location>
#</subcellularLocation>
#</comment>
#</entry>
#</uniprot>
#
#tab file format:
#Entry	Subcellular location [CC]
#Q8C6P4	SUBCELLULAR LOCATION: Secreted, extracellular space, extracellular matrix {ECO:0000256|RuleBase:RU003500}.
#
#Wolf psort file format
#
# k used for kNN is: 32
#tr|Q8CDK1|Q8CDK1_MOUSE nucl 23.5, cyto_nucl 16, cyto 7.5
#tr|E9PWT2|E9PWT2_MOUSE nucl 30.5, cyto_nucl 16.5
#tr|Q3UJR4|Q3UJR4_MOUSE extr 31
#tr|Q3UWG2|Q3UWG2_MOUSE cyto 11.5, cyto_mito 7, extr_plas 5.5, plas 5, extr 4, nucl 4, E.R. 3, pero 2, E.R._golg 2

from xml.dom.minidom import parse
import xml.dom.minidom
import sys
import re
import collectfeatures_ppi

def main():
	genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/genenames_mgiids_uniprotids_ensembleids.csv'
	#genesid_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/4ids.csv'
	xml_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/subcell loc/subcell_loc.xml'
	#xml_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/subcell loc/4.xml'
	psort_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/subcell loc/subcell_loc.wolfpsort'
	#psort_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/subcell loc/4.wolfpsort'
	out_file = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/subcell loc/uniprot_wolfpsort_loc.csv'
	genesidL = [line.strip() for line in open(genesid_file)]
	uniprotids = collectfeatures_ppi.get_uniprotids(genesidL)
	uniprot_locs_hash = {} #hash table: key=uniprot id, value = set of subcell locations of protein in Uniprot
	uniprot_locs_hash = get_subcell_locs(uniprotids,xml_file,uniprot_locs_hash)
	print(str(len(uniprot_locs_hash)))
 	#count no. of subcell locations of the longest proteins
	all_locs_uniprot = get_all_uniprot_subcell_locs(uniprot_locs_hash)
	#psort_locs_hash: key=uniprot id, value = list of (subcell location,score) tuples
	(psort_locs_hash,all_locs_psort) = wolf_psort_locs(uniprotids,psort_file)
	print(str(len(psort_locs_hash)))
	print('no. of single subcellular location scores output by Wolf psort: '+str(len(all_locs_psort)))
	#write subcell locations to a csv file
	fw = open(out_file,"w")
	printdata(uniprotids,all_locs_uniprot,uniprot_locs_hash,all_locs_psort,psort_locs_hash,fw)
	fw.close()

def get_all_uniprot_subcell_locs(subcell_locs_hash):
	uniprot_ids = list(subcell_locs_hash.keys())
	all_locs_uniprot = set()
	for uniprot_id in uniprot_ids:
		all_locs_uniprot = all_locs_uniprot.union(subcell_locs_hash[uniprot_id])
	print('no. of subcell locations from Uniprot: '+str(len(all_locs_uniprot)))
	return all_locs_uniprot

def wolf_psort_locs(uniprotids,psort_file):
	psort_fileL = [line.strip() for line in open(psort_file)]
	if len(psort_fileL)==0:
		print(psort_fileL+" is empty.")
		sys.exit(-1)	
	all_psort_locs = set()	
	locs_scores_hash={} #locs hash table: key=uniprot id, value = [(location 1, score of location 1),(location 2, score of location 2),...,(location k, score of location k)]
	for line in psort_fileL[1:len(psort_fileL)]:#skip the first line i.e. # k used for KNN is: 32
		cols = line.split('|')
		uniprotid = cols[1]
		locs_scores = cols[2]
		m = re.match('^[\w]+_MOUSE(.+)$',locs_scores)
		if m:
			locs_scores = m.group(1)
			locs_scoresL = locs_scores.split(',')
		else:
			print(locs_scores+' does not match 1st subcell location pattern.')
			sys.exit(-1)
		if uniprotid in uniprotids:
			for loc_score in locs_scoresL:
				loc_score = loc_score.strip()#remove leading and ending white spaces			
				m = re.match('^([^_]+)\s+([\d\.]+)$',loc_score)#loc_score represents a single subcell location and its score e.g. nucl 16 i.e. protein is at one nucl only.
				if m:					      
					loc = m.group(1)
					score = m.group(2)
					(locs_scores_hash,all_psort_locs) = add_locs(uniprotid,loc,score,all_psort_locs,locs_scores_hash)
				else:
					print(loc_score+' does not match the single subcell location pattern and not retrieved.')
	return (locs_scores_hash,all_psort_locs)

def add_locs(uniprotid,loc,score,all_psort_locs,locs_scores_hash):
	if locs_scores_hash.get(uniprotid) == None:
		locs_scores_hash[uniprotid] = [(loc,score)]
	else:
		locsL = locs_scores_hash[uniprotid]
		locsL.append((loc,score))
		locs_scores_hash[uniprotid] = locsL
	all_psort_locs.add(loc)
	return (locs_scores_hash,all_psort_locs)

def get_subcell_locs(uniprotids,xml_file,subcell_locs_hash):
	#input: uniprot ids of proteins whose subcellular locations are to be collected
	#	an uniprot xml file
	#output: hash table of (uniprot id, subcell locations) pairs
	#subcell_locs_hash = {} #hash table of the proteins in xml file: key=uniprot id, value = set of subcell locations of protein
	parser = xml.dom.minidom.parse(xml_file)
	uniprot = parser.documentElement
	entrysL = uniprot.getElementsByTagName("entry")
	for entry in entrysL:#an entry is a record in the xml file
		accessionsL = entry.getElementsByTagName("accession")#accession is an uniprot id
		commentsL = entry.getElementsByTagName("comment")
		if commentsL==[]:#protein has no subcellular location
			subcell_locs = set()
		else:
			subcell_locs = set()
			for comment in commentsL:
				comment_type = comment.getAttribute("type")
				if comment_type == "subcellular location":#if a protein has no subcell location information in Uniprot, the xml file does not have <subcellularLocation> element for the protein
					subcell_locsL = comment.getElementsByTagName("subcellularLocation")
					for subcell_loc in subcell_locsL:
						locsL = subcell_loc.getElementsByTagName("location")#locsL is a list of <location> elements
						for loc in locsL: #loc is an <location> element
							subcell_locs.add(loc.childNodes[0].data)								
		for accession in accessionsL:#get an accession i.e. an uniprot id
			uniprot_id = accession.childNodes[0].data
			if uniprot_id in uniprotids:
				subcell_locs_hash[uniprot_id] = subcell_locs
 	return subcell_locs_hash

def printdata(uniprotids,all_locs_uniprot,uniprot_locs_hash,all_locs_psort,psort_locs_hash,fw):
	#format:
	#uniprot_id,uniprot_loc1,uniprot_loc2,...,uniprot_locK,psort_loc1,psort_loc2,...,psort_loc9
	#Q819DK,1,1,...,0,4.5,8.9,...,6.0	
	fw.write('UniProt_ID,')
	all_locs_uniprotL = list(all_locs_uniprot)
	all_locs_psortL = list(all_locs_psort)
	#print(all_locs_psortL)
	for loc in all_locs_uniprotL:
		fw.write(loc+',')
	for loc in all_locs_psortL:		
		fw.write(loc+',')
	fw.write('\n')
	for uniprotid in uniprotids:
		fw.write(uniprotid+',')
		#print uniprot subcell locations
		if uniprot_locs_hash.get(uniprotid) != None:
			locsSet = uniprot_locs_hash[uniprotid]
			for loc in all_locs_uniprotL:
				if loc in locsSet:
					fw.write('1,')
				else:
					fw.write('0,')	
		else:#Uniprot does not have subcell location information of the protein i.e. the protein is not located at any of the subcell locations
			for loc in all_locs_uniprotL:
				fw.write('0,')
		#print wolf psort subcell locations
		if psort_locs_hash.get(uniprotid) != None:
			locs_scoresL = psort_locs_hash[uniprotid]
			for loc2 in all_locs_psortL:
				protein_has_score_at_loc2 = False					
				for (loc,score) in locs_scoresL:#search for this location loc2 in the wolfpsort prediction file
					if loc == loc2:#loc2 is in the wolfpsort prediction file. Write its score to file
						protein_has_score_at_loc2 = True
						fw.write(score+',')
						break
				if protein_has_score_at_loc2 == False:#wolfpsort does not predict a score for any location, every location has 0 score		
					fw.write('0,')
		else:#wolfpsort does not predict a score for any location, every location has 0 score
			for loc2 in all_locs_psortL:
				fw.write('0,')
		fw.write('\n')	
				 	
if __name__ == "__main__":
	main()


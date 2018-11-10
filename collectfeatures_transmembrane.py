#program to collect transmembrane count (no. of transmembrane segments of a protein) of unknown genes
#A transmembrane protein spans the whole membrane from one side of it to the other side of it.
#membrane proteins are targets of over 50% of the drugs on the market.
# 
#input: transmembrane_count-uniprot.tab (downloaded from Uniprot)
#output: data file
#format of transmembrane_count-uniprot.tab
#
#Entry	Transmembrane
#Q8CDK1	
#E9PWT2	
#Q3UJR4	
#Q3UWG2	
#Q8C6P4	
#Q9D563	TRANSMEM 7 27 Helical. {ECO:0000255}.; TRANSMEM 105 125 Helical. {ECO:0000255}.; TRANSMEM 133 153 Helical. {ECO:0000255}.; TRANSMEM 188 208 Helical. {ECO:0000255}.
#Q9Z110	
#Q8JZZ0	TRANSMEM 488 508 Helical. {ECO:0000255}.
#Q9DAP7

import sys
import re

def main():
	#infile = sys.argv[1]
	#outfile = sys.argv[2]
	infile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/transmembrane count/transmembrane_count.tab'
	outfile = '/home/david/Dropbox/datasets/essential genes prediction/unknown essentiality genes/transmembrane count/transmembrane_count.csv'
	infileL = [line.strip() for line in open(infile)]	
	fw = open(outfile,'w')
	fw.write('uniprot_id,transmembrane_count\n')
	for line in infileL[1:len(infileL)]:
		count_segments=0;
		m = re.match('^([\w]+)\s*$',line)#line does not contain transmembrane segments	
		if m:
			fw.write(m.group(1)+',0'+'\n')
		else:	
			m2 = re.match('^([\w]+)\s+(.+)$',line)#line contains 1 or more transmembrane segments
			if m2:
				uniprot_id = m2.group(1)
				segments = m2.group(2)
				l = segments.split('.;')
				for segment in l:
					m3 = re.match('^\s*TRANSMEM.+$',segment)
					if m3 != None:
						#print(segment+' match pattern')
						count_segments +=1
					#else:
					#	print(segment+' does not match pattern')
				fw.write(uniprot_id+','+str(count_segments)+'\n')
	fw.close()

if __name__ == "__main__":
	main()										

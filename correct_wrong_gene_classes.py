#correct wrong genes classes in known genes data set i.e. training set
#
#31 genes are Viable in newly known gene set, but they are Lethal in training set.
#
#the 31 genes:
#Syt1,Rest,Flii,Lama4,Atp2b1,Pcsk1,Il6st,Inpp5e,Sirt1,Pdpn,Traf6,Arpc3,Mad1l1,Tcf4,Fstl1,Fadd,Satb2,Kif1b,Gata2,Ctnnbip1,Dnm1l,Resp18,Prpf31,Hipk2,Vax1,Nes,Scx,Jup,Cttn,Afmid,Setd2
#
#the 16 genes are lethal in newly known gene set, but they are viable in training set.
#the16 lethal genes:
#Pcgf2,Fyn,Slc6a5,Socs1,Ift88,Fam20c,Eya4,Abca1,Pthlh,Lmna,Clcf1,Bag3,Casq2,Cbx2,Bbs7,Satb1,

import sys
import re

def main():
	infile = './known essentiality genes/training_all_genesinfo.csv'
	outfile = './known essentiality genes/training_all_genesinfo_corrected_classes.csv'

	viable = 'Syt1,Rest,Flii,Lama4,Atp2b1,Pcsk1,Il6st,Inpp5e,Sirt1,Pdpn,Traf6,Arpc3,Mad1l1,Tcf4,Fstl1,Fadd,Satb2,Kif1b,Gata2,Ctnnbip1,Dnm1l,Resp18,Prpf31,Hipk2,Vax1,Nes,Scx,Jup,Cttn,Afmid,Setd2'
	lethal = 'Pcgf2,Fyn,Slc6a5,Socs1,Ift88,Fam20c,Eya4,Abca1,Pthlh,Lmna,Clcf1,Bag3,Casq2,Cbx2,Bbs7,Satb1'
	viableL = viable.split(",")
	lethalL = lethal.split(",")
	
	infileL = [line.strip() for line in open(infile)]

	if len(infileL)==0:
		print(infile+" is empty.")
		sys.exit(-1)

	fw=open(outfile,"w")
	line = infileL[0]
	fw.write(line+"\n")
	for line in infileL[1:len(infileL)]:
		valsL = line.split(',')
		gene_name = valsL[0]
		if gene_name in viableL:
			valsL[len(valsL)-1]='Viable'
		elif gene_name in lethalL:
			valsL[len(valsL)-1]='Lethal'
		fw.write(valsL[0])
		for val in valsL[1:len(valsL)]:
			fw.write(','+val)
		fw.write("\n")
	fw.close()


if __name__ == "__main__":
	main()


			

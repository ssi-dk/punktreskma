#!/usr/bin/env python3
#06-06-18 updated to parse configfile and initiate different species parsing
#25-06-18 Added nucleotide recognition
#17-12-18 Began working to set it up for KMA instead
#Final output should add "<organism>_<gene>" columns to the regular resfinder output.

import argparse
import sys
import datetime
import os
import time
import subprocess
sys.path.append("/srv/data/tools/git.repositories/SSI-scripts/maos")

parser = argparse.ArgumentParser(
	description="Finds information on the point mutation positions.")
parser.add_argument("--frompunkt", help="Flag to indicate it's been called from punktres.sh", action="store_true")
parser.add_argument("-d", help="directory containing analysis files.")
parser.add_argument("-s", help="Split args.d flag", action="store_true")
parser.add_argument("--organism", help="Sets a specific organism")
parser.add_argument("--R1", help="Forward read if supplied")
parser.add_argument("--R2", help="Reverse read if supplied")
parser.add_argument("--id", help="Id to use for same-folder running")
parser.add_argument("-r", "--fromotherres", help="Run as part of regular res, if no species can be found, run for all species.", action="store_true")
parser.add_argument("-o", "--outputfile", help="Desired output file", default="punktres.tsv")
args = parser.parse_args()

if args.d:
	os.chdir(args.d)
dir=os.getcwd()
isolate=dir.split("/")[-1]
organisms=["Salmonella", "Ecoli", "Cjejuni"]#Replace to start as ALL OF THEM "Cdifficile", 
if args.organism:
	organisms=[args.organism]



codontable = dict()
codontable[11] = {
	'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
	'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
	'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
	'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
	'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
	'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
	'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
	'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
	'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
	'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
	'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
	'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
	'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
	}
		
def translate(dnaseq, table=11):
	proteinseq = list()
	for n in range(0, len(dnaseq)-2, 3):
		if "-" in dnaseq[n:n+3]:
			proteinseq.append("-")
		else:
			proteinseq.append(codontable[table][dnaseq[n:n+3]])
	return "".join(proteinseq)

def dividecopynumbers (copies, lineelements):
	depth=int(lineelements[1]) + int(lineelements[2]) + int(lineelements[3]) + int(lineelements[4])
	avgdepth=int(depth)/int(copies) #How many reads we assume pr. copy
	alloted=0						#currently assigned copies
	outcontent={}					#Will contain the output
	count=0
	order=["REF", "A", "C", "G", "T", "?", "-", "-"]
	mapcontent={}
	for num in range(1, 5):
		if int(lineelements[num]) > 0:
			mapcontent[order[num]]=lineelements[num]
	for nucleotide in mapcontent:	#Each observed nucleotide
		outcontent[nucleotide]=round(int(mapcontent[nucleotide])/avgdepth)	#mapcontent[nucleotide]=readdepth for current, the rounded rate is the assumed copies
		alloted+=round(int(mapcontent[nucleotide])/avgdepth)
	if len(mapcontent)==0:
		outcontent["-"]=int(copies)
		alloted=int(copies)
	while not alloted==int(copies):
		count+=1
		errors={}
		abserrors={}
			
		for nucleotide in mapcontent:
			errors[nucleotide]=float(mapcontent[nucleotide])-avgdepth*float(outcontent[nucleotide])		#The error here is how the depth of the nucleotide compares to it's number of copies multiplied by the average depth pr. copy
			abserrors[nucleotide]=abs(float(mapcontent[nucleotide])-avgdepth*float(outcontent[nucleotide]))
		try:
			editnucl=max(abserrors, key=lambda key: errors[key])		#Finds the nucleotide with the largest error
		except:
			print(copies)
			print(mapcontent)
			print(lineelements)
			print(abserrors)
			raise
		direction=int(errors[editnucl]*abs(errors[editnucl])**-1)	#Normalizes to 1 but keeps the sign
		alloted+=direction
		outcontent[editnucl]+=direction
		if count >=10:
			alloted=int(copies)
		
	outlist=[]
	for nucl in outcontent:
		outlist.extend([nucl]*outcontent[nucl])
	return(outlist)


def matrixparse(infile):
	os.system("".join(["gunzip ", infile]))
	order=["REF", "A", "C", "G", "T", "?", "-", "-"]
	inputfile=open(infile[:-3], 'r')
	outputseq=[]
	outthingorder={"R":0, "S":1, "u":2}
	outthing={}
	nuclout=[]
	presgenes=[]
	for line in inputfile:
		elements=line.rstrip().split("\t")
		elements.append("10")
		if len(elements) < 3:
			continue
		if elements[0]=="-":
			continue
		else:
			best=7
			for num in range(1,6):
				
				if int(elements[num]) > int(elements[best]):
					best=num
			outputseq.append(order[best])
			modnum=len(outputseq)
			if modnum in nuclposdict[organism]:
				
				thing=dividecopynumbers(int(nuclposdict[organism][modnum][-1]), elements)
				
				if not "_".join([nuclposdict[organism][modnum][2], str(len(outputseq)), "1"]) in headout[organism][2]:
					for n in range(1, len(thing)+1):
						headheadout[organism][2].append("_".join([organism, nuclposdict[organism][modnum][3]]))
						headout[organism][2].append("_".join([nuclposdict[organism][modnum][2], str(len(outputseq)-genes[organism][nuclposdict[organism][modnum][2]][0]), str(n)]))
				
				for nucl in thing:
					if nucl.lower() in nuclposdict[organism][modnum][0]:
						nuclout.append("/".join([nuclposdict[organism][modnum][0][nucl.lower()], nucl.lower()]))
						if nuclposdict[organism][modnum][0][nucl.lower()] == "R":
							presgenes.append(nuclposdict[organism][modnum][2])	#Indicate which genes have been shown as resistant
					else:
						nuclout.append("/".join(["U", nucl.lower()]))
				#if order[best].lower() in nuclposdict[organism][modnum][0]:
				#	print(nuclposdict[organism][modnum][0][order[best].lower()])
					#if not n in outthing:
					#	outthing[n]=[0,0,0]
					#outthing[n][outthingorder[nuclposdict[organism][modnum][0][order[best].lower()]]]+=1
					
				#else:
				#	print("unknown")
					#if not n in outthing:
					#	outthing[n]=[0,0,0]
					#outthing[n][outthingorder["u"]]+=1
						
			
	return(outputseq, nuclout, presgenes)
	
	
	
def output(organism, id, AA, nucl):
	pass
	
	
	
posdict={}
nuclposdict={}
genes={}

for line in open("/srv/data/DB/punktres/orginfo.config", 'r'):
	line=line.rstrip()
	lineelements=line.split("\t")
	if not lineelements[0] in genes:
		genes[lineelements[0]]={}
	for gene in lineelements[1:]:	#Each element is gene:start:end
		genethings=gene.split(":")
		genes[lineelements[0]][genethings[0]]=[int(genethings[1]), int(genethings[2]), {}, {}]

#sys.exit()
header=[]
headout={}
headheadout={}
genetoABclass={}
for organism in organisms:
	mutposinfo=open("".join(["/srv/data/DB/punktres/", organism, "_mutinfo.tsv"]))	
	posdict[organism]={}
	nuclposdict[organism]={}
	bygene={}
	orgheader=[[],[],[],[]]
	if not organism in headout:
		headheadout[organism]=[[],[],[]]
		headout[organism]=[[],[],[]]
	for line in mutposinfo:
		line=line.rstrip()
		if line.startswith("#"):
			continue
		elements=line.split("\t")
		#print(elements)
		if not elements[1] in orgheader[0]:
			orgheader[0].append(elements[1])
		if elements[3] in ["a", "c", "g", "t"]:
			if not "_".join([elements[1], elements[2]]) in orgheader[3]:
				orgheader[3].append("_".join([elements[1], elements[2]]))
				if len(elements) < 7:
					orgheader[2].append("_".join([elements[1], elements[2], "1"]))
				else:
					for n in range(1, int(elements[6])+1):
						orgheader[2].append("_".join([elements[1], elements[2], str(n)]))
				#orgheader[2].append("_".join([elements[1], elements[2]]))
		elif not "_".join([elements[1], elements[2]]) in orgheader[1]:
			orgheader[1].append("_".join([elements[1], elements[2]]))
		
		if elements[1] in genes[organism]:
			genetoABclass[elements[1]]=elements[5]
			actualpos=int(elements[2])+int(genes[organism][elements[1]][0])
			if elements[3] in ["a", "c", "g", "t"]:
				elements.append("1") #a bit of a cheat, appends a "copynumber" in case there isn't one currently.
				genes[organism][elements[1]][3][actualpos]="nucl"
				if actualpos in nuclposdict[organism]:
					nuclposdict[organism][actualpos][0][elements[3]]=elements[4]#.append([actualpos, elements[1], elements[3], elements[4], elements[5], elements[6]]) #pos to look for, gene, observed nucl, fenotype, ABclass, copynumer
				else:
					nuclposdict[organism][actualpos]=[{elements[3]:elements[4]},actualpos, elements[1], elements[5], elements[6]] #pos to look for, gene, observed nucl, fenotype, ABclass, copynumer
			elif int(elements[2]) in posdict[organism]:
				genes[organism][elements[1]][3][elements[2]]="AA"
				#genes[organism][elements[1]][2][elements[2]][elements[3]]=elements[4]
				posdict[organism][actualpos][0][elements[3]]=elements[4]
			else: #initiate the position in posdict
				#genes[organism][elements[1]][2][elements[2]]={elements[3]:elements[4]}
				posdict[organism][actualpos]=[{elements[3]:elements[4]},actualpos, elements[1], elements[5]]	#Note this is one element shorter, no copy numbers for amino		
				genes[organism][elements[1]][3][elements[2]]="AA"
			if elements[2] in genes[organism][elements[1]][2]:
				genes[organism][elements[1]][2][elements[2]][elements[3]]=elements[4]
			else:
				genes[organism][elements[1]][2][elements[2]]={elements[3]:elements[4]}
	

#The isolate is checked against all kept organisms. Later to be summed up for what is present.		
seqs={}
nucloutput={}
os.system("".join(["rm ", args.id, "*mat*"]))
headheadoutthing=["id"]
headoutthing=["id"]
outthing=[args.id]
for organism in organisms:
	os.system("".join(["/tools/git.repositories/kma/kma_index -i /srv/data/DB/punktres/", organism, "_punktgenes.fa -o /srv/data/DB/punktres/", organism, "_kma "]))
	print("".join(["/tools/git.repositories/kma/kma -t_db /srv/data/DB/punktres/", organism, "_kma -ipe ", args.R1, " ", args.R2, " -matrix -o ", args.id, "_", organism ]))
	os.system("".join(["/tools/git.repositories/kma/kma -t_db /srv/data/DB/punktres/", organism, "_kma -ipe ", args.R1, " ", args.R2, " -matrix -o ", args.id, "_", organism ]))
	(seqs[organism], nucloutput[organism], presgenes)=matrixparse("".join([args.id, "_", organism, ".mat.gz"]))	
	AAout=[]
	nuclout=[]
	genepresabsence=[]
	for gene in genes[organism]:
		if not gene in headout[organism][0]:
			headheadout[organism][0].append("_".join([organism, genetoABclass[gene]]))
			headout[organism][0].append(gene)
		if gene in presgenes:
			genepresabsence.append("1")
		else:
			genepresabsence.append("0")
		#if not genes[organism][gene][3]=="AA":
		#	continue
		AAseq=translate("".join(seqs[organism])[genes[organism][gene][0]:genes[organism][gene][1]])
		for pos in genes[organism][gene][2]:
			if pos not in genes[organism][gene][3] or  not genes[organism][gene][3][pos]=="AA":
				continue
			#MAOS note, the database is 1 indexed and the script is 0 indexed, thus a -1 modifier is needed
			if not "_".join([gene, pos]) in headout[organism][1]:
				headheadout[organism][1].append("_".join([organism, genetoABclass[gene]]))
				headout[organism][1].append("_".join([gene, pos]))
			if AAseq[int(pos)-1] in genes[organism][gene][2][pos]:
				AAout.append("/".join([genes[organism][gene][2][pos][AAseq[int(pos)-1]], AAseq[int(pos)-1]]))
				if genes[organism][gene][2][pos][AAseq[int(pos)-1]] == "R":
					genepresabsence[-1]="1"
			else:
				AAout.append("/".join(["U", AAseq[int(pos)-1]]))
			
	#outout=[]
	outthing.extend(genepresabsence)
	outthing.extend(AAout)
	outthing.extend(nucloutput[organism])
	
	output(organism, args.id, AAout, nuclout)		

	for line in headheadout[organism]:
		headheadoutthing.extend(line)
	for line in headout[organism]:
		headoutthing.extend(line)
		
"""		
if not os.path.isfile(args.outputfile):
	os.system("".join(["echo ", "\t".join(headheadoutthing), "\n", "\t".join(headoutthing), "\n", "\t".join(outthing), " > ", args.outputfile]))
os.system("".join(["echo ", "\t".join(headheadoutthing), "\n", "\t".join(headoutthing), "\n", "\t".join(outthing)]))

proc=subprocess.Popen("".join(["grep ", args.id," ", args.outputfile]), stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
while len(out.decode()) < 1:
	time.sleep(2)
	os.system("".join(["echo ", "\t".join(outthing), " >> ", args.outputfile]))
	proc=subprocess.Popen("".join(["grep '",args.id,"' ", args.outputfile]), stdout=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
"""

outputfile=open("_".join([args.id, "punktout.tsv"]),'w')
print("\t".join(headheadoutthing), file=outputfile)
print("\t".join(headoutthing),file=outputfile)
print("\t".join(outthing),file=outputfile)
#print(seqs)
#sys.exit()
#Time to filter and map the reads
	
#Load the database of positions.
punktconfig=open("/srv/data/DB/punktres/orginfo.config", 'r') #File containing the reading info for the different organisms
config={"filter":"", "mutpath":"", "genelocations":""}

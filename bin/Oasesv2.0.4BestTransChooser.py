#!/usr/bin/env python

# Ex : Oasesv2.0.4BestTransChooser.py -l 0.8

import sys
from Bio import SeqIO
from numpy import mean
from numpy import median
from scipy.stats import *
from optparse import *

usage  = "%prog -l <transcript length cutoff fraction>"
desc   = """This script will read the output files from Oases and generate a file containing only the best transcript for each Locus.

You need to navigate to the folder that contains the output from your Oases assembly and you need to have not renamed or modified the following files from Oases:
transcripts.fa, contig-ordering.txt, and stats.txt.

See: https://sites.google.com/a/brown.edu/bioinformatics-in-biomed/velvet-and-oases-transcriptome
"""
parser = OptionParser(usage = usage, description = desc)
igroup = OptionGroup(parser, "Parameter","")
igroup.add_option("-l", "--length", dest="trans_l", type="float", help="The transcript length cutoff fraction. Default 0.8")
parser.add_option_group(igroup)
(options, args) = parser.parse_args()

try:
	TransProblem=open('transcripts.fa','r')
	TransProblemString=''
except IOError:
	TransProblem=True
	TransProblemString='transcripts.fa'
try:
	ContigProblem=open('contig-ordering.txt','r')
	ContigProblemString=''
except IOError:
	ContigProblem=True
	ContigProblemString='contig-ordering.txt'
try:
	StatsProblem=open('stats.txt','r')
	StatsProblemString=''
except IOError:
	StatsProblem=True
	StatsProblemString='stats.txt'
if TransProblem==True or ContigProblem==True or StatsProblem==True:	# or SpliceProblem==True:
	print 'This program requires the unadultered raw output from an Oases assembly. You appear to have renamed or moved or deleted some of these file(s).'
	print 'The following file(s) have issues:',TransProblemString,ContigProblemString,StatsProblemString	#,SpliceProblemString
	print 'Once you have fixed the file(s), please come back. Goodbye!'
	sys.exit(0)
	
SMALLEST_LENGTH_FRACTION = options.trans_l or float(0.8)

if SMALLEST_LENGTH_FRACTION<=0 or SMALLEST_LENGTH_FRACTION>1.0:
	print 'Please enter a transcript length cutoff fraction between 0.01 and 1.00 (converted to 100% scale)'

#print 'Thank you for your input, your opinion is valuable to us. Please stand by...'


ErrorCount=0
BabyErrorCount=0
#ErrorOut=open('BestTransChooser_ErrorLog.txt',"w")
FaRead = open('transcripts.fa',"r")
FaLine=FaRead.readline().strip()	#enters into first line of file
Long_Trans={}#key is locus, value is the length of longest transcript
All_Trans={}#key is locus+Transcript, value is the length of transcript
Trans_Per_Locus={}#key is locus, value is list of all transcript numbers
#print 'Reading transcripts.fa. Measuring length of all transcripts and the longest transcript in each locus.'
PrevLocus=''
LongestTrans=0
AllTransLengthsList = []
AllTransLengthTotal = 0
while FaLine != '':
	if FaLine[0:7] == '>Locus_':
		FaFields=FaLine.split('_')
		Locus=int(FaFields[1])
		TransFields=FaFields[3].split('/')
		Transcript=int(TransFields[0])
		LocusTransID=str(str(Locus)+'T'+str(Transcript))
		if Locus not in Trans_Per_Locus:
			TransLocusList=list()
		else:
			TransLocusList=list(Trans_Per_Locus[Locus])
		TransLocusList.append(Transcript)
		Trans_Per_Locus[Locus]=TransLocusList
		FaLine=FaRead.readline()
		while FaLine=='\n':
			FaLine=FaRead.readline()
		FaLine=FaLine.strip()
	FaSeq = []
	while FaLine[0:7] != '>Locus_' and FaLine != '':	#loop through fasta file to generate fasta sequence
		FaSeq.append(FaLine)
		FaLine=FaRead.readline()
		while FaLine=='\n':
			FaLine=FaRead.readline()
		FaLine=FaLine.strip()
	FaSeq=''.join(FaSeq)	#convert sequence to single line
	AllTransLengthTotal += len(FaSeq)
	AllTransLengthsList.append(len(FaSeq))
	All_Trans[LocusTransID] = len(FaSeq)
	if PrevLocus==Locus and len(FaSeq) > LongestTrans:
		LongestTrans=len(FaSeq)
		Long_Trans[Locus]=len(FaSeq)
	elif PrevLocus!=Locus and Locus not in Long_Trans:
		LongestTrans=len(FaSeq)
		Long_Trans[Locus]=len(FaSeq)
	PrevLocus=Locus
FaRead.close()
AllTransLengthsList.sort(reverse=True)
AllTransN50Target = AllTransLengthTotal * float(.50)
RunningTotal = 0
for i in AllTransLengthsList:
	RunningTotal += i
	if RunningTotal > AllTransN50Target:
		AllTransN50=str(i)
		break

#print "Reading stats.txt. Measuring the fold coverage of each node (parts of each transcript)."
stats = open('stats.txt',"r")
first = True
node_dict = {} #key is nodeid, value is exact coverage
for i in stats:
	if first == True:
		first = False
		continue
	spls = i.strip().split("\t")
	if spls[6] != 'Inf':
		node_dict[int(spls[0])] = float(spls[6])
stats.close()

#print "Reading contig-ordering.txt. Measuring coverage of all transcripts and the transcript with the highest coverage in each locus."
ContigOrdRead = open('contig-ordering.txt',"r")
High_Fold={}#key is locus, value is the highest fold coverage of locus
All_Fold={}#key is locus+Transcript, value is the geometric mean of fold coverage of transcript
PrevCOrdLocus=''
for record in SeqIO.parse(ContigOrdRead, "fasta"):
	if "Transcript" in record.id:
		COrdFields=record.id.split('_')
		COrdLocus=int(COrdFields[1])
		TransFields=COrdFields[3].split('/')
		COrdTranscript=int(TransFields[0])
		COrdLocusTransID=str(str(COrdLocus)+'T'+str(COrdTranscript))
		NodeFields=str(record.seq).split("->")
		Count=0
		nums = []
		GeoMult=float(1.0)
		for j in NodeFields:
			tid = j.split(":")[0]
			if tid[0] == "-":
				tid = tid[1:]
			if int(tid) in node_dict:
				Count+=1
				nums.append(float(node_dict[int(tid)]))
		if Count>0:
			GeoMean=gmean(nums)
		else:
			GeoMean=float(2.0)	#this represents transcripts in which all nodes have 'inf' coverage
		All_Fold[COrdLocusTransID] = float(GeoMean)
		if All_Trans[COrdLocusTransID] >= Long_Trans[COrdLocus]*SMALLEST_LENGTH_FRACTION:	#only transcripts that are larger than SMALLEST_LENGTH_FRACTION of longest can be considered for highest fold coverage
			if PrevCOrdLocus==COrdLocus and GeoMean > HighestFold:
				HighestFold=GeoMean
				High_Fold[COrdLocus]=GeoMean
			elif PrevCOrdLocus!=COrdLocus and COrdLocus not in High_Fold:
				HighestFold=GeoMean
				High_Fold[COrdLocus]=GeoMean
			PrevCOrdLocus=COrdLocus
ContigOrdRead.close()

SmallFractionString='%.f' % float(SMALLEST_LENGTH_FRACTION*100)
ChosenNucSeqs = open('transcripts.best.'+SmallFractionString+'%Longest.fa',"w")
FaRead = open('transcripts.fa',"r")
Chosen_Len_Precent={}#key is locus, value is the percent length of the chosen transcript in the locus 
FaLine=FaRead.readline().strip()	#enters into first line of file
FractionString='%.2f' % float(SMALLEST_LENGTH_FRACTION*100)
#print 'Reading transcripts.fa for output.'
print 'Copying the transcripts with the highest fold coverage that are longer than '+FractionString+'% of the longest transcript in each locus.'
First=True
PrevLocus=''
ChosenTransLengthsList = []
ChosenTransLengthTotal = 0
ChosenTransCount=0
while FaLine != '':
	if FaLine[0:7] == '>Locus_':
		FaFields=FaLine.split('_')
		Locus=int(FaFields[1])
		TransFields=FaFields[3].split('/')
		Transcript=int(TransFields[0])
		LocusTransID=str(str(Locus)+'T'+str(Transcript))
		FullTitle=str(FaLine)
		FaLine=FaRead.readline()
		while FaLine=='\n':
			FaLine=FaRead.readline()
		FaLine=FaLine.strip()
	FaSeq = []
	if PrevLocus!=Locus:
		if First==True:
			First=False
			LocusDone=False
		else:
			LocusDone=False
			###BEGINNING OF DATABASE OUTPUT - ALL BUT LAST ENTRY - COPIED BELOW FOR LAST ENTRY
			ChosenNucSeqs.write(OutString+'\n')
			BigFields=OutString.split('\n')
			OutFields=BigFields[0].split('_')
			TempLocus=int(OutFields[1])
			TempTransFields=OutFields[3].split('/')
			Confidence=str(OutFields[5])
			TempTranscript=int(TempTransFields[0])
			TotalTrans=str(TempTransFields[1])
			TempLocusTransID=str(str(TempLocus)+'T'+str(TempTranscript))
			RNAseq=str(BigFields[1])
			RNALen=str(len(RNAseq))
			if int(RNALen)==Long_Trans[TempLocus]:
				LongestTrans='True'
			else:
				LongestTrans='False'
			OtherTransList=list(Trans_Per_Locus[TempLocus])
			OtherTransList.remove(TempTranscript)
			TempTransList=list(OtherTransList)
			for i in OtherTransList:
				if int(i) <  Long_Trans[TempLocus]*SMALLEST_LENGTH_FRACTION:	#we don't want ALL the transcripts, only those that are longer than the minimum cutoff
					TempTransList.remove(i)
			OtherTransList=list(TempTransList)
			if OtherTransList==[]:
				OtherTrans='False'
			else:
				OtherTrans=str(OtherTransList).strip('[]').replace(',','')
			Chosen_Len_Precent[TempLocus]='%.2f' % float((All_Trans[TempLocusTransID]*100)/Long_Trans[TempLocus])
			###END OF DATABASE OUTPUT - ALL BUT LAST ENTRY - COPIED BELOW FOR LAST ENTRY
			
	if All_Fold[LocusTransID]==High_Fold[Locus] and All_Trans[LocusTransID]>=Long_Trans[Locus]*SMALLEST_LENGTH_FRACTION and LocusDone==False:
		while FaLine[0:7] != '>Locus_' and FaLine != '':	#loop through fasta file to generate fasta sequence
			FaSeq.append(FaLine)
			FaLine=FaRead.readline()
			while FaLine=='\n':
				FaLine=FaRead.readline()
			FaLine=FaLine.strip()
		FaSeq=''.join(FaSeq)	#convert sequence to single line
		ChosenTransLengthTotal += len(FaSeq)
		ChosenTransLengthsList.append(len(FaSeq))
		OutString=FullTitle+'\n'+FaSeq
		LocusDone=True
		ChosenTransCount+=1
	else:
		while FaLine[0:7] != '>Locus_' and FaLine != '':	#loop through fasta file to get to next entry
			FaLine=FaRead.readline()
			while FaLine=='\n':
				FaLine=FaRead.readline()
			FaLine=FaLine.strip()
	PrevLocus=Locus

###BEGINNING OF DATABASE OUTPUT - LAST ENTRY ONLY
ChosenNucSeqs.write(OutString+'\n')
BigFields=OutString.split('\n')
OutFields=BigFields[0].split('_')
TempLocus=int(OutFields[1])
TempTransFields=OutFields[3].split('/')
Confidence=str(OutFields[5])
TempTranscript=int(TempTransFields[0])
TotalTrans=str(TempTransFields[1])
TempLocusTransID=str(str(TempLocus)+'T'+str(TempTranscript))
RNAseq=str(BigFields[1])
RNALen=str(len(RNAseq))
if int(RNALen)==Long_Trans[TempLocus]:
	LongestTrans='True'
else:
	LongestTrans='False'
OtherTransList=list(Trans_Per_Locus[TempLocus])
OtherTransList.remove(TempTranscript)
TempTransList=list(OtherTransList)
for i in OtherTransList:
	if int(i) <  Long_Trans[TempLocus]*SMALLEST_LENGTH_FRACTION:	#we don't want ALL the transcripts, only those that are longer than the minimum cutoff
		TempTransList.remove(i)
OtherTransList=list(TempTransList)
if OtherTransList==[]:
	OtherTrans='False'
else:
	OtherTrans=str(OtherTransList).strip('[]').replace(',','')
Chosen_Len_Precent[TempLocus]='%.2f' % float((All_Trans[TempLocusTransID]*100)/Long_Trans[TempLocus])
###END OF DATABASE OUTPUT - LAST ENTRY ONLY

FaRead.close()
ChosenNucSeqs.close()
ChosenTransLengthsList.sort(reverse=True)
ChosenTransN50Target = ChosenTransLengthTotal * float(.50)
RunningTotal = 0
for i in ChosenTransLengthsList:
	RunningTotal += i
	if RunningTotal > ChosenTransN50Target:
		ChosenTransN50=str(i)
		break


print 'You have a total of '+str(ChosenTransCount)+' transcripts in your file.'
print 'The N50 of the transcriptome in the file is: '+ChosenTransN50
#if ErrorCount==0:
#	print 'Congratulations! There were no serious errors that occured during the script, enjoy the data.'
#	ErrorString='No serious errors were detected during the processing of this script.'
#	ErrorOut.write(ErrorString+'\n')
#else:
#	print 'WARNING! There were '+str(ErrorCount)+' serious errors that occured in the processing of this script. Please read the BestTransChooser_ErrorLog.txt file.'
#print 'There were '+str(BabyErrorCount)+' minor errors that also occured but these can be safely ignored. For details read the BestTransChooser_ErrorLog.txt file.'
#ErrorOut.close()

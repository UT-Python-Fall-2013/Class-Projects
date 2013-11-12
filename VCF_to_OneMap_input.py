#!/usr/bin/env python
##VCF_to_OneMap_input.py
##written 11/6/13 by Groves Dixon
ProgramName = 'VCF_to_OneMap_input.py'
LastUpdated = '11/12/13'
By = 'Groves Dixon'
VersionNumber = '1.0'

##Assign Global User Help Variables

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
This script is intended to conver VCF files into input files for OneMap
Genotype data for parents from an input VCF file is used to assign 
Cross types. Then Genotype data for each sample in the VCF file is output
with the appropriate cross type and format for input into OneMap.
'''

AdditionalProgramInfo = '''
Additional Program Information:\n
Add info here
'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
Start_time = time.time() ##keeps track of how long the script takes to run

##Set Up Argument Parsing

parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-input', '-i', required = True, metavar = 'inputName.vcf', dest = 'In', help = 'The name of the input vcf file')
parser.add_argument('-output', '-o', default='vcf2Onemap.out', metavar = 'outputName.out', dest = 'Out', help = 'The desired name for the output file')
parser.add_argument('-parent1', '-p1', required = True, metavar = 'first parent', dest = 'P1', help = 'The Sample Name for parent 1')
parser.add_argument('-parent2', '-p2', required = True, metavar = 'second parent', dest = 'P2', help = 'The Sample Name for parent 2')
parser.add_argument('-version','-v', action='version', version=VersionString)
#parser.print_help()
args = parser.parse_args()


#Assign Arguments
InfileName = args.In
OutfileName = args.Out
Parent1 = args.P1
Parent2 = args.P2
ParentList = [Parent1, Parent2]

#Assign some global Variables
RequiredHeaderParts = ['#CHROM', 'POS', 'REF', 'ALT', 'INFO', 'FORMAT']


#Functions

def read_file(infileName, delim):
	'''Function to read in a file as a list of lists
	'''
	lineNumber = 0
	fileList = []
	with open(infileName, 'r') as infile:
		for line in infile:
			lineNumber += 1
			line = line.strip('\n').split(delim)
			##skip lines starting with ##
			if line[0][0:2] == '##':
				continue
			#record the header
			if line[0][0] == '#':
				print(line)
				header = line
			#record all other lines as data
			else:
				fileList.append(line)
	return header, fileList

def parse_heading(header, required):
	'''Function to check the header and gather metadata.
	Looks through the header and make sure
	all the necessary header elements are present. Returns the 
	indices for the required header elements.
	'''
	headerIndices = []
	for i in required:
		if i not in header:
			exit('\nERROR:\nThe required component {} was not found in header. Please check input file.\n'.format(i))
		else:
			headerIndices.append(header.index(i))
	print('\nChecking that header has necessary components\n')
	print('\nAll necessary header components have been found\n')
	print('\nHeader indices are as follows:')
	for x in range(len(required)):
		print('{}: {}'.format(required[x], str(headerIndices[x])))
	chromIndex = header.index('#CHROM')
	infoIndex = header.index('INFO')
	posIndex = header.index('POS')
	refIndex = header.index('REF')
	altIndex = header.index('ALT')
	formatIndex = header.index('FORMAT')
	return chromIndex, infoIndex, posIndex, refIndex, altIndex, formatIndex

def get_sample_list(header, parentList):
	'''The samples are included in the header in VCF files.
	This function pulls a list of just the samples from the header.
	Also checks the sample list to make sure parents are present.
	'''
	start = header.index('FORMAT')+1 ##in VCF samples always start after FORMAT
	sampleList = header[start:]
	print('\nNumber of Samples Found in VCF: {}\n'.format(len(sampleList)))	
	print('\nChecking the parents can be found among the samples\n')
	##check for parents and exit if not there
	for i in parentList:
		if i not in sampleList:
			exit('\nERROR:\nCould not find parent labeled {} among the samples\
			please check VCF file and parent entries'.format(i))
	print('\nBoth Parent Samples Found\n')
	return sampleList
		
def get_variant_list(filelist, header):
	'''Function pulls two lists from the vcf filelist:
	a list of all the variants recorded as CHROM_POS
	and a coindexed list with the INFO column for each variant
	Indexing in this list will be used to keep track of which variants
	the sample genotypes are connected to
	'''
	#create two empty lists to store info in
	variantList = []
	infoList = []
	#being loop through data
	for line in filelist:
		chrom = line[ChromIndex] ; pos = line[PosIndex] ; info = line[InfoIndex] ##using indices assinged by parse_header()
		#reformat the variant as 'CHROM_POS'
		variant = '{}_{}'.format(chrom,pos)
		##add the variant and info for this line to the two lists
		infoList.append(info)
		variantList.append(variant)
	numberVariants = len(variantList)
	print('\nReformatting genotypes for {} variant sites\n'.format(numberVariants))
	return variantList, infoList
	
def gather_genotypes(fileList, sampleList, variantList, header):
	'''This function iterates through the file list and records
	the genotype data for each sample in a dictionary. The genotype data for each 
	variant position is recorded as a list for the sample. The list for
	each sample coordinates with the VariantList. SNP calls 
	are stored in dictionary as A,T,G,and C. The base calls are retrieved using
	metadata gathered by parse_heading().
	
	Dictionary return format:
	{sample1Name: [['A','T'],['G','G'],['-','-']...], sample2Name: [[etc.]]}
	'''
	#empty dictionary to store genotype data
	genotypeDict = {}
	#make a key for each sample with list to store genotype data
	for sample in sampleList:
		genotypeDict[sample] = [] ##list to store genotypes for each variant
	lineNumber = 0
	##begin loop through each line of data
	for line in fileList:
		lineNumber += 1
		##get ref and alt based on indices output by parse_heading()
		ref = line[RefIndex] #single base
		alt = line[AltIndex].split(',')##can have multiple alternates delimited by ','
		##get the index for the genotype data from the metadata provided in FORMAT.
		genotypeIndex = line[FormatIndex].split(':').index('GT')##FORMAT is colon delimited. GT marks index for genotype.
		samplesStart = FormatIndex+1 ##samples always begin after FORMAT column in VCF
		##get slice of the line with the genotype data
		genotypeData = line[samplesStart:]
		#set up counter to keep track of index
		sampleIndex = -1
		##iterate through the genotype data in the line
		for i in genotypeData:
			i = i.split(':') ##the genotype entries (like the FORMAT elemnt) is colon delimited
			#keep track of index
			sampleIndex += 1
			#get sample name from sampleList using index
			sampleName = sampleList[sampleIndex]
			##use the genotypeIndex to grab the genotype data
			genotype = i[genotypeIndex]
			
			if '|' in genotype:
				genotype = genotype.split('|') ##genotypes can be split with '|' or '/' for phased and unphased data repsectively
			elif '/' in genotype:
				genotype = genotype.split('/')
			else:
				print('/nunrecognized genotype delimiter\
				for sample {} for variant number {}'.format(sampleName, variantList[lineNumber-1]))
				continue
			#genotype is now recorded as a list of integers
			#convert the integers into A T G C or -
			convertedGenotype = [] #variable to store the genotype in converted form
			#loop through the two alleles for this genotype
			for allele in genotype:
				#convert for missing data
				if allele == '.': ## . marks missing data in VCF file
					convertedGenotype.append('-') ## - marks missing data in OneMap
				#convert for present data
				else:
					allele = int(allele)##alleles are stored as numbers in VCF. Make integer for indexing
					#convert if it is reference allele
					if allele == 0:
						convertedGenotype.append(ref)## zero indicates the allele is the reference
					#convert if it is an alternate allele 
					else:
						convertedGenotype.append(alt[(allele-1)]) # nonzero gives tells which of the alternate alleles it is (indexed from 1 in VCF so subtract for python)
			#done converting the alleles for the genotype
			#add the converted genotype to the correct sample's list in the dictionary
			genotypeDict[sampleName].append(convertedGenotype)
	return genotypeDict
	
def test_heterozygosity(genotype):
	'''Function to test if a genotype is heterozygous
	Returns either True or False
	Used to check parents for assigning cross types in 
	convert_genotypes()
	'''
	if genotype[0] == genotype[1]:
		answer = 'False'
	else:
		answer = 'True'
	return answer

def convert_genotypes(genotypeDict, parent1, parent2, variantList, sampleList):
	'''This funciton converts the genotypes recorded in genotypeDict
	into the proper format for OneMap. It returns converted genotype
	data in three forms: a dictionary of converted genotypes for all samples,
	a list of the variants that were converted, and a list of the
	cross type for each converted variant. 
	In each returned item, variants have the same index.
	'''
	##create dictionary to store converted genotypes in
	cGenotypeDict = {}
	noHetCount = 0 ##counter to keep track of markers with no heterozygote parent (uninformative for linkage and not converted)
	##add key to dictionary for each sample with empty list to store converted genotypes
	for i in sampleList:
		cGenotypeDict[i] = []
	#set up list to record converted Variants in
	cVariantList = []
	#set up list to record crossTypes in
	crossTypeList = []
	#create list with the four possible OneMap alleles
	oneMapAlleles = ['a','b','c','d']##will be used for converting allele names
	weirdOffspringList = [] ##list to store any offspring in case some have alleles not found in parents
	##set up the two parents' genotype lists for easy access
	p1List = genotypeDict[parent1]
	p2List = genotypeDict[parent2]
	##begin loop from 0 to total number of variants
	for i in range(len(p1List)):
		
		#list to store alleles coindexed with oneMapAlleles
		#this will be used during conversion
		alleleList = ['','','','']
		
		##grab the parental genotype for this marker
		p1Genotype = p1List[i]
		p2Genotype = p2List[i]
		
		#grab the individual alleles
		p1a1 = p1Genotype[0] ; p1a2 = p1Genotype[1]
		p2a1 = p2Genotype[0] ; p2a2 = p2Genotype[1]
		
		#build list of alleles
		parentAlleleList = [p1a1, p1a2, p2a1, p2a2]
		
		##designate whether p1 and p2 are heterozygotes
		p1Het = test_heterozygosity(p1Genotype)
		p2Het = test_heterozygosity(p2Genotype)
		
		##if both parents are homozygotes, the marker is uninformative and is skipped
		if p1Het == 'False' and p2Het == 'False':
			noHetCount += 1
			continue
		
		##deal with both parents being heterozygotes
		##This can return cross types A.1, A.2, or B3.7
		##which are (ab x cd), (ab x ac), and (ab x ab)
		
		elif p1Het == 'True' and p2Het == 'True':
				
			#if both alleles are shared then cross type is B3.7 (ab x ab)
			#record alleles in the alleleList
			if p1a1 in p2Genotype and p1a2 in p2Genotype:
				crossType = 'B3.7'
				alleleList[0] = p1a1 ; alleleList[1] = p1a2
				
			##if one allele is shared then cross type is A.2 (ab x ac)
			##the shared allele is called as 'a'
			##record the alleles in alleleList appropriately
			elif p1a1 in p2Genotype and p1a2 not in p2Genotype: ##if statement detects this case: ab x ac
				crossType = 'A.2'
				##now we have to find which allele for parent2 is not in parent 1
				##assign that allele as 'c'
				if p2a1 in p1Genotype:
					c = p2a2
				else:
					c = p2a1
				alleleList[0] = p1a1 ; alleleList[1] = p1a2 ; alleleList[2] = c
			elif p1a1 not in p2Genotype and p1a2 in p2Genotype: ##if statement detects this case: ba x ac
				crossType = 'A.2'
				##now we have to find which allele for parent2 is not in parent 1
				##assign that allele as 'c'
				if p2a1 in p1Genotype:
					c = p2a2
				else:
					c = p2a1
				alleleList[0] = p1a2 ; alleleList[1] = p1a1 ; alleleList[2] = c
				
			##if none of the alleles are shared the cross type is A1 (ab x cd)
			elif p1a1 not in p2Genotype and p1a2 not in p2Genotype:
				crossType = 'A1'
				##assignment is simply left to right for the 4 alleles
				alleleList[0] = p1a1 ; alleleList[1] = p1a2 ; alleleList[2] = p2a1 ; alleleList[3] = p2a2
			else:
				exit('Did not find a cross type for double Heterozygote!!')
	
		##Deal with cases when parent1 as heterozygote and parent2 as homozygote
		##This can return cross type D1.9 or D1.10
		## which are (ab x cc) and (ab x aa)
		elif p1Het == 'True' and p2Het == 'False':
			#if the homozygote allele is shared then it is labeled as 'a'
			if p2a1 in p1Genotype: ##detects this case: ab x aa
				crossType = 'D1.10'
				#assign the unshared allele as b
				if p1a1 not in p2Genotype:
					b = p1a1
				else:
					b = p1a2
				alleleList[0] = p2a1 ; alleleList[1] = b
			else: ##it must be ab x cc
				crossType = 'D1.9'
				#the unshared homozygous allele in parent 2 labeled as 'c'
				alleleList[0] = p1a1 ; alleleList[1] = p1a2 ; alleleList[2] = p2a1
				
		##Deal with cases when parent1 is homozygote and parent2 is heterozygote
		##This can return cross type D2.14 or D2.15
		## which are (c x cc) and (ab x aa)
		##This section is the reciprocal case of the precious section
		elif p1Het == 'False' and p2Het == 'True':
			#if the homozygote allele is shared then it is labeled as 'a'
			if p1a1 in p2Genotype:
				crossType = 'D2.15'
				#assign the unshared allele as b
				if p2a1 not in p2Genotype:
					b = p2a1
				else:
					b = p2a2
				alleleList[0] = p1a1 ; alleleList[1] = b
			else:
				crossType = 'D2.14'
				alleleList[0] = p2a1 ; alleleList[1] = p2a2 ; alleleList[2] = p2a1
		else:
			exit('Did not sort the marker based on the homo and heterozygosity!')
		
		##Finished solving the cross types based on parental genotypes
		##At this point the 'alleleList' variable is populated with 
		##alleles from the parents so their indices 
		##correctly match their conversions for OneMap (stored in the oneMapAlleles variable)
		
		#use this list pair to convert the alleles for this variant for each offspring
		##begin loop through all the samples
		for sample in sampleList: 
			##grab the genotype for the sample
			sampleGenotype = genotypeDict[sample][i]
			##create variable to store the converted genotype under
			convertedGenotype = []
			##begin loop through the 2 alleles for sample's genotype
			##to get the converted allele's
			for allele in sampleGenotype:
				##the sample could possess an allele absent from parents
				##this indicates a mutation, or more likely sequencing error
				##the sample will be recorded as having no data
				if allele not in alleleList:
					weirdOffspringList.append([sample, variantList[i]])
					convertedGenotype = ['-','-']
				else:
					#get the index for the allele
					index = alleleList.index(allele)
					#use the oneMapAlleles list to pull the correct allel
					convertedGenotype.append(oneMapAlleles[index])
					#sort genotype to be in alphabetical order
			#reformat the genotype from a list to a string in OneMap format
			convertedGenotype.sort()
			if convertedGenotype[0] == convertedGenotype[1]:
				convertedGenotype = convertedGenotype[0]
			else:
				convertedGenotype = ''.join(convertedGenotype)
			##record the genotype in OneMap format for the sample
			cGenotypeDict[sample].append(convertedGenotype)
			##record the variant in cVariantList
		cVariantList.append(variantList[i])
		crossTypeList.append(crossType)
	return cGenotypeDict, cVariantList, crossTypeList
					
		
def output(outfileName, cGenotypeDict, cVariantList, crossTypeList, sampleList):
	with open(outfileName, 'w') as out:
		##build the header line
		##which is simply number of samples and number of variants
		numberIndividuals = (len(sampleList) - 2)
		numberVariants = len(cVariantList)
		header = '{} {}'.format(numberIndividuals, numberVariants)
		out.write(header)
		##begin loop through all variants
		for i in range(len(cVariantList)):
			genotypeList = []
			#To build list of genotypes
			#loop through all samples
			for sample in sampleList:
				#Do not include the parents
				if sample == Parent1:
					continue
				if sample == Parent2:
					continue
				else:
					genotypeList.append(cGenotypeDict[sample][i])
			genotypeString = ','.join(genotypeList)
			outString = '\n*{} {}\t{}'.format(cVariantList[i], crossTypeList[i], genotypeString)
			out.write(outString)

##Initiate each function in turn to build OneMap file		
Header, FileList = read_file(InfileName, '\t')
#print(FileList[0])
SampleList = get_sample_list(Header, ParentList)
ChromIndex, InfoIndex, PosIndex, RefIndex, AltIndex, FormatIndex = parse_heading(Header, RequiredHeaderParts)
VariantList, InfoList = get_variant_list(FileList, Header)
GenotypeDict = gather_genotypes(FileList, SampleList, VariantList, Header)
OneMapGenoDict = convert_genotypes(GenotypeDict, Parent1, Parent2, VariantList, SampleList)
CGenotypeDict, CVariantList, CrossTypeList = convert_genotypes(GenotypeDict, Parent1, Parent2, VariantList, SampleList)
output(OutfileName, CGenotypeDict, CVariantList, CrossTypeList, SampleList)

#return time to run
Time = time.time() - Start_time
print('Time took to run: {}'.format(Time))

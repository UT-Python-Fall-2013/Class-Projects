#!/bin/env python3

import sys, os, argparse, logging, re, string
from utils import *

def main ():
	usage=("usage: %s <output of run_einverted_parallel> <genome file> <output file>" % sys.argv[0])
	parser = argparse.ArgumentParser(description='Folds einvertd', usage=usage)

	# Add optional switches
	parser.add_argument('-v', '--verbose', action='store_true', help='produce verbose output [%(default)s]')
	parser.add_argument('-l', '--log', action='store', dest='log', default=None, help='dump messages in log file [%(default)s]')
	parser.add_argument('-f', '--flanking_nt', action='store', type=int, default=10, help='number of flanking_nt added to each sides of inverted repeat [%(default)s]')
	opts,args = parser.parse_known_args()
	
	logLevel = logging.WARN
	if opts.verbose:
		logLevel = logging.INFO
	Logger(logLevel, opts.log)
	logger = logging.getLogger()

	flankingNT = opts.flanking_nt
	invertedRepeatsFileName = args[0]
	genomeFileName = args[1]
	outputFileName = args[2]
	if not os.path.isfile(invertedRepeatsFileName):
		logger.error("input file %s not found" % invertedRepeatsFileName)
		parser.print_help()
		sys.exit(1)
	if not os.path.isfile(genomeFileName):
		logger.error("input file %s not found" % genomeFileName)
		parser.print_help()
		sys.exit(1)
	try:
		tmpOutput = open(outputFileName+"_tmp", "w")
	except Exception as e:
		logger.error("unable to open temporary file %s" % outputFileName+"_tmp")
		logger.exception(e)
		parser.print_help()
		sys.exit(1)

	sequences = FastaParser(genomeFileName).hashMap() # create name to genome map
	input = open(invertedRepeatsFileName) # open input 
	dataRegex = re.compile(r'(\S+) (\d+) \d+ (\d+)') # this is the format
	lineCount = 0
	for line in input:
		lineCount += 1
		if lineCount % 1000000 == 0: # print progress message
			logger.info("processed %dK lines" % (lineCount/1000))
		line = line.strip()
		searchRes = dataRegex.search(line)
		if searchRes: # data matched
			chr = searchRes.group(1)
			x = int(searchRes.group(2)) - flankingNT
			y = int(searchRes.group(3)) + flankingNT

			seq = getSubseq(sequences[chr], x, y)# gets + strand
			seq = seq.upper()
			seq = seq.translate(bytes.maketrans(b'MRWSYKVHDBX', b'NNNNNNNNNNN'))
			print(">%s_%d_%d\n%s" % (chr, x, y, seq), file=tmpOutput) 
			
			seq = getSubseq(sequences[chr], y, x)# gets + strand
			seq = seq.upper()
			seq = seq.translate(bytes.maketrans(b'MRWSYKVHDBX', b'NNNNNNNNNNN'))
			print( ">%s_%d_%d\n%s" % (chr, y, x, seq), file=tmpOutput) 

	input.close()
	tmpOutput.close()
	logger.info("launching RNAfold with input %s and output %s" % (outputFileName+"_tmp", outputFileName))
	os.system("RNAfold -noPS < %s > %s" % (outputFileName+"_tmp", outputFileName));

if __name__ == "__main__":
	main()

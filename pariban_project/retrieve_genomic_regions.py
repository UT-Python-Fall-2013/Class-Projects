#!/bin/env python3
import sys, argparse, logging, re, utils, string

def main():
	usage=("usage: %s <window size> <miRNA match file> <input genome file> <output file>" % sys.argv[0])
	parser = argparse.ArgumentParser(description='Retrieves genomic regions', usage=usage)
	# Add optional switches
	parser.add_argument('-v', '--verbose', action='store_true', default=False, help='produce verbose output')
	parser.add_argument('-l', '--log', action='store', dest='log', default=None, help='dump messages in log file')
	opts,args = parser.parse_known_args()
	logLevel = logging.WARN
	if opts.verbose:
		logLevel = logging.INFO
	utils.Logger(logLevel, opts.log)
	logger = logging.getLogger()
	if len(args) < 4:
		logger.error("too few args")
		parser.print_help()
		sys.exit(1)
	try:
		window = int(args[0])
	except Exception as e:
		logger.exception(e)
		parser.print_help()
		sys.exit(1)
	hitFileName = args[1]
	genomeFileName = args[2]
	outputFileName = args[3]
	try:
		hitFile = open(hitFileName)
		outputFile = open(outputFileName, "w")
	except Exception as e:
		logger.exception(e)
		parser.print_help()
		sys.exit(1)
	hashes = utils.FastaParser(genomeFileName).hashMap()
	logger.info("Collected %d genomes from input %s" % (len(hashes.keys()), genomeFileName))
	dataRegex = re.compile(r'query:(\S+)\s+qseq:(\S+)\s+hit:(\S+)\s+sense:(\S+)\s+beg:(\d+)\s+end:(\d+)\s+hseq:(\S+)')
	for line in hitFile:
		line = line.strip()
		match = dataRegex.search(line)
		if match:
			try:
				id = match.group(1)
				chr = match.group(3)
				sense = match.group(4)
				x = int(match.group(5))
				y = int(match.group(6))
				hseq = match.group(7)
			except Exception as e:
				logger.error(e)
				continue
			chr = chr.replace('>', '').strip()
			if sense == 'sense':
				beg = max(1, x - window)
				end = min(len(hashes[chr]), y + window)
				s = "+"
				pos = min(window, x-1)
			else:
				beg = min(len(hashes[chr]), x + window)
				end = max(1, y - window)
				s = "-"
				pos = min(window, len(hashes[chr]) - x -1)
			name = "%s_%s_%d%s" %(id, chr, x, s)
			try:
				subseq = utils.getSubseq(hashes[chr], beg, end)
			except Exception as e:
				logger.error(e)
				continue
			start = pos
			stop = start + abs(x - y)
			print(">%s\t%s\t%d\t%d\t%s\t%d\t%d\n%s" %(name, chr, beg, end, hseq, start, stop, subseq), file=outputFile)
	outputFile.close()
	hitFile.close()

if __name__ == "__main__":
	main()

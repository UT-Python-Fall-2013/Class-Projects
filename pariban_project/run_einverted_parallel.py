#!/bin/env python3
import sys, os, time, argparse, logging, re, signal
from utils import *
from multiprocessing import Pool

def init_worker():
	signal.signal(signal.SIGINT, signal.SIG_IGN)  # use this to prevent hanging on a ctrl-c. source: http://noswap.com/blog/python-multiprocessing-keyboardinterrupt

def runEinverted(x_uniq, x_key, x_subseq, x_gap, x_threshold, x_match, x_mismatch, x_dist, x_x, x_y):
	logger = logging.getLogger() # get logger
	fileName = os.path.join("/tmp", str(x_uniq)+x_key) # create temp file name
	try:
		file = open(fileName, "w")
	except Exception as e:
	 	logger.exception(e)
	 	return None
	file.write(">%s\n%s\n" % (x_key, x_subseq)) # write sequence into the file
	file.close() 
	# now launch einverted (emboss) 
	os.system("einverted -sequence %s -outfile %s -gap %d -threshold %d -match %d -mismatch -%d -maxrepeat %d -auto" % (fileName, fileName+"_out", x_gap, x_threshold, x_match, x_mismatch, x_dist))
	try:
		# output written to a temporary file
		output = open(fileName + "_tmp", "w")
		file = open(fileName+"_out")
		output.write(">%s %d %d\n" % (x_key, x_x, x_y))
		output.write(file.read())
		file.close()
		output.close()
	except Exception as e:
	 	logger.exception(e)
	 	return None
	if x_uniq % 1000 == 0:
		logger.info("+1K Complete: %s %d %d" % (x_key, x_x, x_y))
	logger.debug(fileName+"_tmp")
	return fileName+"_tmp" # return the name of the file

def parseEinverted (x_input, x_output, x_minPct, x_minArm):
	logger = logging.getLogger() # get logger
	try:
		input = open(x_input)
	except Exception as e:
		logger.exception(e)
		return
	data = {}
	lineCnt = 5
	headerRegex = re.compile(r'^>(?P<CHR>\S+)\s+(?P<POS>\d+)\s+(?P<LEN>\d+)$') # regular expression, copied from miRcheck
	scoreRegex = re.compile(r'Score\s*(\d+):\s*(\d+)\/(\d+)')
	dataRegex = re.compile(r'^(\d+)\s+\S+\s+(\d+)')
	chr = None
	pos = None
	counter = 0
	for line in input:
		counter += 1
		if counter % 100000 == 0:
	 		logger.info("Processed %dK lines" %(counter/1000))
		lineCnt += 1
		line = line.strip() # strip newline
		match = headerRegex.search(line) # check for header
		if match:
			chr = match.group('CHR') # get chr
			pos = int(match.group('POS')) # get pos
			logger.debug("found chr = %s pos = %d" % (chr, pos))
			continue
		match = scoreRegex.search(line) # check if score match found
		if match:
			score = match.group(1)
			try:
				ratio = float(match.group(2))/float(match.group(3)) # get match ratio
				if ratio >= x_minPct:
					lineCnt = 0
					logger.debug("found match with ratio %f" % ratio)
			except Exception as e:
				logger.exception(e)
				continue
		match = dataRegex.search(line) # check if line contains data
		if match and lineCnt == 1:
			try:
				logger.debug("begin:found data match with %d %d" % (int(match.group(1)), int(match.group(2))))
				beg = int(match.group(1)) + pos -1
				end = int(match.group(2)) + pos -1
			except Exception as e:
				logger.exception(e)
				continue
		if match and lineCnt == 3:
			try:
				logger.debug("end:found data match with %d %d" % (int(match.group(1)), int(match.group(2))))
				if ((end - beg >= x_minArm) and (int(match.group(1)) - int(match.group(2)) >= x_minArm)):
					beg2 = int(match.group(1)) + pos -1
					end2 = int(match.group(2)) + pos -1;
					data["%s %d %d %d %d %s %f" %(chr, beg, end, beg2, end2, score, ratio)] = 1
			except Exception as e:
				logger.exception(e)
				continue
	try:
		output = open(x_output, "w")
	except Exception as e:
		logger.exception(e)
		return
	for line in sorted(data): # sort data and then print into line
		print(line, file=output)
	output.close()

def main():
	# create a parser
	usage=("usage: %s <input file> <output file>" % sys.argv[0])
	parser = argparse.ArgumentParser(description='A wrapper to run_einverted', usage=usage)

	# Add optional switches
	parser.add_argument('-v', '--verbose', action='store_true', help='produce verbose output [%(default)s]')
	parser.add_argument('-l', '--log', action='store', dest='log', default=None, help='dump messages in log file [%(default)s]')
	parser.add_argument('-g', '--gap', action='store', type=int, default=6, help='Please add description [%(default)s]')
	parser.add_argument('-t', '--threshold', action='store', type=int, default=40, help='Please add description [%(default)s]')
	parser.add_argument('-m', '--match', action='store', type=int, default=3, help='Please add description [%(default)s]')
	parser.add_argument('-M', '--mismatch', action='store', type=int, default=3, help='Please add description [%(default)s]')
	parser.add_argument('-d', '--dist', action='store', type=int, default=240, help='Please add description [%(default)s]')
	parser.add_argument('-w', '--window', action='store', type=int, default=2000, help='size of fragments sent to einverted [%(default)s]')
	parser.add_argument('-s', '--step', action='store', type=int, default=1000, help='overlap of fragments sent to einveted [%(default)s]')
	parser.add_argument('-p', '--min_pct', action='store', type=float, default=0.35, help='minimum percent complementarity in inverted repeats [%(default)s]')
	parser.add_argument('-a', '--min_arm', action='store', type=float, default=15, help='minimum size in inverted repeats [%(default)s]')
	parser.add_argument('-D', '--degree', action='store', type=int, default=10, help='Degree of parallelism [%(default)s]')
	opts,args = parser.parse_known_args()

	logLevel = logging.WARN
	if opts.verbose:
		logLevel = logging.INFO
	Logger(logLevel, opts.log)
	logger = logging.getLogger()

	# below are parameters passed to einverted 
	gap = opts.gap
	threshold = opts.threshold
	match = opts.match
	mismatch = opts.mismatch
	dist = opts.dist
	win = opts.window  # size of fragments sent to einverted
	step = opts.step # overlap of fragments sent to einveted
	minPct = opts.min_pct  #mimimum percent complementarity in inverted repeats        
	minArm = opts.min_arm #minimum size in inverted repeat

	if len(args) < 2:
		logger.error("Too few args")
		parser.print_help()
		sys.exit(1)

	inputFileName = args[0]
	outputFileName = args[1]
	
	# create a pool of workers
	pool = Pool(opts.degree, init_worker)

	# read fasta file
	try:
		logger.info("creating hashes")
		hashes = FastaParser(inputFileName).hashMap()
	except Exception as e:
		logger.exception(e)
		pool.close()
		pool.terminate()
		pool.join()

	results = []
	uniq = 0
	try:
		for key in sorted(hashes):
			val = hashes[key]
			size = len(val)
			logger.info("%s %d" % (key, size))
			for x in range(1, (size - (win - step)+1), step):
				uniq += 1
				y = x + win -1
				subseq = getSubseq(val, x, y) # get subsequence from hash
				result = pool.apply_async(runEinverted, [uniq, key, subseq, gap, threshold, match, mismatch, dist, x, y]) # pass the data to pool of workers to run
				results.append(result) # append the result
		pool.close()
		pool.join()
		logger.info("Finished running Einverted")
	except KeyboardInterrupt:
		logger.warn("Received KeyboardInterrupt")
		pool.close()
		pool.terminate()
		pool.join()
	logger.info("Processing output files")
	output = open(outputFileName+"_tmp", "w") # open output file
	counter = 0
	for result in results:
	 	counter += 1
	 	if counter % 1000 == 0:
	 		logger.info("Processed %dK files" %(counter/1000))
	 	if not result.ready():
	 		logger.warn("result not ready")
	 		continue
	 	if not result.successful(): 
	 		logger.warn("result not successful")
	 		continue
	 	fileName = result.get()
	 	if not fileName:
	 		logger.warn("output none")
	 		continue
	 	logger.debug("processing file: %s" % fileName)
	 	file = open(fileName)
	 	output.write(file.read())
	 	file.close()
	output.close()
	logger.info("Parsing temp file to generate final output")
	parseEinverted(outputFileName+"_tmp", outputFileName, minPct, minArm) # parse einverted output

if __name__ == "__main__":
	main()


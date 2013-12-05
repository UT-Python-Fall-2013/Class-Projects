# contains some utility classes
import os,sys,traceback,logging

class Assert:
  def __init__ (self, x_eval, x_msg, x_file=sys.stderr):
    if not x_eval: # evaluate the expression and if it evaluates to false
      x_file.write(x_msg + '\n') # write the message
      traceback.print_exc(file=x_file) # print backtrace
      sys.exit(1) # exit

class FastaParser:
  def __init__ (self, x_file):
    Assert(os.path.isfile(x_file), "ERROR: File '%s' not found" % x_file) # check that the file exists
    self.sequences = {} # create a map to store the sequence
    with open(x_file) as fileHandler: # open the fasta file
      header = ""
      for line in fileHandler: # for all lines
        line = line.strip() # get rid of newline
        if line.startswith(">"): # header line
          header = line.lstrip(">") # get rid of starting ">"
          self.sequences[header] = "" # initialize the sequence
        elif line.startswith("#"): # comment lines, ignore
          pass
        else:
          Assert(len(header), "ERROR: reading fasta file "+x_file) # header should be set
          self.sequences[header] += line # join the line to already read sequence

  def hashMap (self): # return the map itself
    return self.sequences

  def getHeaders (self): # get all headers
    return self.sequences.keys()

  def getSequences (self): # get all sequences
    return self.sequences.values()
  
  def getSequence (self, x_name): # get a particular sequence
    return self.sequences[x_name]

  def __iter__ (self): # iterator
    return self.sequences.__iter__()

  def next (self): # so that we can call next() on the parser, it will return the next one until end
    return self.sequences.next()
  
class c_Formatter( logging.Formatter ): # required to setup logger
    """Extension to Formatter class for customized print messages"""
    outputFormat = "%(message)s" #create formatter
    def __init__(self, fmt="%(levelname)s: %(message)s"):
        logging.Formatter.__init__(self, fmt) # set formatter for logging

    def format(self, record):
        # Save the original format configured by the user when the logger
        # formatter was instantiated
        l_origFormat = self._fmt
        # Replace the original format with one customized by logging level
        if record.levelno > logging.CRITICAL or record.levelno < logging.WARNING: # Print only critical, error and warning messages
            self._fmt = c_Formatter.outputFormat
        # Call the original formatter class to do the grunt work
        l_result = logging.Formatter.format(self, record)
        # Restore the original format configured by the user
        self._fmt = l_origFormat # back up old format
        return l_result

class Logger:
  def __init__ (self, x_level, x_logfile=None, x_name = None):
    """Sets up logging object"""
    # get a logger and set its level
    l_logger = logging.getLogger(x_name) if x_name else logging.getLogger() # if a name is provided, create the logger, else get root logger
    l_logger.setLevel(x_level) # set logging level
    l_logger.propagate = False # parent logger does not replicate the message

    l_formatter = c_Formatter() #create formatter
    # if a log file was specified, add a handler
    if x_logfile:
        l_fh = logging.FileHandler(x_logfile, "w") # open logfile
        l_fh.setLevel(x_level) # set level to the file handler
        l_fh.setFormatter(l_formatter) #add formatter to file handler
        l_logger.addHandler(l_fh) # add the file handle to logger

    l_ch = logging.StreamHandler() # output also to terminal
    l_ch.setLevel(x_level) # set level
    l_ch.setFormatter(l_formatter) #add formatter to ch
    l_logger.addHandler(l_ch) #add ch to logger
    logging.__dict__['OUTPUT'] = sys.maxsize # add a new level OUTPUT in logger that is higher than all existing levels

def revcomp(x_val):
  seq = x_val[::-1] # reverse
  seq = seq.translate(bytes.maketrans(b'Uu', b'Tt')) # replace U with T and u with t
  seq = seq.translate(bytes.maketrans(b'ACGTacgt', b'TGCAtgca')) # repace A->T,C->G, etc.
  return seq

def getSubseq(x_val, x_x, x_y): # gets the sub sequence 
  subseq = None
  if (x_y >= x_x): # if end >= start
    beg = max(1,x_x)
    end = min(x_y,len(x_val))
    subseq = x_val[beg-1:end]
  else:
    beg = min(len(x_val),x_x)
    end = max(x_y,1)
    subseq = x_val[end-1:beg]
    subseq = revcomp(subseq)
  return subseq


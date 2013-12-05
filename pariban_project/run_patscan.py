#!/bin/env python3
import os, sys, re, subprocess, argparse

g_packageDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(g_packageDir)
from setupEnv import PackageRequirements
from utils import Assert,FastaParser

class RunPatScan:
  def __init__ (self, x_fastaFile, x_db, x_mismatches, x_deletions, x_insertions, x_outfile):
    req = PackageRequirements(True)
    self.fastaFile = x_fastaFile
    self.db = x_db
    self.mismatches = x_mismatches
    self.deletions = x_deletions
    self.insertions = x_insertions
    self.outfile = x_outfile
    Assert(self.outfile, "Please provide an outfile")
    Assert(self.db and os.path.isfile(self.db), "Please provide a reference genome file")
    self.fastaData = FastaParser(self.fastaFile)

    hitSearcher = re.compile(r'^>(?P<chr>\S+):\[(?P<beg>\d+),(?P<end>\d+)\]')
    with open(self.outfile, "w") as outFileHandler:
      for name in iter(self.fastaData):
        with open(self.outfile+"_pattern", "w") as patFileHandler:
          print("%s[%d,%d,%d]" % (self.fastaData.getSequence(name), self.mismatches, self.deletions, self.insertions), file=patFileHandler )
        with open(self.db) as dbHandler , open(self.outfile+"_hits", "w") as hitFile:
          subprocess.check_call(["scan_for_matches", "-c", self.outfile+"_pattern"], stdin=dbHandler, stdout=hitFile, stderr=subprocess.STDOUT)
          print("scan_for_matches -c "+name+"_pattern", "<", self.db, ">", name+"_hits")
        with open(self.outfile+"_hits") as hitFile:
          for line in hitFile:
            line = line.strip()
            if line.startswith(">"):
              matched = hitSearcher.match(line)
              if matched:
                chr = matched.group('chr')
                end = int(matched.group('end'))
                beg = int(matched.group('beg'))
                if end > beg:
                  sense = "sense"
                else:
                  sense = "antisense"
            else:
              outFileHandler.write("query:%s qseq:%s hit:%s sense:%s beg:%d end:%d hseq:%s\n" % (name, self.fastaData.getSequence(name), chr, sense, beg, end, line))


if __name__ == "__main__":
  descr = "FIXME ********** Description of run_patscan"
  parser = argparse.ArgumentParser(description=descr)
  parser.add_argument("-p", "--pat", metavar="FILE", help="pattern fasta file", default="", required=True)
  parser.add_argument("-r", "--db", metavar="FILE", help="reference genome fasta file", default="", required=True)
  parser.add_argument("-o", "--out", metavar="FILE", help="output hit file", default="", required=True)
  parser.add_argument("-m", "--mismatches", metavar="NUM", help="number of allowed mismatches", type=int, default=0)
  parser.add_argument("-d", "--deletions", metavar="NUM", help="number of allowed deletions", type=int, default=0)
  parser.add_argument("-i", "--insertions", metavar="NUM", help="number of allowed insertions", type=int, default=0)
  args = parser.parse_args()
  RunPatScan(args.pat, args.db, args.mismatches, args.deletions, args.insertions, args.out)

### The script separate R1 and R2 pair-end reads into two files
### command line: "python separate_reads.py -h" will print the usage

import argparse
import fileinput

parser = argparse.ArgumentParser()
parser.add_argument("-r1","--read1", type=argparse.FileType('w'),help="r1 reads will write to this file")
parser.add_argument("-r2","--read2", type=argparse.FileType('w'),help="r2 reads will write to this file")
    
parser.add_argument("-in","--file",help="the mixed read file to be separated")
args = parser.parse_args()


group = 4   # four lines per reads in fastq format
count = 0
write_to_r1 = 1 

for line in fileinput.input(args.file):     # loop through the merged read file
    if write_to_r1:                         # if write_to_r1 is True
        args.read1.write(line)              # write the line into R1 file
        count +=1                   
        if count == group:                  # after loop though the fourth line
            write_to_r1=0                   # set write_to_r1 False
            count = 0                       # reset the accumulator
        else:
            pass
    elif not write_to_r1:                   # if write_to_r1 is False
        args.read2.write(line)              # # write the line into R2 file
        count +=1
        if count == group:                  # after loop though the fourth line
            write_to_r1 =1                  # set write_to_r1 False
            count = 0                       # reset the accumulator
        else:
            pass

args.read2.close()
args.read1.close()


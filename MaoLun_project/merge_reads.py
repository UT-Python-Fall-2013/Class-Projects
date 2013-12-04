
### This script merge two pair-end reads into one fastq file 


import argparse
import fileinput


parser = argparse.ArgumentParser()
parser.add_argument("-r1","--read1",help="r1 reads file to open")
parser.add_argument("-r2","--read2",help="r2 reads file to open")
    
parser.add_argument("-out","--file", help="the mixed read will write to this file")
args = parser.parse_args()

    
with open(args.read2) as r2, open (args.file,"w") as out:
    N =4
    cout =1 
    
    for line in fileinput.input(args.read1):   # loop through R1 read file
        
        if cout <N:                            # when accumulator smaller than 4 
            cout += 1
            out.write(line)                    # print line into merged file
        
        elif cout == N:                        # when loop to fourth line 
            out.write(line)                    # print the fourth line of R1 file into merged file
            out.write(r2.readline())           # print first line of R2 file into merged file
            out.write(r2.readline())           # print second line of R2 file into merged file
            out.write(r2.readline())           # print third line of R2 file into merged file
            out.write(r2.readline())           # print fourth line of R2 file into merged file
            cout =1                            # reset the accumulater
            

r2.close()
out.close()

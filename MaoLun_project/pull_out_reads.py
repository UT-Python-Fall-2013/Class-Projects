#! /usr/bin/env python


import subprocess
import argparse
import random


def parse_arg():
    """parse out the arguments from commend line"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-in","--in_file",help="enter the r1 to open")
    parser.add_argument("-out","--out_file", help="enter the out1 to write")
    parser.add_argument("-t","--type", help="PE:pair-end reads, or SE:single-end read")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-p","--percent",type=float,help="percentage of reads to pull out")
    group.add_argument("-n", "--number", type=int,help="the number of reads to pull out")
    args = parser.parse_args()
    
        
    return ([args.in_file,args.out_file,args.type,args.percent,args.number])

def draw_reads(N,infile,outfile,percentage):
    """radomly pull out certain percent of reads from a fastq file"""
    L=[]                             # a list to store N lines from the read file    
    for x in range (N):            
        L.append(infile.readline())  # append N lines in a list

    if random.random() < percentage: # if random number smaller than the percentage specified by user
        for i in range(len(L)):      # write the N lines into output file  
            outfile.write(L[i])
            
    else:                            # if random number not smaller than the percentage specified by user
        L1 =[]                       # empty the list     

def count_line(Infile):
    # subprocess bash command 'wc -l' to get the total number of lines of the read file
    process = subprocess.Popen(['wc','-l',Infile],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr = process.communicate()
    length = stdout.split()[0]        # total number of lines in read file
    
    return (length)

def main():
    
  
    Infile,Outfile,read_type,percent,number = parse_arg()    # get arguments 
    
    if read_type == 'PE':      # eight lines per pair-end read in fastq format
        N = 8            
    elif read_type == 'SE':   # four lines per single-end read in fastq format
        N = 4               
    

    length = count_line(Infile) # get total number of lines in the read file

    with open (Outfile,"w") as out,open(Infile,"r") as r1:
        if percent:
            for i in range(0,int(length),N):    # loop over the read file every N line
                draw_reads(N,r1,out,percent)    # draw random reads
             
        elif number:
            ratio = float(number*N)/int(length) # ratio = (the number of reads * N, which are lines needed)/total lines
            for i in range(0,int(length),N):    # loop over the read file every N line
                draw_reads(N,r1,out,ratio)      # draw radon reads
 
    out.close()
    r1.close()
main()



"""This script take a reference fasta file and a bam file
   to plot the coverage along the reference"""

import pysam
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse

def parse_arg():
    """parse out the arguments from commend line"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-in","--fasta",help="enter the fasta file to open")
    parser.add_argument("-b","--bam",help="enter the bam file")
    parser.add_argument("-w","--window", type=int,
                        choices=([100,200,500,1000,2000]),
                        help="Scan window size in bp")
    args = parser.parse_args()
    
    return (args.fasta,args.bam,args.window)


def Fasta_Length_and_N(infile):
    """Get the length of fasta seq, the number of N, and the position of N's"""
    N_count = 0
    N_position = []
    for record in SeqIO.parse(open(infile),"fasta"):
        seq_length = len(record.seq)
        reference = record.id
        for position, nul in enumerate(record.seq):
            if nul =="N":
                N_position.append(position)
                N_count+=1
    
    return (reference, seq_length,N_count, N_position)



def Mapped_count(reference,ref_length,bam,window):
    """Use pysam module to count mapped reads along a reference in a given scan window size
       It estimate coverage for reads that both pair-end reads mapped onto the reference"""
    
    samfile = pysam.Samfile(bam,"rb")
    count_list=[]
    for x in range(1,ref_length,window):
        count=0
        for alignedread in samfile.fetch(reference,x,x+window):
            if alignedread.is_paired:
                count += 1
        count_list.append(count)
        count = 0
    
    return count_list


def Coverage_plot(ref_length,window,count_list,N_position):
    """Use matplotlib.pyplot to plot coverage"""
    
    Avg = sum(count_list)/len(count_list)

    p1,=plt.plot(range(1,ref_length,window),count_list)     # plot the coverage 
    p2,=plt.plot(N_position,[Avg]*len(N_position),'ro')     # plot the N's (gaps)
    plt.legend([p1,p2],[str(window)+'bp',"Gap"])    
    
    plt.show()

def main():
    
    # get arguments from command line
    fasta,bam,window = parse_arg()
    
    # get reference seq ID, reference length, the number of N's, and the position of these N's
    reference,ref_length, N_count, N_position = Fasta_Length_and_N(fasta)
    
    # get a list of mapped read counts bin by a given window size
    count_list = Mapped_count(reference,ref_length,bam,window)
    
    # plot the coverage
    Coverage_plot(ref_length,window,count_list,N_position)
    
    
main()
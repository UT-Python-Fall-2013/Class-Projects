#! /usr/bin/env python

"""The velvet pipeline"""


import subprocess
import os
import argparse
import tempfile
from Bio import SeqIO


def parse_arg():
    """parse out the arguments from commend line"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-o","--out_dir",help="The velvet output dir")
    parser.add_argument("-f","--read",help="The merged read file name without extension")
    parser.add_argument("-k","--kmer",type=int,help="Starting kmer size")
    parser.add_argument("-in","--increment",type=int,help="kmer size increment")
    parser.add_argument("-N","--times",type=int,help="Number of increments")
    parser.add_argument("-c","--Coverage",type=int,help="Expected coverage")
    parser.add_argument("-g","--genome", choices=(["CP","MT"]),help="CP:plastid contigs, or MT:mitoshondrial contigs")
    
    args = parser.parse_args()
    
    return (args.out_dir,args.read,args.kmer,args.increment,args.times,args.Coverage,args.genome)

def velveth_args(read,kmer_start,increment,N):
    """"Generate velveth arguments without specifing an output folder"""
    kmer_end = kmer_start + (increment * N)    
    velveth_args_list = ['velveth',str(kmer_start)+','+str(kmer_end)+','+str(increment),
                         '-fastq','-shortPaired',read]
    return velveth_args_list

    # velveth_args_list.insert(1,output_folder)
    
def velvetg_args(exp_cov,scaf='yes'):
    """Generate velvetg arguments without specifing an output folder"""
    cov_cutoff = exp_cov / 10
    max_coverage = exp_cov * 10
    
    velvetg_args_list=['velvetg','-cov_cutoff', str(cov_cutoff), '-exp_cov',str(exp_cov),
                       '-max_coverage',str(max_coverage), '-scaffolding', scaf]
    return velvetg_args_list
    # velvetg_args_list.insert(1,output_folder+'/_'+i)
    
def Run(arg_list):
    """subprocess given a list of arguments"""
    process = subprocess.Popen(arg_list,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr = process.communicate()
    return (stdout)
    
def blast_arg(query, db, out_format = '6',program='blastn'):
    """Generate arguments for blast"""
    blast_arg_list = [program,'-query',query, '-db',db,'-outfmt',out_format]
    return blast_arg_list


def blast_hit(blast_arg_list):
    """return a list of blast hits"""
    stdout = Run(blast_arg_list)     # Run blast

    temp = tempfile.TemporaryFile()  # write blast ouput into a temp file
    temp.write(stdout)                          
    temp.seek(0)
    
    hit_list=[]                      # extract the first column from the blast output,
    for line in temp:                # which are the hits
        ele = line.split('\t')
        if ele[0] not in hit_list:
            hit_list.append(ele[0])        
    temp.close()
    
    return hit_list

def fasta(hit_list,seq,genome):
    """Extract blast hits and write to new files"""
    if genome == "CP":
        file_name = 'CP'+ seq
    elif genome == "MT":
        file_name = 'MT'+ seq
    output_file = open(file_name,"w")
    fa = SeqIO.to_dict(SeqIO.parse(open(str(seq),'r'),"fasta"))
    for j in hit_list:
        SeqIO.write(fa[str(j)],output_file,"fasta")


def countN(fasta):
    """count the number of N's and the lenght of a contig"""
    fa = SeqIO.to_dict(SeqIO.parse(open(str(fasta),'r'),"fasta"))
    info =[]
    for i in fa.keys():
        count = 0
        length = len(fa[i])
        for j in fa[i]:
            if j =='N':
                count+=1
        info.append([i,count,length])
    return info

def main():
    
    # Get arguments from command line
    out_dir, read_name, Kmer_start, increment, N, Cov,genome = parse_arg()
    read = read_name+".fq"
    
    # set up path for blast db 
    if genome == 'CP':
        blastdb = '../../../cp_blast/vex'   # relative path to blastdb 
    elif genome == 'MT':
        blastdb = '../../../cp_blast/vex'
    
    # set up CP or MT contigs folder and seq names
    Genome_contigs_folder = genome+"contigs"
    Genome_contigs_name = genome+read_name
    
    
    ########## Run velveth
    velveth_args_list = velveth_args(read,Kmer_start,increment,N)   # generate arguments
    os.mkdir(out_dir)                                               # mkdir output folder
    velveth_args_list.insert(1,out_dir+'/')                         # insert output folder into argument list
    Run(velveth_args_list)                                          # run velvet
    
    
    ######### Run velvetg in each kmer folder
    for i in range(Kmer_start,Kmer_start+increment*N,increment):    # loop through kmer folders
        velvetg_args_list = velvetg_args (Cov)                      # generate arguments
        velvetg_args_list.insert(1,out_dir+'/_'+str(i))             # insert output folder and kmer folder into argument list
        Run(velvetg_args_list)                                      # run velvetg
    
    
    
    ######### Copy contigs from each kmer folder into contigs folder 
    os.mkdir(out_dir+"/contigs")
    for i in range(Kmer_start,Kmer_start+increment*N,increment):
        src = "kmer/_"+str(i)+"/contigs.fa"
        dst = "kmer/contigs/"+read_name+"Cov"+str(Cov)+"K"+str(i)+".fa"
        os.rename(src ,dst)         # the same as "mv" in bash
    
    
    
    ######## Run blast and write hitted contigs into new folder
    contigs = os.listdir("kmer/contigs")
    os.chdir("kmer/contigs")
    
    for i in contigs:                               # loop though contigs.fa in contigs folder
        blast_arg_list = blast_arg(i,blastdb)       # generate blast argument for each contig file
        hit= blast_hit(blast_arg_list)              # get a list of hitted sequence for each contig file
        fasta(hit,i,genome)                         # write hitted sequences into new organelle contig files
     
    
    
    ######### Copy organelle contigs to organelle contigs folder
    os.chdir("..")
    os.mkdir(Genome_contigs_folder)
    for i in range(Kmer_start,Kmer_start+increment*N,increment):
        src = "contigs/"+Genome_contigs_name+"Cov"+str(Cov)+"K"+str(i)+".fa"
        dst = Genome_contigs_folder+"/"+Genome_contigs_name+"Cov"+str(Cov)+"K"+str(i)+".fa"
        os.rename(src ,dst)
    
    
    ######## For each contig file print seq name, number of N's, seq length
    Cpcontigs_list = os.listdir(Genome_contigs_folder)
    os.chdir(Genome_contigs_folder)
    for x in Cpcontigs_list:
        information = countN(x)
        count_N=[]
        length_list=[]
        name_list=[]
        print (x)
        print ("contig\tN\t\tLength")
        for i in range(len(information)):
            count_N.append(information[i][1])
            length_list.append(information[i][2])
            print ("{0}\t{1}\t{2}".format(information[i][0],information[i][1],information[i][2]))

    

main()
        
    
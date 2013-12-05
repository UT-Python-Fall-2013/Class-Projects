./run_patscan.py -p data/sample_miRNAs.fa -r data/sample_genome.fa -m 2 -d 0 -i 0 -o sample_miRNA_matches

./retrieve_genomic_regions.py 10 sample_miRNA_matches data/sample_genome.fa sample_miRNA_matches_genomic

./run_einverted_parallel.py data/sample_genome.fa sample_IRs                                                                                                                                         
./fold_inverted_repeats.py sample_IRs data/sample_genome.fa sample_IRs_f 

# Usage 
#usage: run_patscan.py [-h] -p FILE -r FILE -o FILE [-m NUM] [-d NUM] [-i NUM]
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -p FILE, --pat FILE   pattern fasta file
#  -r FILE, --db FILE    reference genome fasta file
#  -o FILE, --out FILE   output hit file
#  -m NUM, --mismatches NUM
#                        number of allowed mismatches
#  -d NUM, --deletions NUM
#                        number of allowed deletions
#  -i NUM, --insertions NUM
#                        number of allowed insertions

# usage: ./retrieve_genomic_regions.py <window size> <miRNA match file> <input genome file> <output file>

# usage: run_einverted_parallel.py sample_genome.fa sample_IRs                                                                                                                                         
# usage: fold_inverted_repeats.pl sample_IRs sample_genome.fa sample_IRs_f  

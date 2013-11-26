# Annotates a genomic contig, directly from an assembler
# such as 'velvet'

from Bio import SeqFeature
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
import os



# Name of the contig as input
contig_name = input('Please enter the name of the contig: ')



# Run Glimmer3 on the contig to predict the best ORFs, it will
# produce several files but the one that contains the coordinates
# of the predicted genes is 'prefix.predict'
prefix = contig_name.split('.')[0]

print('\n##### PREDICTING GENES WITH GLIMMER 3.02 #####\n')

command = 'g3-iterated.csh ' + contig_name + ' ' + prefix
os.system(command)

print('\n##### PREDICTION FINISHED #####\n')



# The dictionary 'annotations' will contain as keys 
# the names of the ORFs and as values the start and
# end positions, the strand, and a field for the
# name of the best BLAST hit
annotations = {}

predict_file = prefix + '.predict'
file_in = open(predict_file, 'r')
for line in file_in:
    if line[0] != '>':
        items = line.strip('\n').split()
        annotations[items[0]] = items[1:3]
        if int(items[3]) < 0:
            annotations[items[0]].append('-1')
        else:
            annotations[items[0]].append('+1')
        annotations[items[0]].append('NO PREDICTION')
file_in.close()        



# Create a multiple-sequence fasta with the ORFs, based on the coordinates
# given by the Glimmer prediction, using 'extract' tool that is part of Glimmer
orfs_file = prefix + '_orfs.fa'

command = 'extract ' + contig_name + ' ' + predict_file + ' >' + orfs_file
os.system(command)

print('\n##### EXTRACTION OF ORFs FINISHED #####')



# BLASTX using as queries the sequences in _orfs.fa file, to target proteins 

print('\n##### BLASTX STARTED #####\n')

file_in = open(orfs_file, 'rU')
queries = list(SeqIO.parse(file_in, 'fasta', ))
file_in.close()
for query in queries:
    result_handle = NCBIWWW.qblast('blastx', 'nr', query.format('fasta'), hitlist_size=5)
    xml_out = query.id + '.xml'
    print('Saving {0}...'.format(xml_out))
    save_xml = open(xml_out, 'w')
    save_xml.write(result_handle.read())
    save_xml.close()
    result_handle.close()

print('\n##### BLASTX DONE #####\n')



# Parse the xmls from BLASTX output
all_files = os.listdir('.')
all_xmls = []
for file in all_files:
    if file.endswith('.xml'):
        all_xmls.append(file)

for xml_file in all_xmls:
    result_handle = open(xml_file)
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
    orf_name = xml_file.strip('.xml')
    if len(blast_record.alignments) > 0:
        gene_name = blast_record.alignments[0].title.split('|')[4].strip('>gi').strip()
        annotations[orf_name][3] = gene_name



# Create a SeqRecord from the contig, so we can add features (annotations)
handle = open(contig_name, 'rU')
contig = SeqIO.read(handle, 'fasta', alphabet=IUPAC.ambiguous_dna)
handle.close()



# Add annotations to the contig from the ORFs predicted by glimmer3.02 and from
# their BLASTX results
for orf in annotations:
    start = int(annotations[orf][0])
    end = int(annotations[orf][1])
    orf_strand = int(annotations[orf][2])
    feature_type = 'ORF Glimmer'
    new_feature = SeqFeature(FeatureLocation(start, end), type=feature_type, strand=orf_strand)
    new_feature.qualifiers['locus_tag']=orf
    new_feature.qualifiers['gene']=annotations[orf][3]
    contig.features.append(new_feature)



# Write the annotated sequence in GenBank format
annotated_name = prefix + '_annotated.gb'
SeqIO.write(contig, annotated_name, 'genbank')

print('\n##### ANNOTATION FINISHED, THE FILENAME IS {0} #####\n'.format(annotated_name))

##########################################
#
# This program finds putative transcription
# factor binding sites in the sequence of
# the promiscuous ancestor and monogamous
# ancestor using regular expressions and
# some biopython commands- The input is an
# alignment of New world voles and the output
# is a matrix of transcription factor binding
# sites. The program also informs the position
# respect the first nucleotide.
#
#       By: Alejandro Berrio
#       For Python class
#
###########################################


def openfasta(fasta):   # This module helps to parse a fasta alignment

    from Bio import AlignIO
    file = input("Enter the name of your fasta alignement file (please indicate extention): ")   # To ask for input

    try:    #to handle errors in opening the file

        for seqs in AlignIO.parse(file, "fasta"):   #for parsing fasta file
                alignment = seqs                    #to convert sequence in string

    except NameError:
        print("\nThe file doesn't seem to exist in home folder")            # Error excemption

    return alignment

def sequencing(sequencia):                              # This module make sure that the sequences are correct style
    dna = []
    dna1= ""
    for i in sequencia:
        for x in i:
            if x == 'A':
                dna.append(x)
            elif x == 'G':
                dna.append(x)
            elif x == 'C':
                dna.append(x)
            elif x == 'T':
                dna.append(x)
    
     
    return dna

def reversing(secuencia):                               # this formula reverses the sequence because Regulatory elements can also be in minus
    secuencia.reverse()

    return secuencia

# this is a formula to find a the complementary sequence
def complement(secuencia):
    cdna = []
    for i in secuencia:
        if i == 'A':
            x = 'T'
            cdna.append(x)
        elif i == 'T':
            x = 'A'
            cdna.append(x)
        elif i == 'G':
            x = 'C'
            cdna.append(x)
        elif i == 'C':
            x = 'G'
            cdna.append(x)

    return cdna

def transfac(secuencia):            # this module finds transcription factor binding sites using regular expressions
    import re
    count = 0
    ## just for reference, this is the ASCII code for DNA sequence ambiguities
        #K = [G|T]
        #M = [A|C]
        #R = [A|G]
        #Y = [C|T]
        #S = [C|G]
        #W = [A|T]
        #B = [C|G|T]
        #V = [A|C|G]
        #H = [A|C|T]
        #D = [A|G|T]
        #N = [A|C|G|T]
    ForwardTFBSs = []           #output forward
    ReverseTFBSs = []           #output reverse


    #Dictionary of transcription factors that are assumed to be highly expressed in the ventral pallidum vs Striatum

    
    TFBS = {"PU1" : r"[A|C|G|T][A|C|G|T][A|T]A[A|G|T][C|G][A|C|G]GGAA[C|G]T[A|C|G|T][A|C|G|T][A|C|G|T][A|C|G|T]",
        "STAT3" : r"[A|C|G|T][A|C|G|T][A|C|G|T]TTCCC[A|C|G|T]",
        "PAX5" : r"[A|G|T][A|G][A|C|G|T][C|G|T][A|C|T][A|C|T][A|C|G|T][C|T][C|G][A|G|T][A|T]G[A|C]G[G|T][A|G][A|G]C[C|G][A|G]",
        "GABP" : r"[A|C|G]CCGGAAG[A|C|G|T]G[C|G][A|G]",
        "FOXM1" : r"A[A|G][A|C|T]T[A|G|T]GA[C|G|T]T",
        "EGRI" : r"[A|C|T]GCGTGGG[C|T]G[G|T]",
        "E2F1" : r"[A|C|T][A|C|G|T]TTTC[A|C|T][A|C|G|T]",
        "MEF2" : r"[A|G|T]G[C|T]TAT[A|T]TT[A|T]A[A|G|T]",
        "REST" : r"[A|C|G|T][A|C|G|T][A|C|G|T][A|G|T]G[C|G][A|C|T][A|G|T]CTGTCC[A|G|T][A|C|G|T]GGT[C|G]CT",
        "AP2" : r"[C|G][C|G]C[C|G][C|G]C[A|G]GGC[A|C|G|T][A|G][A|C|G|T][A|G][A|C|G|T][A|C|G|T]",
        "AP1" : r"[A|G][C|G]TGAC[A|C|G|T][A|C|G|T][A|C][A|C][A|C|G|T][A|C|T]",
        "MYC" : r"[A|C|G|T][C|G]CACGTGG[A|C|G|T]",
        "USF" : r"[A|C|G|T][A|C|G|T][A|G][A|C|T]CACGTG[A|G|T][C|T][A|C|G|T][A|C|G|T]",
        "SPZ1" : r"[A|G|T][A|C|G|T][A|C|G|T]GG[A|G|T]GG[G|T][A|T][A|C|G|T][A|C|G|T][A|G|T][A|C|G|T][A|C|G]",
        "NRF1" : r"CGC[A|G]TGCGC[A|G]",
        "GATA1" : r"[A|C|G|T][A|C|G|T][A|T]GATAA[A|C|G][A|C|G][A|C|G|T]",
        "YY1" : r"[A|C|G|T][A|C|G][A|C|G|T][A|C|G|T][A|C|G|T]CCAT[A|C|G|T]T[A|T][A|C|G|T][A|C|G|T][A|C|G|T][A|C|G|T][A|C|G|T]",
        "NRF1" : r"CGC[A|G]TGCGC[A|G]",
        "CTCF" : r"CCGCG[A|C|G|T]GG[A|C|G|T]GGCAG",
        "NF1" : r"[C|T][A|G][C|G]C[A|C|G|T][C|G]T[C|G|T][C|G|T][A|C|G|T][A|C|G|T][A|T]TTGGC[A|G][A|C|G|T][C|G|T][C|G][A|C|G|T]GCCA[A|G][A|C|G|T]"}      

    #to indicate if it is forward or reverse

    for i in secuencia:
        count += 1
        if count == 1:
            r = 'forward'
        elif count == 2:
            r = 'reverse'

    #to apply regular expressions for finding patterns and positions
        for key, value in TFBS.items():
            if re.search(value,i) == None:      # to indicate if the TFBS was present
                print('\n\nNo transcription factor binding site ', key,' is found in ',r, ' sequence')
                

            else:
                motif = re.findall(value, i)    # to find all possible TFBSs


                print('\n\nThe motif ', motif, ' for ',key,' is present in the ', r, 'sequence') # to indicate if the TFBS was found


                print('The position in the sequence is: ')


                motifs = re.finditer(value, i)     # to indicate the position of the TFBS 
                for match in motifs:
                    motifs_start = match.start()
                    motifs_end = match.end()
                    print("TFBS from " + str(motifs_start) + " to " + str(motifs_end))

                if r == 'forward':
                    ForwardTFBSs.append(motif)
                elif r == 'reverse':
                    ReverseTFBSs.append(motif)

    print('\nThe forward transcription factor binding sites are: \n\n', ForwardTFBSs,'\n\nand reverse TFBSs are: \n\n', ReverseTFBSs)


    table = [ForwardTFBSs, ReverseTFBSs]
                    
    return table

 
def main():
    import Bio
    from Bio import Motif
    
    from Bio.Alphabet import IUPAC # to import IUPAC to specify alphabet

    fasta= []
    try:
        alignment = openfasta(fasta) # to apply function to open fasta file of interest
        print(alignment)
        n=-1
        for record in alignment:    # to print the records od each sequence in the alignment
            print(record, '\n')     

        code = eval(input("Enter the ID number of the monogamous ancestor of the new world voles: ")) # to ask for the most recent common ancestor between monogamous voles
        print(code)

        
        print('Your monogamous ancestor is: ', alignment[code].id) #to print a confirmation

        forward =""
        for i in sequencing(alignment[code].seq):
            forward += i

        reverse = ""
        for v in reversing(complement(alignment[code].seq)): # to generate a reverse complement
            reverse += v


        dna = [forward, reverse]


        x = transfac(dna)    #to apply the module for transcription factors in the monogamous ancestor

        


        promiscuouscode = eval(input("Enter the ID number of the promiscuous ancestor of the new world voles: ")) # to ask for the most recent common ancestor between promiscuous voles
        print(code)

        
        print('Your promiscuous ancestor is: ', alignment[promiscuouscode].id) #to print a confirmation

        forward =""
        for i in sequencing(alignment[promiscuouscode].seq):
            forward += i

        reverse = ""
        for v in reversing(complement(alignment[promiscuouscode].seq)):   # to generate a reverse complement
            reverse += v


        dna = [forward, reverse]


        y = transfac(dna) #to apply the module for transcription factors in the promiscuous ancestor


        print('\n The monogamous TFBS are: \n',x)    # prints the output
        print('\n The promiscuous TFBS are: \n',y)   # prints the output
            
    except UnboundLocalError:
        print('something went wrong')
main()

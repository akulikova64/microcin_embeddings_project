from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys

# Find all potential open reading frames in both strands.

genome_path = "../../genome_files/genomes_10_known_microcins/ecoli_NCTC11128_N.fna"
output_path = "../../ORF_files/ORFs_ecoli_N/protein_ORFs_ecoli_N.fasta"

genome = SeqIO.parse(genome_path, "fasta")

gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP',
    'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W' }

noncannonical = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']

orf_count = 0
orf = ""
strands = ["(+)", "(-)"]

with open(output_path, 'w') as output_file:
    for rec in genome:
        seq = rec.seq
        name = rec.description.split(" ")
        type = name[5][0:-1]

        for strand in strands:  
            if strand == "(-)": # if we are looking in the (-) strand, our sequence is the reverse complement
                seq = seq.reverse_complement()
            for i in range(0,3):   # 3 possible frames
                found_start = False  # Initally, we have not found a start
                start_loc = 0
                end_loc = 0
                for j in range(i, len(seq), 3):  # moving along the genome sequence every codon
                    codon = seq[j:j+3]

                    if len(codon) == 3:
                        if found_start == False and any(x in noncannonical for x in codon):
                            continue

                        if found_start == True and any(x in noncannonical for x in codon):
                            found_start = False
                            end_loc = j + 3
                            location = "[" + str(start_loc) + "-" + str(end_loc) + "]"
                            length = len(orf)
                            start = 'ATG'
                            stop = 'Noncannonical base in sequencing'
                            frame = i + 1
                            continue

                        aa = gencode[codon]

                        if aa == 'M' and found_start == False:  # If we found a start for the first time after last ORF
                            found_start = True
                            start_loc = j + 1
                            orf += 'M'
                            orf_count += 1

                        elif found_start == True and aa != 'STOP': # Translating/elongating our ORF
                            orf += aa

                        elif found_start == True and aa == 'STOP':  # if we run into stop after start
                            found_start = False
                            end_loc = j + 3
                            location = "[" + str(start_loc) + "-" + str(end_loc) + "]"
                            length = len(orf)
                            start = 'ATG'
                            stop = codon
                            frame = i + 1

                            label = str(">" + name[0] + "_ORF." + str(orf_count) + " " + location + strand + " type:" + type + " length:" + str(length) + " frame:" + str(frame) + " start:" + start + " stop:" + stop)
                            output_file.write(label + "\n")
                            output_file.write(orf + "\n")
                            orf = ""
                            
                    else:
                        break
        
print("Finished extracting ORFs! :-)")

   

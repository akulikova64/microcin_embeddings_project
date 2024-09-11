import os 
import sys
import time
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime

def timestamp():
	return str(datetime.now().time())

def get_composition_vector(seq):
    # order of vector: 
    aa_dict = {'A':0,'R':0,'N':0,'D':0,'C':0,'Q':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}

    for aa in seq:
        aa_dict[aa] += 1
        
    aa_list = aa_dict.values()
    aa_sum = sum(aa_list)
    
    new_aa_list = []
    for aa_count in aa_list:
        new_aa_list.append(aa_count/aa_sum)
    
    return new_aa_list


#--------------------------main--------------------------------------------------

# list of genomes here:
genome_path = '../../ORF_files/ORFs_10_known_microcins/'
file_list = os.listdir(genome_path)

output_path = '../../ORF_files/ORFs_10_known_microcins_composition/'
isExist = os.path.exists(output_path)
if not isExist:
    os.makedirs(output_path)

for file in file_list:

    input_path = '../../ORF_files/ORFs_10_known_microcins/' + file + ''
    orf_list = list(SeqIO.parse(input_path, "fasta"))

    with open(output_path + file, 'w') as output_file:
        for rec in orf_list:
            
            # extracting ORF information from the label/name
            orf_sequence = str(rec.seq)
            orf_description = str(rec.description)

            composition_vector = get_composition_vector(orf_sequence)

            output_file.write(">" + orf_description + "\n")
            output_file.write(str(composition_vector) + "\n")
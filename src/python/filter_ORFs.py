from Bio import SeqIO
from Bio.Seq import Seq
import sys
import re
import string


# filter ORFs by size and find smaller frames

#Note: Alternative start codons can occur within the larger open reading frames. 
# Since we are filtering by size, we might remove all ORFs where the larger frame is too big. 
# The smaller alternative sequences might be under the size limit and we don't want to miss them.
# This script deals with this issue.


input_path = "../../ORF_files/ORFs_ecoli_N/protein_ORFs_ecoli_N.fasta"
output_path = "../../ORF_files/ORFs_ecoli_N/protein_ORFs_ecoli_N_filtered.fasta"

orf_list = list(SeqIO.parse(input_path, "fasta"))

with open(output_path, 'w') as output_file:
    for rec in orf_list:
        
        # extracting ORF information from the label/name
        orf_sequence = str(rec.seq)
        orf_description = str(rec.description)
        split_orf = re.split(r'[_.() ]', orf_description)
        location = split_orf[5] # use to be 4
        split_location = re.split(r'[\[\-\]]', location)
        start_loc = int(split_location[1])
        end_loc = int(split_location[2])

        if len(orf_sequence) <= 150 and len(orf_sequence) > 30: #if the ORF is already within size limit, we include it.
            output_file.write(">" + orf_description + "\n")
            output_file.write(orf_sequence + "\n")
            continue # if the ORF is already within limits, we do not have to search further

        if len(orf_sequence) > 30: # excluding ORFs that are too small right away
            sub_ORF_count = 0 # how many sub-ORFs I find
            for position, aa in enumerate(orf_sequence):
                if position != 0 and aa == 'M': # if "M" is found within an ORF (not the first residue)
                    sub_ORF_count += 1
                    new_start = start_loc + position + 1 # recalculate the new start (end will be the same)
                    new_location = "[" + str(new_start) + "-" + str(end_loc) + "]" 

                    # next five lines create a new fasta label/name for the new sub-ORF:
                    split_orf = re.split(r'[\[\]]', orf_description)
                    split_orf[0] = split_orf[0].replace(' ', '-' + str(sub_ORF_count)) 
                    split_orf[1] = new_location 
                    new_label = ' '.join(split_orf)
                    new_label = new_label.replace('] (', '](')
                    
                    # getting the sub-ORF sequence:
                    new_orf = "M" + orf_sequence[position + 1:]
                    new_length = len(new_orf)

                    # if sub-ORF is within size limit, then add it to the new fasta:
                    if len(new_orf) <= 150 and len(new_orf) > 30:
                        output_file.write(">" + new_label + "\n")
                        output_file.write(new_orf + "\n")
                        break # onece we add an ORF, we do not need to keep looking for subORFs

print("Finished filtering ORFs! :-)")
print("Next, use the following command to generate esm1-b embeddings:")
print("python3 extract.py esm1b_t33_650M_UR50S ../../ORF_files/ORFs_HUW04/protein_ORFs_HUW04_filtered.fasta ../../embeddings/genome_embeddings_HUW04/ --include mean per_tok")
          
#python3 extract.py esm1b_t33_650M_UR50S ../../ORF_files/ORFs_ecoli_L/protein_ORFs_ecoli_L_filtered.fasta ../../embeddings/genome_embeddings_ecoli_L/ --include mean per_tok
           

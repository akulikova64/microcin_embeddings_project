import gemmi
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import csv
import re

# align sequences


def get_alignment(seq1, seq2):
    """
    Returns an alignment and its score
    """

    alignment = gemmi.align_string_sequences(list(seq1), list(seq2), [])

    #print(alignment.formatted(seq1, seq2), end='') 
    score = alignment.score

    return alignment.formatted(seq1, seq2), score

output_path = "../../alignments/ORF26669_alignments_vibrio_cholerae_HUW04.txt"

orf_seq = "MAYSLLNLGADTRKRALAGMQESAQREEQRNQTNQSLKDAQRTKRLSSVTTGAGMGMMAGMQAGSVGGPMGAAMGAAAGLILGELF"
input_path_microcins =  "../../microcin_files/45_microcins.fasta"
microcin_list = list(SeqIO.parse(input_path_microcins, "fasta"))

with open(output_path, 'w') as output_file:
    i = 0
    for rec in microcin_list:
        i += 1
        microcin_sequence = str(rec.seq)
        microcin_description = str(rec.description)
        microcin_name = microcin_description.split("|")[0][0:-3]

        alignment, score = get_alignment(orf_seq, microcin_sequence)
            
        output_file.write(str(i+1) + ")\n" + "score: " + str(score) + ", microcin: " + str(microcin_name) + "\nseq_1: best_hit \nseq_2: microcin" + "\n\n" + alignment + "\n\n")
        
print("Complete! :-)")
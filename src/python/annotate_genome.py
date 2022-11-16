import gemmi
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import csv
import re

# this script annotates genome for known microcins. 

def get_alignment_score(orf_seq, microcin_seq):
    """
    Returns alignment score
    """

    alignment = gemmi.align_string_sequences(list(orf_seq), list(microcin_seq), [])

    return alignment

input_path_microcins =  "../../microcin_files/Microcins_Known.fasta"
input_path_orfs = "../../ORF_files/ORFs_ecoli_H47/protein_ORFs_ecoli_H47_filtered.fasta"
output_path = "../../genome_annotation_csvs/ecoli_H47_annotation.csv" # output path contains the orf ID and whether or not it's a microcin. 

#input_path_orfs = "../../ORF_files/ORFs_ecoli_nissle/protein_ecoli_nissle_orfs.fasta"
#output_path = "../../genome_annotation_csvs/ecoli_nissle_annotation.csv" # output path contains the orf ID and whether or not it's a microcin. 


microcin_list = list(SeqIO.parse(input_path_microcins, "fasta"))
orf_list = list(SeqIO.parse(input_path_orfs, "fasta"))

with open(output_path, "w") as CSV_file:
    writer = csv.writer(CSV_file)
    writer.writerow(['orf_number', 'orf_strand', 'orf_location', 'is_microcin', 'closest_microcin', 'alignment_score'])

    # found microcins genome-wide
    found_microcins = {}
    count_found = 0


    for rec in orf_list:
        higheset_score = -1000000
        microcin_winner = ""
        microcin_found = "" 

        # orf info
        orf_sequence = str(rec.seq)
        orf_description = str(rec.description)
        split_orf = re.split(r'[_.() ]', orf_description)
        # print(orf_description)
        # print(split_orf)
        # print()
        orf_number = split_orf[4]
        orf_strand = "(" + split_orf[6] + ")"
        orf_location = split_orf[5]

        for rec in microcin_list:
            microcin_sequence = str(rec.seq)
            microcin_description = str(rec.description)
            microcin_name = microcin_description.split("|")[0]#[0:-3]

            # get alignment of orf and microcin, update winning microcin 
            alignment = get_alignment_score(orf_sequence, microcin_sequence)
            alignment_score = alignment.score

            if alignment_score > higheset_score:
                higheset_score = alignment_score
                microcin_winner = microcin_name

            # see if orf matches microcin sequence
            if orf_sequence == microcin_sequence:
                found_microcins[orf_number] = microcin_name
                count_found += 1
                microcin_found = microcin_name
                break  
            elif alignment_score > 30:
                found_microcins[orf_number] = microcin_name
                count_found += 1
                microcin_found = microcin_name 
                break  
            else:
                microcin_found = "non-microcin"


        writer.writerow([orf_number, orf_strand, orf_location, microcin_found, microcin_winner, higheset_score])

print("Number of microcins found:", count_found)
print("Found microcins:", found_microcins)
print()
print("Job complete! :-)")

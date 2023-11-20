from Bio import SeqIO
import numpy as np
import torch
import gemmi
import math
import csv
import sys
import os
import re

# This scrip will calc the distance between a new microcin and the 10 knwon microcins as well as the percent 
# sequence similarity to the 10 known microcins.
def get_distance(emb_1, emb_2):
    sum = 0
    for x1, x2 in zip(emb_1, emb_2):
        sum += pow((x2 - x1), 2)
    
    distance = round(math.sqrt(sum), 4)
    distance = float(distance)
    
    return distance

def get_percent(seq1, seq2):
  count_similar = 0
  for i, j in zip(seq1, seq2):
    assert len(seq1) == len(seq2)
    if i == j: # gaps also considered a mismatch
      count_similar += 1

  perc_sim = float(count_similar/len(seq1)) * 100
  return perc_sim

def get_alignment(seq1, seq2):
    """
    Returns an alignment and its score
    """

    alignment = gemmi.align_string_sequences(list(seq1), list(seq2), [])

    return alignment.formatted(seq1, seq2)

input_path_microcins =  "../../microcin_files/Microcins_Known_emb_esm1b/"
input_path_microcin_seqs = "../../microcin_files/Microcins_Known.fasta"

input_path_new_microcin_seq = "../../microcin_files/new_microcin.fasta" 
input_path_new_microcin_emb = "../../microcin_files/new_microcin_emb.pt" 
output_path = "../../microcin_files/new_microcin_analysis_homologs.csv" 

microcin_emb_list = os.listdir(input_path_microcins)
microcin_seq_list = list(SeqIO.parse(input_path_microcin_seqs, "fasta"))

microcin_seq_dict = {}
for i in range(0, len(microcin_seq_list)):
            name = microcin_seq_list[i].name
            name = name.split("_")[0]
            seq = microcin_seq_list[i].seq
            microcin_seq_dict[name] = seq


with open(output_path, "w") as CSV_file:
    writer = csv.writer(CSV_file)
    writer.writerow(['microcin', 'distance', 'sequence_divergence'])

    new_microcin_seq_fasta = list(SeqIO.parse(input_path_new_microcin_seq, "fasta"))
    new_microcin_seq = new_microcin_seq_fasta[0].seq

    for microcin in microcin_emb_list:

        microcin_emb_name = microcin.split("_")[0]
        print("microcin:", microcin_emb_name)

        seq_known = microcin_seq_dict[microcin_emb_name]

        # get sequence div:
        alignment = get_alignment(str(seq_known), str(new_microcin_seq))
        print("alignment:")
        print(alignment)

        print()
        seq1 = alignment.split("\n")[0]
        seq2 = alignment.split("\n")[2]
        print(seq1)
        print(seq2)

        percent_sim = get_percent(str(seq1), str(seq2))
        print("percent sim:", percent_sim)
        percent_div = 100 - percent_sim
        print("percent div:", percent_div)
        print()
        print("-------------------------------------------------------------------")

        # get distance:
        micr_embedding = torch.load(input_path_microcins + microcin)
        mean_emb_microcin = micr_embedding['mean_representations'][33]

        new_microcin_embedding = torch.load(input_path_new_microcin_emb)
        new_microcin_mean_emb = new_microcin_embedding['mean_representations'][33]

        distance = get_distance(mean_emb_microcin, new_microcin_mean_emb)
                
        writer.writerow([microcin_emb_name, distance, percent_div])
           
print("Job complete! :-)")
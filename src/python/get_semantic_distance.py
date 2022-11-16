import numpy as np
import torch
import math
import csv
import sys
import os
import re

# this script computes the distance between two embeddings and saves all data to csv file
# the embedding used here are mean representations (average across residues)

def get_distance(emb_1, emb_2):
    sum = 0
    for x1, x2 in zip(emb_1, emb_2):
        sum += pow((x2 - x1), 2)
    
    distance = round(math.sqrt(sum), 4)
    distance = float(distance)
    
    return distance

input_path_microcins =  "../../microcin_files/Microcins_Known_emb_esm1b/"
input_path_orfs = "../../embeddings/genome_embeddings_ecoli_N/" 
genome = "ecoli_N"

microcin_list = os.listdir(input_path_microcins)

for microcin in microcin_list:
    microcin_name = microcin.split("|")[0]
    output_path = "../../distance_csvs/" + genome + "/distances_" + microcin_name + "_vs_" + genome + ".csv"

    with open(output_path, "w") as CSV_file:
        writer = csv.writer(CSV_file)
        writer.writerow(['genome', 'orf_number', 'orf_strand', 'microcin', 'orf_seq', 'distance', 'orf_location'])

        micr_embedding = torch.load(input_path_microcins + microcin)
        mean_emb_microcin = micr_embedding['mean_representations'][33]

        orf_list = os.listdir(input_path_orfs)
        for orf in orf_list:
            try:
                orf_embedding = torch.load(input_path_orfs + orf)
                mean_emb_orf = orf_embedding['mean_representations'][33]
                distance = get_distance(mean_emb_microcin, mean_emb_orf)
                split_orf = re.split(r'[_.() ]', orf)
                orf_number = split_orf[3]
                orf_strand = "(" + split_orf[5] + ")"
                orf_location = split_orf[4]
                orf_seq = "NA"
                writer.writerow([genome, orf_number, orf_strand, microcin, orf_seq, distance, orf_location])
            except: 
                continue

print("Job complete! :-)")
            

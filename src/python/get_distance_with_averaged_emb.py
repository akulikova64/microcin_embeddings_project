import numpy as np
import torch
import math
import csv
import sys
import os
import re

# get distances between average microcin embedding (10 microcins) and all the ORFs. 

def get_distance(emb_1, emb_2):
    sum = 0
    for x1, x2 in zip(emb_1, emb_2):
        sum += pow((x2 - x1), 2)
    
    distance = round(math.sqrt(sum), 4)
    distance = float(distance)
    
    return distance

genome = sys.argv[1]
dataset = sys.argv[2]

#------------paths-------------------------
input_path_microcin_avg =  "../../microcin_files/average_emb_10_esm1b/average_10_microcin_emb.pt"

input_path_orfs = "../../embeddings/" + dataset + "/genome_embeddings_" + genome + "/" 
isExist = os.path.exists(input_path_orfs)
if not isExist:
   os.makedirs(input_path_orfs)

output_path = "../../distance_csvs/" + dataset + "/" + genome + "/distances_average_vs_" + genome + ".csv"
output_folder = "../../distance_csvs/" + dataset + "/" + genome + "/"
isExist = os.path.exists(output_folder)
if not isExist:
    os.makedirs(output_folder)
#-------------------------------------------


with open(output_path, "w") as CSV_file:
    writer = csv.writer(CSV_file)
    writer.writerow(['genome', 'orf_number', 'orf_strand', 'microcin', 'orf_seq', 'distance', 'orf_location'])

    mean_emb_microcin = torch.load(input_path_microcin_avg)

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
            writer.writerow([genome, orf_number, orf_strand, 'average', orf_seq, distance, orf_location])
        except: 
            continue

print("Job complete! :-)")

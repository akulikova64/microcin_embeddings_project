import numpy as np
import torch
import math
import csv
import sys
import os
import re


# get distance between the 10 microcin embeddings.
# the embedding used here are mean representations (average across residues)

def get_distance(emb_1, emb_2):
    sum = 0
    for x1, x2 in zip(emb_1, emb_2):
        sum += pow((x2 - x1), 2)
    
    distance = round(math.sqrt(sum), 4)
    distance = float(distance)
    
    return distance

input_path_microcins =  "./Microcins_Known_emb_esm1b/"
output_path = "./distance_csvs/distance_bw_microcins.csv"
genome = "ecoli_s88"

microcin_list_1 = os.listdir(input_path_microcins)
microcin_list_2 = os.listdir(input_path_microcins)


with open(output_path, "w") as CSV_file:
    writer = csv.writer(CSV_file)
    writer.writerow(['microcin1', 'microcin2', 'distance'])

    for microcin1 in microcin_list_1:
        microcin1_name = microcin1.split("|")[0]

        micr_embedding1 = torch.load(input_path_microcins + microcin1)
        mean_emb_microcin1 = micr_embedding1['mean_representations'][33]
    
        for microcin2 in microcin_list_2:
            microcin2_name = microcin2.split("|")[0]

            micr_embedding2 = torch.load(input_path_microcins + microcin2)
            mean_emb_microcin2 = micr_embedding2['mean_representations'][33]

            # get distance
            distance = get_distance(mean_emb_microcin1, mean_emb_microcin2)

            writer.writerow([microcin1_name, microcin2_name, distance])

print("Job complete! :-)")
        

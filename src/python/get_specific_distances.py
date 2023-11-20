import numpy as np
import torch
import math
import csv
import sys
import os
import re

def get_distance(emb_1, emb_2):
    sum = 0
    for x1, x2 in zip(emb_1, emb_2):
        sum += pow((x2 - x1), 2)
    
    distance = round(math.sqrt(sum), 4)
    distance = float(distance)
    
    return distance

input_path_embeddings =  "../../embeddings/genome_embeddings_ecoli_L/"
output_path = "../../analysis/distance_bw_nonmicrocins_L_194.csv"
pairs_csv = "../../analysis/194_sampled_ORFs_by_perc_sim.csv"
genome = "ecoli_L"


with open(output_path, "w") as CSV_file:
    writer = csv.writer(CSV_file)
    #writer.writerow(['microcin1', 'microcin2', 'distance'])
    writer.writerow(['seq1', 'seq2', 'distance'])

    with open(pairs_csv, "r") as CSV_file:
        csv_reader = csv.reader(CSV_file, delimiter=',')
        next(csv_reader)

        for row in csv_reader:
            seq1 = row[0]
            seq2 = row[1]

            orf_list_1 = os.listdir(input_path_embeddings)
            for orf1 in orf_list_1:
                orf1_name = orf1.split(" ")[0]
                orf1_name = orf1_name.split(".")[2]

                if seq1 == orf1_name:

                    embedding1 = torch.load(input_path_embeddings + orf1)
                    mean_emb_1 = embedding1['mean_representations'][33]
            
                    orf_list_2 = os.listdir(input_path_embeddings)
                    for orf2 in orf_list_2:
                        orf2_name = orf2.split(" ")[0]
                        orf2_name = orf2_name.split(".")[2]

                        if seq2 == orf2_name:

                            embedding2 = torch.load(input_path_embeddings + orf2)
                            mean_emb_2 = embedding2['mean_representations'][33]

                            # get distance
                            distance = get_distance(mean_emb_1, mean_emb_2)

                            writer.writerow([orf1_name, orf2_name, distance])

print("Job complete! :-)")
        

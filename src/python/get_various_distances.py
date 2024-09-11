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

def get_cosine_distance(emb_1, emb_2):
    """
    Compute the cosine similarity between two vectors.

    Parameters:
    emb_1 (numpy.ndarray): The first vector.
    emb_2 (numpy.ndarray): The second vector.

    Returns:
    float: The cosine similarity between the two vectors.
    """

    dot_product = np.dot(emb_1, emb_2)
    norm_emb_1 = np.linalg.norm(emb_1)
    norm_emb_2 = np.linalg.norm(emb_2)
    cosine_dist = dot_product / (norm_emb_1 * norm_emb_2)

    return 1 - cosine_dist

def get_manhattan_distance(vec1, vec2):
    #print("entered manhattan function")
    
    #manhattan_distance = np.sum(np.abs(vec1 - vec2))
    manhattan_distance = sum(abs(val1-val2) for val1, val2 in zip(vec1, vec2))
    manhattan_distance = manhattan_distance.item()
    return manhattan_distance


# set distance type to "cosine" or "manhattan"
distance_type = "manhattan"

#-------------------------------------------
input_path_microcins =  "../../microcin_files/Microcins_Known_emb_esm1b/" # microcin embeddings are here:

genome_list =  ["ecoli_L", "ecoli_N"] #, "ecoli_pcolv-k30_V"]#, "ecoli_PDI", "ecoli_S", "ecoli_H47", "k_pneumoniae_cluster"] # list of genome names

for genome in genome_list:
    print("Working on:", genome)
    input_path_orfs = '../../embeddings/genome_embeddings_' + genome + '/'  # orf embeddings are here: 
    microcin_list = os.listdir(input_path_microcins)

    for microcin in microcin_list:
        microcin_name = microcin.split("|")[0]
        output_path = "../../distance_csvs/10_microcin_files_" + distance_type + "/" + genome + "/distances_" + distance_type + "_" + microcin_name + "_vs_" + genome + ".csv"
        output_folder = "../../distance_csvs/10_microcin_files_" + distance_type + "/" + genome + "/"
        isExist = os.path.exists(output_folder)
        if not isExist:
            os.makedirs(output_folder)

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

                    if distance_type == "cosine":
                        distance = get_cosine_distance(mean_emb_microcin, mean_emb_orf)
                    elif distance_type == "manhattan":
                        distance = get_manhattan_distance(mean_emb_microcin, mean_emb_orf)

                    split_orf = re.split(r'[_.() ]', orf)
                    #if genome == "ecoli_L" | genome == "ecoli_N":
                    #print(split_orf)
                    orf_number = split_orf[4]
                    #print(orf_number)
                    orf_strand = "(" + split_orf[6] + ")"
                    #print(orf_strand)
                    orf_location = split_orf[5]
                    #print(orf_location)
                    #sys.exit()
                    orf_seq = "NA"

                    #else:
                    '''
                    orf_number = split_orf[3]
                    orf_strand = "(" + split_orf[5] + ")"
                    orf_location = split_orf[4]
                    orf_seq = "NA"
                    '''

                    writer.writerow([genome, orf_number, orf_strand, microcin, orf_seq, distance, orf_location])
                except: 
                    print("error")
                    #sys.exit()
                    continue

print("Job complete! :-)")
            

import os 
import re
import math
import sys
import csv
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

def get_distance(emb_1, emb_2):

    sum = 0
    for x1, x2 in zip(emb_1, emb_2):
        sum += pow((x2 - x1), 2)
    
    distance = round(math.sqrt(sum), 4)
    distance = float(distance)
    
    return distance

import numpy as np

def get_cosine_similarity(emb_1, emb_2):
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
    cosine_sim = dot_product / (norm_emb_1 * norm_emb_2)

    return 1 - cosine_sim

def get_manhattan_distance(vec1, vec2):
    
    manhattan_distance = np.sum(np.abs(vec1 - vec2))

    return manhattan_distance



#-------------------------------------------
# get microcin list here:

input_path_microcins = '../../microcin_files/Microcins_Known.fasta'
microcin_list = list(SeqIO.parse(input_path_microcins, "fasta"))

microcin_names = []
microcin_compositions = []

for microcin in microcin_list:
    microcin_compositions.append(get_composition_vector(microcin.seq))
    microcin_names.append(microcin.description.split("|")[0])

#--------------------------------------------
# taking each of the 10 microcins consecutively:
for name, microcin in zip(microcin_names, microcin_compositions):

    input_path_orfs = '../../ORF_files/ORFs_10_known_microcins_composition/'
    orf_file_list = os.listdir(input_path_orfs)

    for orf_file in orf_file_list:
        orf_file_name = orf_file.split("_")[3]
        output_folder = '../../distance_csvs/10_microcin_files_composition/' + orf_file_name + "/"
        isExist = os.path.exists(output_folder)
        if not isExist:
            os.makedirs(output_folder)
        output_path = output_folder + 'distances_' + name + '_vs_' + orf_file_name + '.csv'
       
        with open(output_path, "w") as CSV_file:
            writer = csv.writer(CSV_file)
            writer.writerow(['genome', 'orf_number', 'orf_strand', 'microcin', 'orf_seq', 'distance', 'orf_location'])

            orf_list = list(SeqIO.parse('../../ORF_files/ORFs_10_known_microcins_composition/' + orf_file, "fasta"))
            for orf in orf_list:
                #orf_genome = orf.description.split(" ")[0]
                orf_seq = str(orf.seq).replace("[", "").replace("]", "") 
                orf_seq = orf_seq.split(",")
                
                
                microcin = [float(i) for i in microcin]
                
                try:
                    orf_seq = [float(i) for i in orf_seq]
                except:
                    print(orf_seq)
                    print(orf_file_name)
        
                genome = orf_file_name
            
                split_orf = re.split(r'[.() ]', orf.description)
                orf_number = split_orf[2]
                orf_strand = "(" + split_orf[4] + ")"
                orf_location = split_orf[3]
                distance = get_distance(microcin, orf_seq)
                orf_seq = "NA"
                microcin_name = name

                writer.writerow([genome, orf_number, orf_strand, microcin_name, orf_seq, distance, orf_location])

print("Finished calculating distances!")
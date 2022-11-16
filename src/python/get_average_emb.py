import numpy as np
import torch
import sys
import os

# this program will find the average embedding from the 10 known microcin embeddings

input_path = "../../microcin_files/Microcins_Known_emb_esm1b/"
output_path =  "../../microcin_files/average_emb_10_esm1b/average_10_microcin_emb.pt"


microcin_list = os.listdir(input_path)
emb_list = []

for microcin in microcin_list:
    microcin_name = microcin.split("|")[0]
    micr_embedding = torch.load(input_path + microcin)
    mean_emb_microcin = micr_embedding['mean_representations'][33]
    vector = mean_emb_microcin.tolist()
    emb_list.append(vector)
    
avg_vector = np.mean(emb_list, axis = 0)
torch.save(avg_vector, output_path)
    




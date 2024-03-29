import os
import re
import csv
import sys
import torch

# get PCA-ready CSV file

output_path = "../../for_PCA/ecoli_L_PCA_data.csv"
microcin_path = "../../microcin_files/Microcins_Known_emb_esm1b/"
orf_path = "../../embeddings/ecoli_L/"

microcin_embeddings = os.listdir(microcin_path)
orf_embeddings = os.listdir(orf_path)

with open(output_path, "w") as CSV_file:
    writer = csv.writer(CSV_file)
    row = ['group', 'name']
    for i in range(1, 1280 + 1):
        row.append(str(i))

    writer.writerow(row)

    for microcin in microcin_embeddings:
        microcin_name = microcin.split("|")[0][0:-3]
        group = "microcin"
        embedding = torch.load(microcin_path + microcin)
        avg_embedding = embedding['mean_representations'][33]
        list = avg_embedding.tolist()
        writer.writerow([group, microcin_name] + list)

    for orf in orf_embeddings:
        split_orf = re.split(r'[_.() ]', orf)
        #print(split_orf)
        #sys.exit()
        #try:
        orf_number = split_orf[4]
        #print(orf_number)
        #sys.exit()
        group = "orf"
        embedding = torch.load(orf_path + orf)
        avg_embedding = embedding['mean_representations'][33]
        list = avg_embedding.tolist()
        writer.writerow([group, orf_number] + list)
        #print(group, orf_number)
        #sys.exit()

        #except:
           # print(split_orf)


print("Finished making CSV! :-)")
#Note: last layer in esm1-b is layer 33.



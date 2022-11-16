from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import sys
import csv


# this script will parse through all the blast and annotation files and extract the blast rankings for the known microcins.
# rankings are based on distance between the query embeddings (microcin embeddings) and orf embeddings with the known microcin. 


annotation_csvs = "../../genome_annotation_csvs/10_microcin_files/"
blastp_files = "../../BLAST_results/10_known_microcins/"
output_path = "../../analysis/rankings/blastp_rankings/"

# create an Empty DataFrame object to store all data (from all genomes)
df_list = []
annot_genome_list = os.listdir(annotation_csvs)

for main_folder_path, genome_folders, blast_files in os.walk(blastp_files):
    for filename in blast_files:
        #print(filename)
        blast_df = pd.read_csv(main_folder_path + "/" + filename, sep='\t', usecols=[0, 1, 2, 10], names=['query', 'subject', 'percent_sim', 'e-value'])
        path_list = main_folder_path.split("/")
        blast_genome = path_list[4][13:-9]

        for genome_file in annot_genome_list:
            annot_genome = genome_file[0:-15]

            if annot_genome == blast_genome:
                annotation_df = pd.read_csv(annotation_csvs + genome_file)
                #print(annotation_df)
                blast_df['orf_number'] = blast_df['subject'].str.split(".").str[2]
                blast_df['query'] = blast_df['query'].str.split("|").str[0].str[0:-3]
                #blast_df['orf_number'] = blast_df['orf_number'].astype('int') 
                annotation_df = annotation_df[['orf_number', 'is_microcin']]
                
                joined_df = pd.merge(left = annotation_df, right = blast_df, left_on = ['orf_number'], right_on=['orf_number'])
                

print(df_list[0:2])
'''             
complete_df = pd.concat(df_list)
complete_df.to_csv(output_path + "10_known_microcin_blast_finds.csv", index=False) 
print("Complete! :-)") 
'''

      
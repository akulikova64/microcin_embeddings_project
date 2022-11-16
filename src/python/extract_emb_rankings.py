from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import sys
import csv


# this script will parse through all the distance and annotation files and extract the rankings for the known microcins.
# rankings are based on distance between the query embeddings (microcin embeddings) and orf embeddings with the known microcin. 


# steps:
# 1) combine annotation with correponding distance csv. 

# resulting filtes will have the following header:
#  genome, query, found_microcin, ORF, rank, distance


annotation_csvs = "../../genome_annotation_csvs/10_microcin_files/"
distance_csvs = "../../distance_csvs/10_microcin_files/"
output_path = "../../analysis/rankings/10_known_microcin_rankings/"

# create an Empty DataFrame object to store all data (from all genomes)
df_list = []
annot_genome_list = os.listdir(annotation_csvs)


for main_folder_path, genome_folders, distance_files in os.walk(distance_csvs):
    for filename in distance_files:
        #print(filename)
        distance_df = pd.read_csv(main_folder_path + "/" + filename)
        
        for genome_file in annot_genome_list:
            annot_genome = genome_file[0:-15]
            distance_df_genome = distance_df.loc[0].at['genome']
            query_name = distance_df.loc[0].at['microcin']
            query = query_name.split("_")[0]

            if annot_genome == distance_df_genome:
                annotation_df = pd.read_csv(annotation_csvs + genome_file)

                joined_df = pd.merge(left = annotation_df, right = distance_df, left_on = ['orf_number','orf_location', 'orf_strand'], right_on=['orf_number','orf_location', 'orf_strand'])
                joined_df = joined_df.sort_values(by=['distance'])
                joined_df.insert(loc = 0, column = 'rank', value = list(range(1, len(joined_df)+1)))
                joined_df =  joined_df[joined_df['is_microcin'] != 'non-microcin']
                joined_df = joined_df.rename(columns={"closest_microcin": "found_microcin"})
                joined_df = joined_df[["genome", "found_microcin", "orf_number", "rank", "distance"]]
                joined_df['query'] = query
                joined_df = joined_df[["genome", "query", "found_microcin", "orf_number", "rank", "distance"]]
                df_list.append(joined_df)


complete_df = pd.concat(df_list)
complete_df.to_csv(output_path + "10_known_microcin_rankings.csv", index=False) 
print("Complete! :-)") 




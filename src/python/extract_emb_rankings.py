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

# resulting files will have the following header:
#  genome, query, found_microcin, ORF, rank, distance


annotation_csvs = "../../genome_annotation_csvs/10_microcin_files/"
distance_csvs = "../../distance_csvs/10_microcin_files_cosine/"
output_path = "../../analysis/rankings/10_known_microcin_rankings_cosine/"
isExist = os.path.exists(output_path)
if not isExist:
    os.makedirs(output_path)


# create an Empty DataFrame object to store all data (from all genomes)
df_list = []
annot_genome_list = os.listdir(annotation_csvs)

for main_folder_path, genome_folders, distance_files in os.walk(distance_csvs):
    for filename in distance_files:
        if "DS_Store" in filename:
            continue
        #print("left", filename)
        try:
            distance_df = pd.read_csv(main_folder_path + "/" + filename)
            
            for genome_file in annot_genome_list:
                #print("right", genome_file)
                annot_genome = genome_file[0:-14]
                distance_df_genome = distance_df.loc[0].at['genome']
                distance_df_genome = distance_df_genome + "_"

                #print(annot_genome)
                #print(distance_df_genome)
                #print(distance_df_genome in annot_genome)
                #sys.exit()

                query_name = distance_df.loc[0].at['microcin']
                query = query_name.split("_")[0]

                if  distance_df_genome in annot_genome:
                    print("current query:", query)

                    #print("distance_df_genome:", distance_df_genome, "annot_genome:", annot_genome)
                    annotation_df = pd.read_csv(annotation_csvs + genome_file)
                    annotation_df["orf_number"] = annotation_df["orf_number"].astype(str)
                    distance_df["orf_number"] = distance_df["orf_number"].astype(str)
                    #print("annotation_df:", annotation_df["orf_number"].dtypes)
                    #print("distance_df:", distance_df["orf_number"].dtypes)
                    #print("annotation_df:", annotation_df)
                    #print("distance_df:", distance_df)
                    #sys.exit()
                    joined_df = pd.merge(left = annotation_df, right = distance_df, left_on = ['orf_number','orf_location', 'orf_strand'], right_on=['orf_number','orf_location', 'orf_strand'])
                    joined_df = joined_df.sort_values(by = ['distance'], ascending = False) # sorting by similarity (descending)
                    joined_df.insert(loc = 0, column = 'rank', value = list(range(1, len(joined_df)+1))) # adding a "rank" column
                    joined_df =  joined_df[joined_df['is_microcin'] != 'non-microcin'] # removing non microcin columns
                    joined_df = joined_df.rename(columns={"closest_microcin": "found_microcin"}) # changing column name to "found microcin"
                    joined_df = joined_df[["genome", "found_microcin", "orf_number", "rank", "distance"]] # keeping only 5 columns
                    joined_df['query'] = query # adding a column- the query microcin
                    joined_df = joined_df[["genome", "query", "found_microcin", "orf_number", "rank", "distance"]] # reordering columns
                    print(joined_df)
                    print()
                    df_list.append(joined_df)

                else:
                    continue
                    #print("first:", annot_genome)
                    #print("second:", distance_df_genome)
                    #sys.exit()
        except Exception as e: 
            print(e)
            print("Error")
            print(main_folder_path + "/" + filename)
            print()
            sys.exit()



complete_df = pd.concat(df_list)
complete_df.to_csv(output_path + "10_known_microcin_rankings_cosine.csv", index=False) 
print("Complete! :-)") 




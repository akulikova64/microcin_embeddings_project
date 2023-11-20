import os 
import sys
import time
import numpy as np
from datetime import datetime

def timestamp():
	return str(datetime.now().time())

#--------------------------main--------------------------------------------------

# list of genomes here:
dataset = "klebsiella_47"
genome_path = '../../genome_files/' + dataset + '/'
genome_list = os.listdir(genome_path)

for genome in genome_list:
    genome = genome.split(".")[0]

    orf_path = '../../ORF_files/' + dataset + '/ORFs_' + genome + '/protein_ORFs_' + genome + '_filtered.fasta'
    orf_folder = '../../ORF_files/' + dataset + '/ORFs_' + genome + '/'
    isExist = os.path.exists(orf_folder)
    if not isExist:
        os.makedirs(orf_folder)

    embedding_path = '../../embeddings/' + dataset + '/genome_embeddings_' + genome + '/'
    isExist = os.path.exists(embedding_path)
    if not isExist:
        os.makedirs(embedding_path)

    list_scripts = ['find_ORFs.py ' + genome + ' ' + dataset,
                    'filter_ORFs.py ' + genome + ' ' + dataset,
                    'annotate_genome.py ' + genome + ' ' + dataset,
                    'extract.py esm1b_t33_650M_UR50S ' + orf_path + ' ' + embedding_path + ' --include mean per_tok',
                    'get_semantic_distance.py ' + genome + ' ' + dataset,
                    'get_distance_with_averaged_emb.py ' + genome + ' ' + dataset]

    for script in list_scripts: 
        print()
        print("Running the following script:", script)
        print("Genome:", genome, "Time:", timestamp(), "\n")
        print()
        os.system('python ' + script)




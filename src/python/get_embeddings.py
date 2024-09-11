import os 
import sys
import time
from datetime import datetime

def timestamp():
	return str(datetime.now().time())

#--------------------------main--------------------------------------------------
genome_list =  ["ecoli_H47", "ecoli_L", "ecoli_N", "ecoli_pcolv-k30_V", "ecoli_PDI", "ecoli_S", "k_pneumoniae_cluster"] # list of genome names

for genome in genome_list:

    orf_path = '../../ORF_files/ORFs_10_known_microcins/protein_ORFs_' + genome + '_filtered.fasta'
    orf_folder = '../../ORF_files/ORFs_10_known_microcins/'

    embedding_path = '../../embeddings/genome_embeddings_' + genome + '/'
    isExist = os.path.exists(embedding_path)
    if not isExist:
        os.makedirs(embedding_path)

    os.system('python ' + 'extract.py esm1b_t33_650M_UR50S ' + orf_path + ' ' + embedding_path + ' --include mean per_tok')

print("Complete!")
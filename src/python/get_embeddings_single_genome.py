import os

# get embeddings from a single genome:

genome = 'L'
orf_path = '../../ORF_files/ORFs_ecoli_L/protein_ORFs_ecoli_' + genome + '_filtered.fasta'
embedding_path = '../../embeddings/ecoli_L/'

script = 'extract.py esm1b_t33_650M_UR50S ' + orf_path + ' ' + embedding_path + ' --include mean per_tok'

os.system('python ' + script)
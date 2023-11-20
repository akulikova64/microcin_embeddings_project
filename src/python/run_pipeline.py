import os 
import sys
import time
import numpy as np
from datetime import datetime
import multiprocessing as mp

def timestamp():
	return str(datetime.now().time())

def run_scripts(genomes):
	for genome in genomes:
		orf_path = '../../ORF_files/Touchon_genomes_25/ORFs_' + genome + '/protein_ORFs_' + genome + '_filtered.fasta'
		orf_folder = '../../ORF_files/Touchon_genomes_25/ORFs_' + genome + '/'
		isExist = os.path.exists(orf_folder)
		if not isExist:
			os.makedirs(orf_folder)

		embedding_path = '../../embeddings/Touchon_genomes_25/genome_embeddings_' + genome + '/'
		isExist = os.path.exists(embedding_path)
		if not isExist:
			os.makedirs(embedding_path)

		list_scripts = [#'find_ORFs.py ' + genome,
						#'filter_ORFs.py ' + genome,
						'annotate_genome.py ' + genome,
						'extract.py esm1b_t33_650M_UR50S ' + orf_path + ' ' + embedding_path + ' --include mean per_tok',
						'get_semantic_distance ' + genome,
						'get_distance_with_averaged_emb' + genome]
		
		for script in list_scripts: 
			print("Running the following script:", script)
			print("Genome:", genome, "Time:", timestamp(), "\n")
			print()
			os.system('python ' + script)


#--------------------------main--------------------------------------------------

DIV = 5 # number of groups into which the data is diveded that run in parallel (parallel processes)

# list of genomes here:
genome_path = "../../genome_files/Touchon_genomes_25/"
genome_list = os.listdir(genome_path)

# spliting and cleaning the data:
split_array = np.array_split(genome_list, DIV)
genome_array = []

for array in split_array:
	new_array = []

	for genome in array:
		genome = genome.split(".")[0]
		new_array.append(genome)

	genome_array.append(new_array)

# running the scripts in "DIV" number of parallel processes:
if __name__ == '__main__':

	for i in range(0, DIV):
		p = mp.Process(target = run_scripts, args = (genome_array[i],))
		p.start()
	




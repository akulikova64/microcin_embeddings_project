from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys


microcin_path = "../../microcin_files/single_microcin_fastas/"
genome_path = "../../ORF_files/ORFs_10_known_microcins/"
output_path = "../../BLAST_results/10_known_microcins/"

microcin_list = os.listdir(microcin_path)
genome_list = os.listdir(genome_path)

#==========================================================================================
# setting parameters:
#==========================================================================================
e_value = "10"
max_target_seqs = "50"
output_format = "6"
#length_cutoff = "" # in percent of query sequence
#==========================================================================================

# running blastp
for genome in genome_list:
    for microcin in microcin_list:

        # get genome name for folder naming:
        genome_name = genome[0:-6]

        # get specific paths for blast results:
        subject = genome_path + genome
        query =  microcin_path + microcin
        output_path_extended = output_path + genome_name + "/"

        if not os.path.exists(output_path_extended):
            os.makedirs(output_path_extended)

        output = output_path_extended + "blastp_" + microcin

        # run blast:
        command = "blastp -subject " + subject + " -query " + query + " -out " + output + " -evalue " + e_value + " -outfmt " + output_format + " -max_target_seqs " + max_target_seqs
        os.system(command)
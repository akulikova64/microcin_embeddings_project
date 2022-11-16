from Bio import SeqIO
from Bio.Seq import Seq
import sys

# this script translates the orfs made by orfipy

# the following commands were used to find orfs in genome:
# orfipy ./genomes/Ecoli_S88.fna --dna e_coli_orfs.fasta --outdir ./ORFs_ecoli_s88
# orfipy ./genomes/klebsiella_pneumoniae_E492.fna --dna k_pneumoniae_orfs.fasta --outdir ./ORFs_kpneimoniae
# orfipy ../../genome_files/klebsiella_pneumoniae_gene_cluster.fna --dna k_pneumoniae_cluster_orfs.fasta --outdir ../../ORF_files/ORFs_kpneimoniae_cluster
# orfipy ../../genome_files/GCA_021559835.1_ASM2155983v1_genomic_genbank.fna --dna ecoli_nissle_genomic_orfs.fasta --outdir ../../ORF_files/ORFs_ecoli_nissle

input_file = "../../ORF_files/ORFs_ecoli_nissle/ecoli_nissle_genomic_orfs.fasta"
output_file = "../../ORF_files/ORFs_ecoli_nissle/protein_ecoli_nissle_genomic_orfs.fasta"

orfs_fasta = list(SeqIO.parse(input_file, "fasta"))
exceptions = []

with open(output_file, "w") as protein_seq_file:
    for rec in orfs_fasta:
        dna_sequence = str(rec.seq)
        dna_sequence = Seq(dna_sequence)
        orf_name = str(rec.description)

        # translating dna sequence into protein sequence:
        protein_seq = dna_sequence.translate(to_stop=True)
        protein_seq = str(protein_seq)

        # writing to new fasta file:
        if len(protein_seq) <= 150 & len(protein_seq) >= 30:
            protein_seq_file.write(">" + orf_name + "\n")
            protein_seq_file.write(protein_seq + "\n")

print("Done!")
print("Next, use the following command to generate esm1-b embeddings:")
print("python3 extract.py esm1b_t33_650M_UR50S ../../ORF_files/ORFs_ecoli_nissle/protein_ORFs_nissle_mine_filtered.fasta ../../embeddings/genome_embeddings_ecoli_nissle/ --include mean per_tok")

# python3 extract.py esm1b_t33_650M_UR50S ../../microcin_files/45_microcins.fasta ../../microcin_files/45_microcins_emb_esm1b/ --include mean per_tok
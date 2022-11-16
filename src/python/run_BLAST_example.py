import pandas as pd

# Get BLAST results

def load_blast(blastout, hasHeaders=False):
	if hasHeaders:
		blastpDF = pd.read_csv(blastout, sep = "\t")
	else:
		blastpDF = pd.read_csv(blastout, sep = "\t", names = ["qseqid","sseqid","pident","length","mismatch",			"gapopen","qstart","qend","sstart","send","evalue","bitscore"])
	return blastpDF
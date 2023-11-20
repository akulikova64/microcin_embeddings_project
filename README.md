# Semantic search using protein large language models detects class II microcins in bacterial genomes
## This repository contains code for microcin search, analysis and additional data.

distance_csvs/ - csv files that contain embedding distances between 10 known microcin embeddings and ORFs

BLAST_results/ - blast results

analysis/figures/ - contains all plots and figures with distances

src/ - R and python scripts for data processing

<br />

### Scripts for data preparation:

1) python/run_pipeline2.py - pipeline for microcin search. It extracts ORFs, filters ORFs, annotates known microcin ORFs, runs esm1-b to generate embeddings, calculates semantic distance with known microins and averaged microcin embedding. 

2) python/run_blastp.py - Runs protein blast using the 10 known microcins sequences as queries.

<br />

### Scripts for analysis:

1) R/genome_template_emb.Rmd - template for data analysis for a single genome/assembly

2) R/microcin_homology_plots.Rmd - code for Figure 5 in paper


---

CC BY-NC 4.0 This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License (see https://creativecommons.org/licenses/by-nc/4.0/legalcode).

Using this work commercially may require a license to intellectual property owned by the Board of Regents of the University of Texas System. Contact licensing@otc.utexas.edu for inquiries.

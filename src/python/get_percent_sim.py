from Bio import SeqIO
import gemmi
import csv
import sys
import os


def get_percent(seq1, seq2):
  count_similar = 0
  for i, j in zip(seq1, seq2):
    assert len(seq1) == len(seq2)
    if i == j: # gaps also considered a mismatch
      count_similar += 1

  perc_sim = float(count_similar/len(seq1)) * 100
  return perc_sim

def get_alignment(seq1, seq2):
    """
    Returns an alignment and its score
    """

    alignment = gemmi.align_string_sequences(list(seq1), list(seq2), [])

    return alignment.formatted(seq1, seq2)

#---------------main-------------------------
input_path = "../../microcin_files/Microcins_Known.fasta"
output_path = "../../analysis/10_microcin_perc_sim.csv" 


with open(output_path, 'w', newline='\n', encoding='utf-8') as csv_file:
  writer = csv.writer(csv_file)
  writer.writerow(['seq1', 'seq2', 'percent_sim'])

  microcin_list = list(SeqIO.parse(input_path, "fasta"))
  num_microcins = len(microcin_list)

  for i in range(0, num_microcins - 1):
    for j in range(i + 1, num_microcins):

        name1 = microcin_list[i].name.split("|")[0][0:-3]
        name2 = microcin_list[j].name.split("|")[0][0:-3]

        alignment = get_alignment(str(microcin_list[i].seq), str(microcin_list[j].seq))
        seq1 = alignment.split("\n")[0]
        seq2 = alignment.split("\n")[2]

        percent_sim = get_percent(str(seq1), str(seq2))
    
        writer.writerow([name1, name2, percent_sim])
        
print("Finished writing to csv! :-)")

import sys
import torch


# prints pytorch file
embedding = torch.load("Microcins_Known_emb_esm1b/E492_sp|Q9Z4N4|MCEA_KLEPN.pt")
#print(embedding)

#Note: last layer in esm1-b is layer 33.

seq = "MREISQKDLNLAFGAGETDPNTQLLNDLGNNMAWGAALGAPGGLGSAALGAAGGALQTVGQGLIDHGPVNVPIPVLIGPSWNGSGSGYNSATSSSGSGS"
print("--Sequence length:",len(seq), "\n")

# representations
print("--Residue embeddings:")
print("Total number of residue embeddings:", len(embedding['representations'][33]))
print("Note: total number of embeddings should be equal to sequence length.")
print("Length of residue embedding vector:", len(embedding['representations'][33][0]), "\n")
   
# mean representations
print("--Averaged embeddings:")
print("Length of averaged embeddings:", len(embedding['mean_representations'][33]))
print("Note: should equal to the length of a single embedding")

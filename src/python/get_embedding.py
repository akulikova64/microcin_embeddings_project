from bio_embeddings.embed import SeqVecEmbedder

embedder = SeqVecEmbedder()
embedding = embedder.embed("MRELDREELNCVGGAGDPLADPNSQIVRQIMSNAAWGAAFGARGGLGGMAVGAAGGVTQTVLQGAAAHMPVNVPIPKVPMGPSWNGSKG")
print(embedding)
import pickle as pkl
import pandas as pd


with open("../../info/pfam_metadata.pkl", "rb") as f:
    file = pkl.load(f)
     
df = pd.DataFrame(file)

df.to_csv(r'../../info/pfam_metadata.csv')


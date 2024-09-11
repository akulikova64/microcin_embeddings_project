import pandas as pd
from sklearn.metrics import silhouette_samples
import numpy as np

def calculate_silhouette_scores(data, group_labels):
    """
    Calculate silhouette scores for each data point based on group labels.

    Parameters:
    data (list of lists or numpy array): Data points.
    group_labels (list): Group labels for each data point (e.g., 'microcins' or 'non-microcins').

    Returns:
    silhouette_scores (list): Silhouette scores for each data point.
    """

    silhouette_scores = silhouette_samples(np.array(data), group_labels)
    return silhouette_scores


#------------------------------------ main ---------------------------------------------------

# Load data from CSV file
df = pd.read_csv("../../analysis/microcin_emb_data.csv")

# Extract group labels and data
group_labels = df["group"].tolist()
data = df.iloc[:, 3:].values  # Columns "X1" through "X1280" contain data

# Calculate silhouette scores
silhouette_scores = calculate_silhouette_scores(data, group_labels)

# Add silhouette scores as a new column to the DataFrame
df["silhouette"] = silhouette_scores

# Save DataFrame back to CSV file
df.to_csv("../../analysis/microcin_emb_data.csv_with_silhouette.csv", index=False)
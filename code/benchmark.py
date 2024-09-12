from sklearn.neighbors import NearestNeighbors
import numpy as np
from numpy import genfromtxt
import pandas as pd
from scipy.spatial import distance_matrix
import os

def mixing(indices, num_cells, k):
    deviation = []
    for i in range(num_cells*2):
        count = sum(1 for neighbor in indices[i, 1:k+1] if neighbor < num_cells)
        mixing_ratio = count / k
        deviation.append(abs(0.5-mixing_ratio))
    # Calculate the average percentage
    deviation_avg = np.mean(deviation)

    return deviation_avg

def foscttm(indices, num_cells):
    dist = distance_matrix(indices[:num_cells, :], indices[num_cells:, :])
    assert dist.shape == (num_cells, num_cells)
    foscttm_x = (dist < np.expand_dims(np.diag(dist), axis=1)).mean(axis=1)
    foscttm_y = (dist < np.expand_dims(np.diag(dist), axis=0)).mean(axis=0)
    fracs = []
    for i in range(len(foscttm_x)):
        fracs.append((foscttm_x[i] + foscttm_y[i]) / 2)
    return np.mean(fracs).round(4)

def get_nn(embed, k):
    embed = embed[1:embed.shape[0],1:]
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='auto').fit(embed)
    distances, indices = nbrs.kneighbors(embed)
    return indices

def read_csv(file_path):
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, index_col=0)
        return df
    else:
        return pd.DataFrame()


"""
embed: low-dimensional embeddings
method: method_name
tissue: tissue origin
output: string of output directory
mix_min: minimum number for top percentage for mix metric, float percentage(0.01)
mix_max: maximum number for top percentage for mix metric
prox_min: min for prox metric
prox_max: max for prox metric
"""
def run_benchmark(embed_filepath, tissue, output):
    
    filepath_output = output + "/" + tissue + "/"
    results_prox = read_csv(filepath_output+"foscttm.csv")
    directory = embed_filepath +"/"+tissue + "/"
    os.chdir(directory)
    # List all files in the directory
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    for file in files:
        method = os.path.splitext(file)[0]
        print(method)

        embed = genfromtxt(file, delimiter=',')
        num_cells = int(embed.shape[0]/2)
        indices = get_nn(embed, int(num_cells*0.01))
        print("FOSCTTM:")
        temp = foscttm(embed[1:embed.shape[0],1:], num_cells)
        results_prox.loc[method, "foscttm"] = temp
        results_prox.loc[method, "num_cells"] = num_cells
        print(temp)
    results_prox.to_csv(filepath_output + "foscttm.csv", index = True)
 
run_benchmark("/users/tfurutan/integration/subsample/coembed/", "pancreas", output = "/users/tfurutan/integration/benchmark/subsample/")
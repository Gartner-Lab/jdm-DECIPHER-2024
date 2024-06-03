# conda activate div-idx
# python3 plot_div_index.py

# TODO: update so that obj_name, dir are arguments
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from collections import Counter
from sklearn.neighbors import NearestNeighbors

import os
#os.chdir("/krummellab/data1/erflynn/premier/jdm/data/sobj_GSE135779/merged_df_proc2")
os.chdir("/krummellab/data1/erflynn/premier/jdm/data/jdm/merged_jdm/")

obj_name = "merged_raw_sample_formatted"
x = pd.read_csv('{}.tsv'.format(obj_name), sep='\t', header=0, index_col=0)
x = x[['LIBRARY', 'UMAP_1', 'UMAP_2']].copy()

y = np.array(list(zip(x['UMAP_1'], x['UMAP_2'])))
nbrs = NearestNeighbors(n_neighbors=int(x.shape[0]**0.5), algorithm='ball_tree').fit(y)
distances, indices = nbrs.kneighbors(y)

threshold = float(np.quantile(distances, q=[0.9]))
libraries = np.array(list(x['LIBRARY']))
library_avgs = Counter(x['LIBRARY'])
library_avgs = {k: v/x.shape[0] for k, v in library_avgs.items()}

def dscore(nbh, avgs):
  temp = Counter(nbh)
  temp = Counter({k: v/len(nbh) for k, v in temp.items()})
  return 0.5 * sum([abs(avgs[k] - temp[k]) for k in avgs])

#x['library_dscore'] = [dscore(libraries[i[d<=threshold]], avgs=library_avgs) for d, i in zip(distances, indices)]


x['ND'] = [len(set(libraries[i[d<=threshold]])) for d, i in zip(distances, indices)]
plt.figure(figsize=(10, 10))
#cm=sns.color_palette("mako", as_cmap=True)
x.plot.scatter(x='UMAP_1', y='UMAP_2', c='ND', s=0.1, vmin=0, vmax=len(set(libraries)), cmap="mako") 
plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
plt.axis('off')
plt.savefig('{}_batch_plot.pdf'.format(obj_name), bbox_inches='tight', pad_inches=0.0)
plt.close()

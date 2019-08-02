import numpy as np

def gene_means(g,genes,X):
    g_idx = np.where(genes == g)[0][0]
    pos_val = X[g_idx, pos]
    neg_val = X[g_idx, neg]
    return
    pos_list.append(np.mean(pos_val))
    neg_list.append(np.mean(neg_val))

import numpy as np
import colorsys
from scipy.stats import spearmanr

def Generate_Color_Dict(labels):
    N = len(np.unique(labels))
    HSV_tuples = [(x * 1.0 / N, 1.0, 0.5) for x in range(N)]
    np.random.shuffle(HSV_tuples)
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    color_dict = dict(zip(np.unique(labels), RGB_tuples))
    return color_dict

def Gene_Corr_func(c,genes,X,val_g):
    idx_c = np.where(genes == c)[0][0]
    val_c = X[idx_c]
    corr, _ = spearmanr(val_g, val_c)
    return corr

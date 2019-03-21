### set environment
import numpy  as np
import pandas as pd
import glob, os, re

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.manifold.t_sne   import _joint_probabilities
dat_dir = "/data/clintko/SMPK"

### function declaration
def get_similarity(X, perplexity = 30, verbose = False):
    """calculate the similarity matrix Pi|j + Pj|i
    ref: https://www.oreilly.com/learning/an-illustrated-introduction-to-the-t-sne-algorithm
    """
    # Pairwise distances between all data points.
    D = pairwise_distances(X, squared=True)    
    # Similarity with variable sigma.
    P = _joint_probabilities(D, perplexity, verbose)
    return P

def get_sample(dat, idx_sample):
    """get the sample from the data"""
    X = dat.loc[(dat["sample"] == idx_sample), ("tsne1", "tsne2")]
    return X

def process_data(dat, idx_sample):
    """get a sample from data and get similarity"""
    print("  Sample:", idx_sample)
    X = get_sample(dat, idx_sample)
    P = get_similarity(X, perplexity = 30, verbose = False)
    print("    Check P dim:", P.shape)
    return P

if __name__ == "__main__":
    ### read in raw / uncalibrated data
    print("==========================")
    print("raw")
    y = pd.read_csv(
        os.path.join(
            dat_dir, 
            "ep8cs_tsne_raw.txt"),
        header=None,
        names=["sample", "tsne1", "tsne2"],
        sep = "\t") 
    
    ### calculate similarity matrix for each sample in raw data
    mat = np.array([process_data(y, idx) for idx in range(1, 18+1)]).T
    print("Check Final dimension:", mat.shape)
    np.save(os.path.join(dat_dir, "ep8cs_tsne_similarity_raw.npy"), mat)
    
    ### process calibrated data
    fnames = glob.glob(os.path.join(dat_dir, "ep8cs_dat_cal_prior*"))
    priors = [re.findall("prior(.*).txt", x)[0] for x in fnames]
    for prior in priors:
        
        ### read in calibrated data
        print("\n==========================")
        print("Prior:", prior)
        #"""
        y = pd.read_csv(
            os.path.join(
                dat_dir, 
                "ep8cs_tsne_cal_prior" + prior + ".txt"),
            header=None,
            names=["sample", "cluster", "tsne1", "tsne2"],
            sep = "\t")
        #"""
        
        ### calculate similarity matrix for each sample in calibrated data
        mat = np.array([process_data(y, idx) for idx in range(1, 18+1)]).T
        print("Check Final dimension:", mat.shape)
        np.save(
            os.path.join(
                dat_dir, 
                "ep8cs_tsne_similarity_cal_" + "prior" + prior + ".npy"
            ), mat)
    
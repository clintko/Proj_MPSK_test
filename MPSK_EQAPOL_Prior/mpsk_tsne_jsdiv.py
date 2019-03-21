### python
import numpy as np
from scipy.stats            import entropy
from scipy.spatial.distance import squareform
import glob, os, re

### cython
### ref: https://github.com/cython/cython/tree/master/pyximport
import pyximport
pyximport.install()
import jsdiv

### file path & name
dat_dir = "/data/clintko/SMPK"
fout    = "ep8cs_tsne_jsdiv"

if __name__ == "__main__":
    ### calculate the pairwise jsdiv for raw
    print("==========================")
    print("get jsdiv matrix: raw")
    print("    Read file")
    p_raw = np.load(os.path.join(dat_dir, "ep8cs_tsne_similarity" + "_raw.npy"))
    DIM = p_raw.shape[1]
    
    print("    Calculate raw x raw")
    #"""    
    mat_jsdiv = np.zeros((DIM, DIM))
    jsdiv.pairwise_jsdiv_cy(mat_jsdiv, p_raw, p_raw)
    mat_jsdiv = mat_jsdiv * 2
    np.save(
        os.path.join(dat_dir, fout + "_raw.npy"), 
        mat_jsdiv)
    #"""
    
    ### get priors from file names
    fnames = glob.glob(os.path.join(dat_dir, "ep8cs_dat_cal_prior*"))
    priors = [re.findall("prior(.*).txt", x)[0] for x in fnames]
    
    ### calculate the pairwise jsdiv for each prior
    fname  = "ep8cs_tsne_similarity_cal_prior"
    for prior in priors:
        print("==========================")
        print("get jsdiv matrix: prior:" + prior)
        
        print("    Read file")
        p_cal = np.load(os.path.join(dat_dir, fname + prior + ".npy"))
        
        print("    Calculate cal x cal")
        #"""    
        mat_jsdiv = np.zeros((DIM, DIM))
        jsdiv.pairwise_jsdiv_cy(mat_jsdiv, p_cal, p_cal)
        mat_jsdiv = mat_jsdiv * 2
        np.save(
            os.path.join(
                dat_dir, 
                fout + "_cal_prior" + prior + ".npy"), 
            mat_jsdiv)
        #"""    
        
        print("    Calculate raw x cal")
        #"""    
        mat_jsdiv = np.zeros((DIM, DIM))
        jsdiv.pairwise_jsdiv_cy(mat_jsdiv, p_raw, p_cal)
        mat_jsdiv = mat_jsdiv * 2
        np.save(
            os.path.join(
                dat_dir, 
                fout + "_raw_cal_prior" + prior + ".npy"), 
            mat_jsdiv)
        #"""    
    
    
    
    
    
    
    
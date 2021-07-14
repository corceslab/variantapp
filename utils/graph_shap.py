import numpy as np
import pandas as pd
import h5py
from matplotlib import pyplot as plt

def load_hdf5(counts_hdf5, profile_hdf5):
    counts = h5py.File(counts_hdf5, 'r')
    profile = h5py.File(profile_hdf5, 'r')
    return counts['coords_chrom'], counts['coords_start'], counts['coords_end'], counts['hyp_scores'], counts['input_seqs'], \
        profile['coords_chrom'], profile['coords_start'], profile['coords_end'], profile['hyp_scores'], profile['input_seqs']

def get_imp(scores, seqs, start, end):
    scores = np.asarray(scores)
    seqs = np.asarray(seqs)
    vals = np.multiply(scores, seqs)
    psA = vals[start:end, 0]
    psC = vals[start:end, 1]
    psG = vals[start:end, 2]
    psT = vals[start:end, 3]
    return psA, psC, psG, psT

def gen_graph(sA, sC, sG, sT, entry, startpos, trim_start, trim_end, chrom):
    plt.switch_backend('Agg')
    ind = np.arange(startpos + trim_start, startpos + trim_end)
    plt.title("Counts Importance Scores")
    s = "Bases on " + chrom.decode("utf-8")
    plt.xlabel(s)
    plt.ylabel("Scores (A=green/C=blue/G=yellow/T=red")
    plt.bar(ind, sA, color = '#4D7F1E')
    plt.bar(ind, sC, color = '#1A17E1')
    plt.bar(ind, sG, color = '#D69824')
    plt.bar(ind, sT, color = '#D40603')
    plt.savefig('../static/images/' + chrom.decode("utf-8")  + "_" + str(startpos) + "_" + str(entry) + '.png')

def gen_delta(counts_hdf5, profile_hdf5, startpos, trim_start, trim_end, chrom):
    c_chrom, c_start, c_end, c_scores, c_seqs, p_chrom, p_start, p_end, p_scores, p_seqs = load_hdf5(counts_hdf5, profile_hdf5)
    start, end = 1037, 1077
    csA1, csC1, csG1, csT1 = get_imp(c_scores[0], c_seqs[0], start, end)
    csA2, csC2, csG2, csT2 = get_imp(c_scores[1], c_seqs[1], start, end)
    plt.switch_backend('Agg')
    ind = np.arange(startpos + trim_start, startpos + trim_end)
    plt.title("Delta Scores Graph")
    s = "Bases on " + chrom.decode("utf-8")
    plt.xlabel(s)
    plt.ylabel("Scores (A=green/C=blue/G=yellow/T=red")
    plt.bar(ind, csA2-csA1, color = '#4D7F1E')
    plt.bar(ind, csC2-csC1, color = '#1A17E1')
    plt.bar(ind, csG2-csG1, color = '#D69824')
    plt.bar(ind, csT2-csT1, color = '#D40603')
    plt.savefig('../static/images/'+ chrom.decode("utf-8") + "_" + str(startpos) + 'DeltaGraph.png')

def vis_shap(counts_hdf5, profile_hdf5, entry):
    c_chrom, c_start, c_end, c_scores, c_seqs, p_chrom, p_start, p_end, p_scores, p_seqs = load_hdf5(counts_hdf5, profile_hdf5)
    start, end = 1037, 1077
    for ind in entry:
        csA, csC, csG, csT = get_imp(c_scores[ind], c_seqs[ind], start, end)
        psA, psC, psG, psT = get_imp(p_scores[ind], p_seqs[ind], start, end)
        gen_graph(csA, csC, csG, csT, ind, c_start[ind], start, end, c_chrom[ind])
    gen_delta(counts_hdf5, profile_hdf5, c_start[0], start, end, c_chrom[0])


if __name__ == '__main__':
    vis_shap('../shap/counts_scores.h5', '../shap/profile_scores.h5', [0, 1])


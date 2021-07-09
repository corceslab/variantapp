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

def gen_graph(sA, sC, sG, sT, entry, startpos, trim_start, trim_end, chrom, isprofile):
    plt.switch_backend('Agg')
    ind = np.arange(startpos + trim_start, startpos + trim_end)
    if isprofile:
        plt.title("Profile Importance Scores")
    else:
        plt.title("Counts Importance Scores")
    s = "Bases on " + chrom.decode("utf-8")
    plt.xlabel(s)
    plt.ylabel("Scores (A=green/C=blue/G=yellow/T=red")
    plt.bar(ind, sA, color = '#4D7F1E')
    plt.bar(ind, sC, color = '#1A17E1')
    plt.bar(ind, sG, color = '#D69824')
    plt.bar(ind, sT, color = '#D40603')
    if isprofile:
        plt.savefig('static/images/ATAC_importances/ATACprofileimp' + str(entry) + '.png')
    else:
        plt.savefig('static/images/ATAC_importances/ATACcountsimp' + str(entry) + '.png')

def vis_shap(counts_hdf5, profile_hdf5, entry):
    c_chrom, c_start, c_end, c_scores, c_seqs, p_chrom, p_start, p_end, p_scores, p_seqs = load_hdf5(counts_hdf5, profile_hdf5)
    start, end = 1007, 1107
    csA, csC, csG, csT = get_imp(c_scores[entry], c_seqs[entry], start, end)
    psA, psC, psG, psT = get_imp(p_scores[entry], p_seqs[entry], start, end)
    gen_graph(csA, csC, csG, csT, entry, c_start[entry], start, end, c_chrom[entry], False)
    gen_graph(psA, psC, psG, psT, entry, p_start[entry], start, end, p_chrom[entry], True)

if __name__ == '__main__':
    for i in range(20):
        vis_shap('../ATAC/shap/counts_scores.h5', '../ATAC/shap/profile_scores.h5', i)
    #vis_shap('../ENCSR000EGM/shap/counts_scores.h5', '../ENCSR000EGM/shap/profile_scores.h5', 0)
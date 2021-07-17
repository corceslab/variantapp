import numpy as np
import pandas as pd
import h5py
from modisco.visualization import viz_sequence

def load_hdf5(counts_hdf5):
    counts = h5py.File(counts_hdf5, 'r')
    return counts['coords_chrom'], counts['coords_start'], counts['coords_end'], \
        counts['hyp_scores'], counts['input_seqs']

def get_imp(scores, seqs, start, end):
    scores = np.asarray(scores)
    seqs = np.asarray(seqs)
    vals = np.multiply(scores, seqs)
    return vals[start:end]

def gen_graphs(peaks_df, one_hot_sequences, hyp_shap_scores):
    c_chrom = peaks_df['chrom']
    c_start = peaks_df['start']
    c_end = peaks_df['end']
    c_seqs = one_hot_sequences
    c_scores = hyp_shap_scores
    start, end = 1042, 1072
    noneffect_scores = get_imp(c_scores[0], c_seqs[0], start, end)
    effect_scores = get_imp(c_scores[1], c_seqs[1], start, end)
    delta_scores = effect_scores-noneffect_scores
    title1 = "Effect: " + c_chrom[0] + " [" + str(c_start[0] + start) + ", " + str(c_start[0]+end)+ "]"
    title2 = "Noneffect: " + c_chrom[1] + " [" + str(c_start[1] + start) + ", " + str(c_start[1]+end)+ "]"
    title3 = "Delta: " + c_chrom[1] + " [" + str(c_start[1] + start) + ", " + str(c_start[1]+end)+ "]"
    viz_sequence.plot_weights(array=effect_scores, title=title1, filepath='static/images/app/effect.png')
    viz_sequence.plot_weights(array=noneffect_scores, title=title2, filepath='static/images/app/noneffect.png')
    viz_sequence.plot_weights(array=delta_scores, title=title3, filepath='static/images/app/delta.png')
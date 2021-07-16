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

def gen_graphs(counts_hdf5):
    c_chrom, c_start, c_end, c_scores, c_seqs = load_hdf5(counts_hdf5)
    start, end = 1042, 1072
    noneffect_scores = get_imp(c_scores[0], c_seqs[0], start, end)
    effect_scores = get_imp(c_scores[1], c_seqs[1], start, end)
    delta_scores = effect_scores-noneffect_scores
    viz_sequence.plot_weights(array=effect_scores, filepath='static/images/app/effect.png')
    viz_sequence.plot_weights(array=noneffect_scores, filepath='static/images/app/noneffect.png')
    viz_sequence.plot_weights(array=delta_scores, filepath='static/images/app/delta.png')
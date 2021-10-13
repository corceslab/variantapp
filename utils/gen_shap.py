import h5py
import json
import numpy as np
import pandas as pd
import pyBigWig
import pysam
import shap
import tensorflow as tf
import tensorflow_probability as tfp
import time
import io
import csv

from modisco.visualization import viz_sequence
from basepairmodels.cli.bpnetutils import *
from basepairmodels.cli.shaputils import *
from basepairmodels.cli.losses import MultichannelMultinomialNLL
from mseqgen.sequtils import one_hot_encode
from mseqgen.utils import gaussian1D_smoothing
from tensorflow.keras.utils import CustomObjectScope
from tensorflow import compat
from tensorflow.keras.models import load_model

from utils.load_model import load


def insert_variant(seq, allele, position):
    left, right = seq[:position-1], seq[position:]
    return left + allele + right


def shap_scores(model, peaks_df):
    # peaks_df['start'] = peaks_df['st'] + peaks_df['summit'] - (2114 // 2)
    # peaks_df['end'] = peaks_df['st'] + peaks_df['summit'] + (2114 // 2)
    num_peaks = peaks_df.shape[0]
    fasta_ref = pysam.FastaFile('reference/hg38.genome.fa')  #..
    control_bigWigs = []
    bias_counts_input = np.zeros((num_peaks, 1))
    bias_profile_input = np.zeros((num_peaks, 1000, 2))
    
    sequences = []
    for idx, row in peaks_df.iterrows():
        start = row['start']
        end = row['end']
        peak_loc = 1057
        allele = row['allele']
        seq = fasta_ref.fetch(row['chrom'], start, end).upper()
        seq = insert_variant(seq, allele, peak_loc)        
        if len(seq) != 2114:
            continue
        sequences.append(seq)

    X = one_hot_encode(sequences, 2114)
    print("X shape", X.shape)
        
    # inline function to handle dinucleotide shuffling
    def data_func(model_inputs):
        rng = np.random.RandomState()
        return [dinuc_shuffle(model_inputs[0], 20, rng)] + \
        [
            np.tile(
                np.zeros_like(model_inputs[i]),
                (20,) + (len(model_inputs[i].shape) * (1,))
            ) for i in range(1, len(model_inputs))
        ]
    
    # shap explainer for the counts head
    profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
        ([model.input[0], model.input[1]], 
         tf.reduce_sum(model.outputs[1], axis=-1)),
        data_func, 
        combine_mult_and_diffref=combine_mult_and_diffref)

    counts_shap_scores = profile_model_counts_explainer.shap_values(
        [X, bias_counts_input], progress_message=100)
    
    #print(peaks_df.head(), X, counts_shap_scores[0])
    print("RETURNING 3 FILES")
    return peaks_df, X, counts_shap_scores[0]

def get_imp(scores, seqs, start, end):
    scores = np.asarray(scores)
    seqs = np.asarray(seqs)
    vals = np.multiply(scores, seqs)
    return vals[start:end]

def get_range(shap1, shap2, delta):
    minval = min(np.amin(shap1), np.amin(shap2), np.amin(delta))
    maxval = max(np.amax(shap1), np.amax(shap2), np.amax(delta))
    buffer = 0.1 * (maxval-minval)
    minval-=buffer
    maxval+=buffer
    return minval, maxval

def gen_graphs(peaks_df, one_hot_sequences, hyp_shap_scores):
    c_chrom = peaks_df['chrom']
    c_start = peaks_df['start']
    print(c_chrom.head())
    print(c_start.head())
    c_seqs = one_hot_sequences
    c_scores = hyp_shap_scores
    center = 1056
    diff = 100
    start, end = center - diff, center + diff + 1
    alt_scores = get_imp(c_scores[0], c_seqs[0], start, end)
    ref_scores = get_imp(c_scores[1], c_seqs[1], start, end)
    delta_scores = alt_scores-ref_scores

    # with open("static/csv/data.csv", "ab") as f:
    #     f.write(b"Alternate Importance Scores: \n")
    #     np.savetxt(f, alt_scores, delimiter=",")
    #     f.write(b"\n")
    #     f.write(b"Reference Importance Scores: \n")
    #     np.savetxt(f, ref_scores, delimiter=",")
    #     f.write(b"\n")
    #     f.write(b"Delta Importance Scores: \n")
    #     np.savetxt(f, delta_scores, delimiter=",")

    export = io.StringIO()
    export.write("Alternate Importance Scores:\n")
    csv.writer(export).writerows(alt_scores.tolist())
    export.write("\n")
    export.write("\n")
    export.write("Reference Importance Scores:\n")
    csv.writer(export).writerows(ref_scores.tolist())
    export.write("\n")
    export.write("\n")
    export.write("Delta Importance Scores:\n")
    csv.writer(export).writerows(delta_scores.tolist())
    export.write("\n")
    export.write("\n")

    alt_allele_seq = c_seqs[0][1056]
    alt = 'N'
    if alt_allele_seq[0] == 1:
        alt = 'A'
    elif alt_allele_seq[1] == 1:
        alt = 'C'
    elif alt_allele_seq[2] == 1:
        alt = 'G'
    elif alt_allele_seq[3] == 1:
        alt = 'T'
    
    ref_allele_seq = c_seqs[1][1056]
    ref = 'N'
    if ref_allele_seq[0] == 1:
        ref = 'A'
    elif ref_allele_seq[1] == 1:
        ref = 'C'
    elif ref_allele_seq[2] == 1:
        ref = 'G'
    elif ref_allele_seq[3] == 1:
        ref = 'T'

    minval, maxval = get_range(alt_scores, ref_scores, delta_scores)
    title1 = "Alternate Importance Scores [allele: " + alt + "]"
    title2 = "Reference Importance Scores [allele: " + ref + "]"
    title3 = "Delta: [alt-ref]"
    altshap = viz_sequence.plot_weights(array=alt_scores, title=title1, filepath='static/images/app/altimp.png', minval=minval, maxval=maxval, color="lightsteelblue", figsize=(30, 4))
    refshap = viz_sequence.plot_weights(array=ref_scores, title=title2, filepath='static/images/app/refimp.png', minval=minval, maxval=maxval, color="lightsteelblue", figsize=(30, 4))
    delshap = viz_sequence.plot_weights(array=delta_scores, title=title3, filepath='static/images/app/delta.png', minval=minval, maxval=maxval, color="lightsteelblue", figsize=(30, 4))
    return altshap, refshap, delshap, export



def shap_scores_main(cell_type, peaks_df, nc):
    # tf.compat.v1.disable_eager_execution()
    # peaks_df2, sequences, scores = shap_scores(model, peaks_df)
    # gen_graphs(peaks_df2, sequences, scores)
    tf.compat.v1.disable_eager_execution()
    model = load(cell_type, nc)
    print(peaks_df.head())
    peaks_df['start'] = peaks_df['st'] + peaks_df['summit'] - (2114 // 2)
    peaks_df['end'] = peaks_df['st'] + peaks_df['summit'] + (2114 // 2)
    peaks_df2, sequences, scores = shap_scores(model, peaks_df)
    return gen_graphs(peaks_df2, sequences, scores)

if __name__ == '__main__':
    tf.compat.v1.disable_eager_execution()
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model('../models/C24/model.h5')
    peaks_df = pd.read_csv('../data/peaks/app.bed', sep='\t', header=None, 
                            names=['chrom', 'st', 'allele', 'summit', 'signalValue'])
    peaks_df['start'] = peaks_df['st'] + peaks_df['summit'] - (2114 // 2)
    peaks_df['end'] = peaks_df['st'] + peaks_df['summit'] + (2114 // 2)
    peaks_df2, sequences, scores = shap_scores(model, peaks_df)
    gen_graphs(peaks_df2, sequences, scores)

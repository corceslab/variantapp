import h5py
import json
import numpy as np
import pandas as pd
import pyBigWig
import pysam
import shap
import tensorflow as tf
import tensorflow_probability as tfp
from modisco.visualization import viz_sequence

from basepairmodels.cli.bpnetutils import *
from basepairmodels.cli.shaputils import *
from basepairmodels.cli.losses import MultichannelMultinomialNLL
from mseqgen.sequtils import one_hot_encode
from mseqgen.utils import gaussian1D_smoothing
from tensorflow.keras.utils import CustomObjectScope


def insert_variant(seq, allele, position):
    left, right = seq[:position-1], seq[position:]
    return left + allele + right


def shap_scores(model, peaks_df):
    peaks_df['start'] = peaks_df['st'] + peaks_df['summit'] - (2114 // 2)
    peaks_df['end'] = peaks_df['st'] + peaks_df['summit'] + (2114 // 2)
    num_peaks = peaks_df.shape[0]
    fasta_ref = pysam.FastaFile('../reference/hg38.genome.fa') 
    control_bigWigs = []
    bias_counts_input = np.zeros((num_peaks, 1))
    
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
        return [dinuc_shuffle(model_inputs[0], 1, rng)] + \
        [
            np.tile(
                np.zeros_like(model_inputs[i]),
                (1,) + (len(model_inputs[i].shape) * (1,))
            ) for i in range(1, len(model_inputs))
        ]
    
    # shap explainer for the counts head
    profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
        ([model.input[0], model.input[1]], 
         tf.reduce_sum(model.outputs[1], axis=-1)),
        data_func, 
        combine_mult_and_diffref=combine_mult_and_diffref)

    # explainer for the profile head
    weightedsum_meannormed_logits = get_weightedsum_meannormed_logits(
        model, task_id=0, stranded=True)
    profile_model_profile_explainer = shap.explainers.deep.TFDeepExplainer(
        ([model.input[0], model.input[2]], weightedsum_meannormed_logits),
        data_func, 
        combine_mult_and_diffref=combine_mult_and_diffref)

    counts_shap_scores = profile_model_counts_explainer.shap_values(
        [X, bias_counts_input], progress_message=100)
    
    #print(peaks_df.head(), X, counts_shap_scores[0])
    return peaks_df, X, counts_shap_scores[0]

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

def shap_scores_main(model, peaks_df):
    tf.compat.v1.disable_eager_execution()
    peaks_df2, sequences, scores = shap_scores(model, peaks_df)
    gen_graphs(peaks_df2, sequences, scores)

if __name__ == '__main__':
    tf.compat.v1.disable_eager_execution()
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model('../models/C24/model.h5')
    peaks_df = pd.read_csv('../data/peaks/app.bed', sep='\t', header=None, 
                            names=['chrom', 'st', 'allele', 'summit', 'signalValue'])
    shap_scores(model, peaks_df)

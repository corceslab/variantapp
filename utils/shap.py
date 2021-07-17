import h5py
import json
import numpy as np
import pandas as pd
import pyBigWig
import pysam
import shap
import tensorflow as tf
import tensorflow_probability as tfp

from basepairmodels.cli.argparsers import shap_scores_argsparser
from basepairmodels.cli.bpnetutils import *
from basepairmodels.cli.exceptionhandler import NoTracebackException
from basepairmodels.cli.shaputils import *
from basepairmodels.cli.logger import *
from basepairmodels.cli.losses import MultichannelMultinomialNLL
from mseqgen.sequtils import one_hot_encode
from mseqgen.utils import gaussian1D_smoothing
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope
from load_model import load


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
    
    print(peaks_df.head(), X, counts_shap_scores[0])
    # save the hyp shap scores, one hot sequences & chrom positions
    # to a HDF5 file
    save_scores(peaks_df, X, counts_shap_scores[0], output_fname)    
        
def shap_scores_main(model, peaks_df):
    with CustomObjectScope({'MultichannelMultinomialNLL': 
                            MultichannelMultinomialNLL}):
        shap_scores(model, peaks_df)

if __name__ == '__main__':
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model('../models/C24/model.h5')
    peaks_df = pd.read_csv('../data/peaks/app.bed', sep='\t', header=None, 
                            names=['chrom', 'st', 'allele', 'summit', 'signalValue'])
    shap_scores(model, peaks_df)

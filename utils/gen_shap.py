import numpy as np
import pandas as pd
import pysam
import shap
import tensorflow as tf

from modisco.visualization import viz_sequence
from basepairmodels.cli.bpnetutils import *
from basepairmodels.cli.shaputils import *

import shap_utils as shap_utils #utils.
from load_model import load_chrombpnet, load_mpra_chrombpnet #utils.


def shap_scores(model, X):
    """ Uses DeepSHAP to calculate importance scores for both alleles
    """
    counts_model_input = [model.input[0], model.input[2]]
    counts_input = [X, np.zeros((X.shape[0], 1))]
    
    profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
            (counts_model_input, tf.reduce_sum(model.outputs[1], axis=-1)),
            shuffle_several_times,
            combine_mult_and_diffref=combine_mult_and_diffref)

    counts_shap_scores = profile_model_counts_explainer.shap_values(
            counts_input, progress_message=10)
    
    return counts_shap_scores[0]

def get_imp(scores, seqs, start, end):
    """ Combines importance scores with the one-hot-encoded sequence to find the
        shap scores for the active bases
    """
    scores = np.asarray(scores)
    seqs = np.asarray(seqs)
    vals = np.multiply(scores, seqs)
    return vals[start:end]

def get_range_chrombpnet(shap1, shap2):
    """ Calculates the y range for the importance score graphs for the individual alleles
        with a buffer of 20% from the min and max values
    """
    minval = min(np.amin(shap1), np.amin(shap2))
    maxval = max(np.amax(shap1), np.amax(shap2))
    buffer = 0.2 * (maxval-minval)
    minval-=buffer
    maxval+=buffer
    return minval, maxval

def get_minmax_chrombpnet(shap1):
    """ Gets the y range for the delta score graph that uses an independent scale
        with a buffer of 20% from the min and max values
    """
    minval = np.amin(shap1)
    maxval = np.amax(shap1)
    buffer = 0.2 * (maxval-minval)
    minval-=buffer
    maxval+=buffer
    return minval, maxval

def gen_graphs(X, hyp_shap_scores, sequences):
    """ Graphs the importance scores for both alleles and the delta values
        Calculates the relevant scores to display using the get_imp method
        Returns the three graphs in binary form

        @Joel - you could store either the hyp_shap_scores and the sequence
        or just the postprocessed alt_scores and ref_scores
    """
    center = 1056
    diff = 250 # CHANGED DIFF (was 40)
    start, end = center - diff, center + diff + 1

    alt_scores = get_imp(hyp_shap_scores[0], X[0], start, end)
    ref_scores = get_imp(hyp_shap_scores[1], X[1], start, end)
    delta_scores = alt_scores-ref_scores

    minval, maxval = get_range_chrombpnet(alt_scores, ref_scores)
    mindelta, maxdelta = get_minmax_chrombpnet(delta_scores)
    title1 = "Alternate Importance Scores [allele: " + sequences[0][1056] + "]"
    title2 = "Reference Importance Scores [allele: " + sequences[1][1056] + "]"
    title3 = "Delta (alt-ref)"
    altshap = viz_sequence.plot_weights(array=alt_scores, title=title1, filepath='static/images/app/altimp.png', minval=minval, maxval=maxval, color="lightsteelblue", figsize=(30, 4))
    refshap = viz_sequence.plot_weights(array=ref_scores, title=title2, filepath='static/images/app/refimp.png', minval=minval, maxval=maxval, color="lightsteelblue", figsize=(30, 4))
    delshap = viz_sequence.plot_weights(array=delta_scores, title=title3, filepath='static/images/app/delta.png', minval=mindelta, maxval=maxdelta, color="lightsteelblue", figsize=(30, 4))
    return altshap, refshap, delshap

def shap_scores_main(cell_type, X, sequences):
    """ The main method to generate DeepSHAP importance scores
    """

    # Config environment and load models
    tf.compat.v1.disable_v2_behavior()
    model, biasmodel = load_chrombpnet(cell_type)

    # Generate imporatnce scores
    scores = shap_scores(model, X)
    return gen_graphs(X, scores, sequences)

def shap_scores_mpra_chrombpnet(model, X):
    """ Uses DeepSHAP to calculate importance scores for both alleles
    """
    counts_model_input = model.input
    counts_input = X
    
    profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
            (counts_model_input, tf.reduce_sum(model.outputs[1], axis=-1)),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref)

    counts_shap_scores = profile_model_counts_explainer.shap_values(
            counts_input, progress_message=10)
    
    return counts_shap_scores

def shap_scores_main_mpra_chrombpnet(cell_type, X, sequences):
    """ The main method to generate DeepSHAP importance scores
    """

    # Config environment and load models
    tf.compat.v1.disable_v2_behavior()
    model = load_mpra_chrombpnet(cell_type)

    # Generate imporatnce scores
    scores = shap_scores_mpra_chrombpnet(model, X)
    return gen_graphs(X, scores, sequences)
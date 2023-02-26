import numpy as np
import pandas as pd
import tensorflow as tf
import pysam
import math
import shap
import time

from utils.sequtils import one_hot_encode
import scipy

from utils.load_model import load_chrombpnet
from utils.bpnetutils import *
from utils.shaputils import *
from scipy.spatial.distance import jensenshannon

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def get_jsd(preds, cts, min_tot_cts=10):
    return np.array([scipy.spatial.distance.jensenshannon(x,y) for x,y in zip(preds, cts) \
                     if y.sum()>min_tot_cts])

def insert_variant(seq, allele, position):
    left, right = seq[:position-1], seq[position:]
    return left + allele + right

def get_preds(model, seqs):
    pred_logits, pred_logcts = model.predict([seqs,
                                              np.zeros((len(seqs), model.output_shape[0][1])),
                                              np.zeros((len(seqs), ))],
                                             batch_size=256, verbose=True)
    # print(pred_logits, pred_logcts)
    counts_profile = softmax(pred_logits) * (np.exp(pred_logcts) - 1)
    return counts_profile, pred_logits, pred_logcts

#sum predicted counts within 200bp of the SNP, take log(ratio) of alt/ref
def sc_lfc(ref_track, alt_track, cutoff, diff = 100):
    #print("ref_track: ", ref_track)
    ref_track = ref_track[len(ref_track) // 2 - diff : len(ref_track) // 2 + diff + 1]
    alt_track = alt_track[len(alt_track) // 2 - diff : len(alt_track) // 2 + diff + 1]
    ref = sum(ref_track)
    alt = sum(alt_track)
    if ref<cutoff and alt<cutoff:
        return 0
    # print("ref: ",ref)
    # print("alt: ",alt)
    lfc = math.log(alt / ref)
    # print("lfc: ",lfc)
    return lfc, ref, alt

def gen_importance(cell_type, peaks_df, variant_names):
    fasta_ref = pysam.FastaFile('reference/hg38.genome.fa')

    tf.compat.v1.disable_eager_execution()
    model, model_bias = load_chrombpnet(cell_type)

    sequences = []
    dist_seqs = []
    peak_loc = 1057

    for idx, row in peaks_df.groupby(peaks_df.index // 2):
        seq = fasta_ref.fetch(row['chrom'].iloc[0], row['start'].iloc[0], row['end'].iloc[0]).upper()
        alt_allele = row['allele'].iloc[0]
        ref_allele = row['allele'].iloc[1]
        sequences.append(insert_variant(seq, alt_allele, peak_loc))
        sequences.append(insert_variant(seq, ref_allele, peak_loc))
    
    numseq = len(sequences)
    # print("number of sequences:", numseq)
    X = one_hot_encode(sequences, 2114)
    # print("X shape", X.shape)

    preds, pred_logits, pred_logcts = get_preds(model, X)
    # print(len(preds))

    #get lfc score
    lfc = []
    abs_lfc = []
    jsd = []
    ref_scores = []
    alt_scores = []
    max_alleles = []
    profiles = softmax(pred_logits)
    # print(profiles.shape)
    # print("jsd:", jensenshannon(profiles[0], profiles[1]))
    for i in range(numseq // 2):
        lfc_score, ref_score, alt_score = sc_lfc(preds[2 * i], preds[2 * i + 1], 0)
        lfc.append(lfc_score)
        abs_lfc.append(abs(lfc_score))
        jsd.append(jensenshannon(profiles[2 * i], profiles[2 * i + 1]))
        ref_scores.append(ref_score)
        alt_scores.append(alt_score)
        max_alleles.append(max(ref_score, alt_score))
    # print("jsd: ", jsd)

    output = pd.DataFrame()
    output['rsID'] = variant_names
    output['lfc'] = lfc
    output['abs_lfc'] = abs_lfc
    output['alt_scores'] = alt_scores
    output['ref_scores'] = ref_scores
    output['max_alleles'] = max_alleles
    output['jsd'] = jsd

    return output

def gen_MPRA_preds(cell_type, seqA, seqB, metadata):

    tf.compat.v1.disable_eager_execution()
    model, model_bias = load_chrombpnet(cell_type)

    dist_seqs = []
    peak_loc = 1057

    numseq = len(seqA)
    XA = one_hot_encode(seqA, 2114)
    XB = one_hot_encode(seqB, 2114)

    print("Done one-hot-encoding")
    print("Local time:", time.ctime(time.time()))
    predsA, pred_logitsA, pred_logctsA = get_preds(model, XA)
    print("Done predsA", "Local time:", time.ctime(time.time()))
    predsB, pred_logitsB, pred_logctsB = get_preds(model, XB)
    print("Done predsB", "Local time:", time.ctime(time.time()))

    #get lfc score
    lfc = []
    abs_lfc = []
    jsd = []
    A_scores = []
    B_scores = []
    max_alleles = []
    profiles_A = softmax(pred_logitsA)
    profiles_B = softmax(pred_logitsB)
    # print(profiles.shape)
    # print("jsd:", jensenshannon(profiles[0], profiles[1]))
    for i in range(numseq):
        lfc_score, A_score, B_score = sc_lfc(predsA[i], predsB[i], 0, 135)
        lfc.append(lfc_score)
        abs_lfc.append(abs(lfc_score))
        jsd.append(jensenshannon(profiles_A[i], profiles_B[i]))
        A_scores.append(A_score)
        B_scores.append(B_score)
        max_alleles.append(max(A_score, B_score))
    # print("jsd: ", jsd)

    output = pd.DataFrame()
    output['lfc'] = lfc
    output['abs_lfc'] = abs_lfc
    output['A_scores'] = A_scores
    output['B_scores'] = B_scores
    output['max_alleles'] = max_alleles
    output['jsd'] = jsd
    output['metadata'] = metadata
    output['NN_cluster'] = cell_type
    return output
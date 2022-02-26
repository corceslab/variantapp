import numpy as np
import pandas as pd
import tensorflow as tf
import pysam
import math
import shap

from mseqgen.sequtils import one_hot_encode
from scipy.special import softmax

from utils.load_model import load, load_chrombpnet
from deeplift.dinuc_shuffle import dinuc_shuffle
from basepairmodels.cli.bpnetutils import *
from basepairmodels.cli.shaputils import *
from scipy.spatial.distance import jensenshannon



def insert_variant(seq, allele, position):
    left, right = seq[:position-1], seq[position:]
    return left + allele + right

def load_sequences(peaks_df):
    fasta_ref = pysam.FastaFile('reference/hg38.genome.fa')
    sequences = []
    for idx, row in peaks_df.iterrows():
        start = row['start']
        end = row['end']
        peak_loc = 1057
        allele = row['allele']
        seq = fasta_ref.fetch(row['chrom'], start, end).upper()
        seq = insert_variant(seq, allele, peak_loc)
        if(len(seq)!= 2114):
            continue
        sequences.append(seq)
    return sequences

def merge(pred_logits, pred_logcts):
    profile = softmax(pred_logits)
    # print(profile.shape)
    counts = np.exp(pred_logcts)[0] - 1
    # print("counts: ", counts)
    return profile*counts

#sum predicted counts within 200bp of the SNP, take log(ratio) of alt/ref
def sc_lfc(ref_track, alt_track, cutoff):
    #print("ref_track: ", ref_track)
    diff = 100
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


def gen_score(model, peaks_df, variant_names):

    #printing model info
    print("\nmodel inputs:")
    [print(i.shape, i.dtype) for i in model.inputs]
    print("\nmodel outputs:")
    [print(o.shape, o.dtype) for o in model.outputs]
    # print("\nmodel layers:")
    # [print(l.name, l.input_shape, l.dtype) for l in model.layers]

    #preprocessing sequences
    sequences = load_sequences(peaks_df)
    numseq = len(sequences)
    # print("number of sequences:", numseq)
    X = one_hot_encode(sequences, 2114)
    
    #control setup
    cl = np.zeros((numseq, 1))
    cp = np.zeros((numseq, 1000, 2))
    
    # print(X.shape)
    # print(cl.shape)
    # print(cp.shape)
    
    #predict
    pred_logits, pred_logcts = model.predict([X, cl, cp], batch_size=numseq, verbose=True)
    #print("prediction output shapes: ", pred_logits.shape, pred_logcts.shape)
    #merge
    preds = []
    for i in range(numseq):
        preds.append(merge(pred_logits[i], pred_logcts[i]))

    #get lfc score
    lfc = []
    alts = []
    refs = []
    d_lfc = []
    for i in range(numseq // 2):
        out_lfc, out_alt, out_ref = sc_lfc(preds[2 * i], preds[2 * i + 1])
        lfc.append(lfc)
        alts.append(out_alt)
        refs.append(out_ref)
        d_lfc.append(abs(lfc[i]))

    output = pd.DataFrame()
    output['rsID'] = variant_names
    output['lfc'] = lfc
    output['alts'] = alts
    output['refs'] = refs
    output['d_lfc'] = d_lfc
    
    return output

def shuffle_seq(seq):
    return dinuc_shuffle(seq, 10)

def gen_importance(cell_type, nc, peaks_df, variant_names):
    fasta_ref = pysam.FastaFile('reference/hg38.genome.fa')
    num_peaks = peaks_df.shape[0]*5
    bias_counts_input = np.zeros((num_peaks, 1))
    bias_profile_input = np.zeros((num_peaks, 1000, 2))

    tf.compat.v1.disable_eager_execution()
    model = load_chrombpnet(cell_type)

    sequences = []
    dist_seqs = []
    peak_loc = 1057

    for idx, row in peaks_df.groupby(peaks_df.index // 2):
        seq = fasta_ref.fetch(row['chrom'].iloc[0], row['start'].iloc[0], row['end'].iloc[0]).upper()
        alt_allele = row['allele'].iloc[0]
        ref_allele = row['allele'].iloc[1]
        sequences.append(insert_variant(seq, alt_allele, peak_loc))
        sequences.append(insert_variant(seq, ref_allele, peak_loc))
        dist_seqs.append(insert_variant(seq, alt_allele, peak_loc))
        dist_seqs.append(insert_variant(seq, ref_allele, peak_loc))
        shuffled_seqs = shuffle_seq(seq)
        #print(shuffled_seqs)
        for seq in shuffled_seqs:
            dist_seqs.append(insert_variant(seq, alt_allele, peak_loc))
            dist_seqs.append(insert_variant(seq, ref_allele, peak_loc))
    
    dist_numseq = len(dist_seqs)
    numseq = len(sequences)
    print("number of distribution sequences:", dist_numseq)
    X = one_hot_encode(sequences, 2114)
    dist_X = one_hot_encode(dist_seqs, 2114)
    print("X shape", X.shape)
    print("HERE")

    #background prediction
    dist_cl = np.zeros((dist_numseq, 1))
    dist_cp = np.zeros((dist_numseq, 1000, 2))
    dist_pred_logits, dist_pred_logcts = model.predict([dist_X, dist_cl, dist_cp], batch_size=dist_numseq, verbose=True)
    dist_counts = np.exp(dist_pred_logcts)
    print(type(dist_counts))
    dist_counts = np.sort(dist_counts[:, 0])
    print("dist_counts: ", dist_counts, dist_counts.size)
    cutoff = dist_counts[int(dist_counts.size/10)]
    print(cutoff)

    #variant prediction
    cl = np.zeros((numseq, 1))
    cp = np.zeros((numseq, 1000, 2))
    pred_logits, pred_logcts = model.predict([X, cl, cp], batch_size=numseq, verbose=True)

    #merge
    preds = []
    for i in range(numseq):
        preds.append(merge(pred_logits[i], pred_logcts[i]))
    print(len(preds))

    #get lfc score
    lfc = []
    d_lfc = []
    jsd = []
    profiles = softmax(pred_logits)
    print(profiles.shape)
    # print("jsd:", jensenshannon(profiles[0], profiles[1]))
    for i in range(numseq // 2):
        lfc.append(sc_lfc(preds[2 * i], preds[2 * i + 1], 0)) #cutoff 0
        d_lfc.append(abs(lfc[i]))
        jsd.append(jensenshannon(profiles[2 * i], profiles[2 * i + 1])[0])
    print("jsd: ", jsd)

    output = pd.DataFrame()
    output['rsID'] = variant_names
    output['lfc'] = lfc
    output['d_lfc'] = d_lfc
    # output['alts'] = alts
    # output['refs'] = 
    # output['max_alleles'] = 
    output['jsd'] = jsd

    return output

    # def data_func(model_inputs):
    #     rng = np.random.RandomState()
    #     return [dinuc_shuffle(model_inputs[0], 20, rng)] + \
    #     [
    #         np.tile(
    #             np.zeros_like(model_inputs[i]),
    #             (20,) + (len(model_inputs[i].shape) * (1,))
    #         ) for i in range(1, len(model_inputs))
    #     ]
    
    # # shap explainer for the counts head
    # profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
    #     ([model.input[0], model.input[1]], 
    #      tf.reduce_sum(model.outputs[1], axis=-1)),
    #     data_func, 
    #     combine_mult_and_diffref=combine_mult_and_diffref)

    # counts_shap_scores = profile_model_counts_explainer.shap_values(
    #     [X, bias_counts_input], progress_message=5)
    # print(counts_shap_scores)
    # output = pd.DataFrame()
    # output['rsID'] = variant_names
    # return output

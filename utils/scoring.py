import numpy as np
import pandas as pd
import tensorflow as tf
import pysam
import math

from mseqgen.sequtils import one_hot_encode
from scipy.special import softmax

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
    print(profile.shape)
    counts = np.exp(pred_logcts)[0]
    print("counts: ", counts)
    return profile*counts

#sum predicted counts within 200bp of the SNP, take log(ratio) of alt/ref
def sc_lfc(ref_track, alt_track):
    diff = 100
    ref_track = ref_track[len(ref_track) // 2 - diff : len(ref_track) // 2 + diff + 1]
    alt_track = alt_track[len(alt_track) // 2 - diff : len(alt_track) // 2 + diff + 1]
    ref = sum(ref_track)
    alt = sum(alt_track)
    lfc = math.log(alt / ref)
    return lfc

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
    print("number of sequences:", numseq)
    X = one_hot_encode(sequences, 2114)
    
    #control setup
    cl = np.zeros((numseq, 1))
    cp = np.zeros((numseq, 1000, 2))
    
    print(X.shape)
    print(cl.shape)
    print(cp.shape)
    
    #predict
    pred_logits, pred_logcts = model.predict([X, cl, cp], batch_size=numseq, verbose=True)
    print("prediction output shapes: ", pred_logits.shape, pred_logcts.shape)
    
    #merge
    preds = []
    for i in range(numseq):
        preds.append(merge(pred_logits[i], pred_logcts[i]))

    #get lfc score
    lfc = []
    d_lfc = []
    for i in range(numseq // 2):
        lfc.append(sc_lfc(preds[2 * i], preds[2 * i + 1]))
        d_lfc.append(abs(lfc[i]))

    output = pd.DataFrame()
    output['rsID'] = variant_names
    output['lfc'] = lfc
    output['d_lfc'] = d_lfc
    output = output.sort_values('d_lfc', ascending=False)
    print(output)

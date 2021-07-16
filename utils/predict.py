import numpy as np
import pandas as pd
import tensorflow as tf
import pysam
from basepairmodels.cli.losses import MultichannelMultinomialNLL
from basepairmodels.cli.losses import multinomial_nll
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope
from mseqgen.sequtils import one_hot_encode
from scipy.special import softmax
from matplotlib import pyplot as plt

def insert_variant(seq, allele, position):
    left, right = seq[:position-1], seq[position:]
    return left + allele + right

def load_sequences():
    fasta_ref = pysam.FastaFile('../reference/hg38.genome.fa')
    peaks_df = pd.read_csv('../data/peaks/app.bed', sep='\t', header=None, 
                            names=['chrom', 'st', 'allele', 'summit', 'signalValue'])
    peaks_df['start'] = peaks_df['st'] + peaks_df['summit'] - (2114 // 2)
    peaks_df['end'] = peaks_df['st'] + peaks_df['summit'] + (2114 // 2)
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

def make_batches(sequences):
    X = one_hot_encode(sequences, 2114)
    X1 = np.reshape(X[0], (1, 2114, 4))
    X2 = np.reshape(X[1], (1, 2114, 4))
    control_profile = np.zeros((1, 1000, 2))
    control_logcount = np.zeros(1)
    batch1 = {
        "sequence": X1,
        "control_profile": control_profile,
        "control_logcount": control_logcount
    }
    batch2 = {
        "sequence": X2,
        "control_profile": control_profile,
        "control_logcount": control_logcount
    }
    return batch1, batch2

def postprocess(pred):
    profile = softmax(pred[0])[0]
    counts = np.exp(pred[1])[0][0]
    return profile*counts

def gen_graphs(pred, filepath):
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(20,2))
    plt.title("Prediction")
    plt.xlabel("Bases")
    plt.ylabel("Predicted TF Binding")
    plt.plot(pred[300:700])
    plt.savefig(filepath)

def predict_main():
    sequences = load_sequences()
    batch1, batch2 = make_batches(sequences)
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model('../models/bpnetv1/model.h5')
    preds1 = model.predict(batch1)
    preds2 = model.predict(batch2)
    prediction1 = postprocess(preds1)
    prediction2 = postprocess(preds2)
    gen_graphs(prediction1, '../static/images/app/noneffectpred.png')
    gen_graphs(prediction2, '../static/images/app/effectpred.png')

if __name__ == '__main__':
    predict_main()
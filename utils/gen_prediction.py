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

def get_range(pred1, pred2, delta, st, en):
    minval = min(np.amin(pred1[st:en]), np.amin(pred2[st:en]), np.amin(delta[st:en]))
    maxval = max(np.amax(pred1[st:en]), np.amax(pred2[st:en]), np.amax(delta[st:en]))
    buffer = 0.1 * (maxval-minval)
    minval-=buffer
    maxval+=buffer
    return minval, maxval
    
def gen_graphs(pred, title, filepath, minval, maxval):
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(20,2))
    plt.title(title)
    plt.xlabel("Bases")
    plt.ylabel("Predicted TF Binding")
    plt.xlim([0, 400])
    plt.ylim([minval, maxval])
    plt.plot(pred)
    plt.plot(np.zeros(400), color="gray")
    plt.savefig(filepath)

def predict_main(model, peaks_df):
    sequences = load_sequences(peaks_df)
    batch1, batch2 = make_batches(sequences)
    preds1 = model.predict(batch1)
    preds2 = model.predict(batch2)
    prediction1 = postprocess(preds1)
    prediction2 = postprocess(preds2)
    delta = np.subtract(prediction2, prediction1)
    st, en = 300, 700
    minval, maxval = get_range(prediction1, prediction2, delta, st, en)
    gen_graphs(prediction1[st:en], 'Noneffect Prediction [allele: '+sequences[0][1056]+']', 'static/images/app/noneffectpred.png', minval, maxval)
    gen_graphs(prediction2[st:en], 'Effect Prediction [allele: '+sequences[1][1056]+']', 'static/images/app/effectpred.png', minval, maxval)
    gen_graphs(delta[st:en], 'Delta Prediction Graph', 'static/images/app/deltapred.png', minval, maxval)

if __name__ == '__main__':
    predict_main('../data/peaks/app.bed')
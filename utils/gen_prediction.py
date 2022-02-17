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
from matplotlib.patches import Rectangle

from PIL import Image
import base64
import io
import csv


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

def postprocess_chrombpnet(logits, logcts):
    print("LOGITS:", logits[0])
    print("LOGCTS:", logcts[0][0])
    profile = softmax(logits[0])
    counts = logcts[0][0]
    print("PROFILE:", profile)
    print("COUNTS:", counts)
    return profile * counts

def get_range(pred1, pred2, st, en):
    minval = min(np.amin(pred1[st:en]), np.amin(pred2[st:en]))
    maxval = max(np.amax(pred1[st:en]), np.amax(pred2[st:en]))
    buffer = 0.1 * (maxval-minval)
    minval-=buffer
    maxval+=buffer
    return minval, maxval
    
def gen_graphs(pred, title, filepath, minval, maxval):
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(30,3))
    plt.title(title)
    plt.ylabel("Predicted Counts")
    plt.xlim([0, 399])
    plt.ylim([minval, maxval])
    plt.plot(pred)
    plt.plot(np.zeros(400), color="gray")
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((199 - .5, minval), 1, maxval-minval, facecolor="grey", alpha=0.5))
    return fig_to_img(plt.gcf())
    #plt.savefig(filepath)

def fig_to_img(fig):
    buf = io.BytesIO()
    data = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = Image.open(buf)
    img.save(data, "PNG")
    encoded_img_data = base64.b64encode(data.getvalue())
    return encoded_img_data

def log_full_change(ref_track, alt_track):
    ref_track = ref_track + 0.001
    alt_track = alt_track + 0.001
    track = np.divide(alt_track, ref_track)
    #multiplier = 10000/np.sum(track)
    #track = track * multiplier
    track = np.log2(track)
    lfcmin, lfcmax = np.amin(track), np.amax(track)
    buffer = 0.1*(lfcmax-lfcmin)
    lfcmin-=buffer
    lfcmax+=buffer
    return track, lfcmin, lfcmax

def predict_main(model, peaks_df):
    sequences = load_sequences(peaks_df)
    batch1, batch2 = make_batches(sequences)
    preds1 = model.predict(batch1)
    preds2 = model.predict(batch2)
    prediction1 = postprocess(preds1)
    prediction2 = postprocess(preds2)
    #delta = np.subtract(prediction2, prediction1)
    lfc, lfcmin, lfcmax = log_full_change(prediction1, prediction2)
    
    export = io.StringIO()
    export.write("Alternate Allele Prediction:\n")
    csv.writer(export).writerows(prediction1.tolist())
    export.write("\n")
    export.write("\n")
    export.write("Reference Allele Prediction:\n")
    csv.writer(export).writerows(prediction2.tolist())
    export.write("\n")
    export.write("\n")
    export.write("Log-Full-Change Prediction:\n")
    csv.writer(export).writerows(lfc.tolist())
    export.write("\n")
    export.write("\n")

    # with open("static/csv/data.csv", "ab") as f:
    #     f.write(b"Alternate Allele Prediction: \n")
    #     np.savetxt(f, prediction1, delimiter=",")
    #     f.write(b"\n")
    #     f.write(b"Reference Allele Prediction: \n")
    #     np.savetxt(f, prediction2, delimiter=",")
    #     f.write(b"\n")
    #     f.write(b"Log-Full-Change Prediction: \n")
    #     np.savetxt(f, lfc, delimiter=",")
    #     f.write(b"\n")
    
    st, en = 300, 700
    minval, maxval = get_range(prediction1, prediction2, st, en)
    altpred = gen_graphs(prediction1[st:en], 'Alternate Prediction [allele: '+sequences[0][1056]+']', 'static/images/app/altpred.png', minval, maxval)
    refpred = gen_graphs(prediction2[st:en], 'Reference Prediction [allele: '+sequences[1][1056]+']', 'static/images/app/refpred.png', minval, maxval)
    lfcpred = gen_graphs(lfc[st:en], 'Log Full Change Graph [alt/ref]', 'static/images/app/lfcpred.png', lfcmin, lfcmax)
    #gen_graphs(delta[st:en], 'Delta Prediction Graph', 'static/images/app/deltapred.png', minval, maxval)
    return altpred, refpred, lfcpred, export


def predict_main_chrombpnet(model_chrombpnet, model_bias, peaks_df):
    sequences = load_sequences(peaks_df)
    X = one_hot_encode(sequences, 2114)
    X1 = np.reshape(X[0], (1, 2114, 4))
    X2 = np.reshape(X[1], (1, 2114, 4))
    control_profile = np.zeros((1, 1000, 2))
    control_logcount = np.zeros(1)
    
    altpredbias_logits, altpredbias_logcts = model_bias.predict(X1, batch_size=1, verbose=True)
    refpredbias_logits, refpredbias_logcts = model_bias.predict(X2, batch_size=1, verbose=True)
    print("INFO:", altpredbias_logits, altpredbias_logcts)
    altpredchrom_logits, altpredchrom_logcts = model_chrombpnet.predict([X1, altpredbias_logits, altpredbias_logcts], batch_size=1, verbose=True)
    refpredchrom_logits, refpredchrom_logcts = model_chrombpnet.predict([X2, refpredbias_logits, refpredbias_logcts], batch_size=1, verbose=True)
    print("INFO:", altpredchrom_logits, altpredchrom_logcts)
    altlogits_wobias = altpredchrom_logits - altpredbias_logits
    reflogits_wobias = refpredchrom_logits - refpredbias_logits
    altlogcts_wobias = np.exp(altpredchrom_logcts)-np.exp(altpredbias_logcts)
    reflogcts_wobias = np.exp(refpredchrom_logcts)-np.exp(refpredbias_logcts)
    #delta = np.subtract(prediction2, prediction1)
    altpred = postprocess_chrombpnet(altlogits_wobias, altlogcts_wobias)
    refpred = postprocess_chrombpnet(reflogits_wobias, reflogcts_wobias)
    print(altpred, refpred)
    lfc, lfcmin, lfcmax = log_full_change(altpred, refpred)
    
    st, en = 300, 700
    minval, maxval = get_range(altpred, refpred, st, en)
    altpred = gen_graphs(altpred[st:en], 'Alternate Prediction [allele: '+sequences[0][1056]+']', 'static/images/app/altpred.png', minval, maxval)
    refpred = gen_graphs(refpred[st:en], 'Reference Prediction [allele: '+sequences[1][1056]+']', 'static/images/app/refpred.png', minval, maxval)
    lfcpred = gen_graphs(lfc[st:en], 'Log Full Change Graph [alt/ref]', 'static/images/app/lfcpred.png', lfcmin, lfcmax)
    return altpred, refpred, lfcpred

if __name__ == '__main__':
    predict_main('../data/peaks/app.bed')
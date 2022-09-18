import numpy as np
import pandas as pd
import subprocess
import tensorflow as tf
import io
from PIL import Image
import base64
import sys

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

sys.path.append('/home/ubuntu/variantapp/utils')

from load_model import load_mpra_chrombpnet
from load_seqs import load_sequences
from query_variant import gen_var_df
from gen_prediction import predict_main_mpra_chrombpnet
from gen_shap import shap_scores_main_mpra_chrombpnet

def load_mpra():
    peaks_df = gen_var_df('chr17', 4914028, 'T', 'C')
    peaks_df = peaks_df.append(gen_var_df('chr19', 2102824, 'A', 'G'))
    peaks_df = peaks_df.append(gen_var_df('chr13', 115075910, 'A', 'G'))
    peaks_df = peaks_df.append(gen_var_df('chr17', 44305186, 'C', 'A'))
    peaks_df = peaks_df.append(gen_var_df('chr6', 160103084, 'T', 'C'))
    peaks_df = peaks_df.append(gen_var_df('chr6', 89873055, 'A', 'T'))
    peaks_df = peaks_df.append(gen_var_df('chr11', 407708, 'C', 'A'))
    peaks_df = peaks_df.append(gen_var_df('chr4', 152217226, 'G', 'C'))
    peaks_df = peaks_df.append(gen_var_df('chr15', 41746625, 'G', 'A'))
    peaks_df.reset_index(inplace=True)
    print(peaks_df)
    print(peaks_df.info())
    gen_PDFs('GM12878', peaks_df)

def gen_PDFs(cell_type, peaks_df):
    # iterates through all variants (2 rows of the peaks_df dataframe at a time)
    # generates all predictions and importances scores, which are saved to the images list
    
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        model_chrombpnet = load_mpra_chrombpnet(cell_type)
        print(model_chrombpnet.outputs)
        print("g:")
        print(g)
        # load the sequences and one-hot-encode from the peaks_df pandas dataframe
        h = g.reset_index(drop=True)
        X, sequences = load_sequences(h)
        print(X.shape)
        # generate predictions and importance scores
        altrefpred, lfcpred = predict_main_mpra_chrombpnet(model_chrombpnet, X, sequences)
        altshap, refshap, delshap = shap_scores_main_mpra_chrombpnet(cell_type, X, sequences)
        # merge the prediction and shap score outputs and save
        variantgraph = merge_graphs(altrefpred, lfcpred, altshap, refshap, delshap)
        variantgraph.save('MPRAModelPDFs/'+g.iloc[0]['chrom']+str(g.iloc[0]['start'])+g.iloc[0]['allele']+g.iloc[1]['allele']+'.pdf')

def merge_graphs(altrefpred, lfcpred, altshap, refshap, delshap):
    """ Helper method to merge the five graphs of predictions and importance scores
        into one large image per variant to be displayed by the variant app
    """
    p1 = Image.open(io.BytesIO(base64.b64decode(altrefpred)))
    p3 = Image.open(io.BytesIO(base64.b64decode(lfcpred)))
    s1 = Image.open(io.BytesIO(base64.b64decode(altshap)))
    s2 = Image.open(io.BytesIO(base64.b64decode(refshap)))
    s3 = Image.open(io.BytesIO(base64.b64decode(delshap)))
    variant = Image.new('RGB', (p1.width, p1.height + p3.height + s1.height + s2.height+ s3.height))
    # print("width", p1.width)
    # print("heights", p1.height, p3.height, s1.height, s2.height, s3.height)
    variant.paste(p1, (0, 0))
    variant.paste(p3, (0, p1.height))
    variant.paste(s1, (0, p1.height + p3.height))
    variant.paste(s2, (0, p1.height + p3.height + s1.height))
    variant.paste(s3, (0, p1.height + p3.height + s1.height+ s2.height))
    return variant

if ('__main__'):
    load_mpra()
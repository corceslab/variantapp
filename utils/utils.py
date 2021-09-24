import numpy as np
import pandas as pd
import subprocess
import tensorflow as tf
import time

from basepairmodels.cli.shap import shap_scores
from numpy.core.fromnumeric import _shape_dispatcher

from utils.load_model import load
from utils.query_variant import query_rsID, query_values
from utils.gen_prediction import predict_main
from utils.gen_shap import shap_scores_main
from utils.query_motif import get_motifs

def generate_output_values(cell_type, chrom, position, alt_allele, ref_allele):
    subprocess.call(['sh' ,'utils/reset.sh'])
    model = load(cell_type)
    peaks_df = query_values(chrom, position, alt_allele, ref_allele)
    predict_main(model, peaks_df)
    shap_scores_main(cell_type, peaks_df)
    get_motifs(peaks_df.iloc[0]['chrom'], peaks_df.iloc[0]['st'])

def generate_output_rsID(cell_type, rsID, nc):
    subprocess.call(['sh' ,'utils/reset.sh'])
    model = load(cell_type, nc)
    peaks_df = query_rsID(rsID)
    altpred, refpred, lfcpred = predict_main(model, peaks_df)
    altshap, refshap, delshap = shap_scores_main(cell_type, peaks_df, nc)
    get_motifs(peaks_df.iloc[0]['chrom'], peaks_df.iloc[0]['st'])
    subprocess.call(['sh' ,'utils/export.sh'])
    return altpred, refpred, lfcpred, altshap, refshap, delshap

if __name__ == '__main__':
    #generate_output_values('abc', 'chr1', 35641660, 'A', 'G')
    generate_output_rsID('C24', 'rs181391313')

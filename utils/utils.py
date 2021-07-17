import numpy as np
import subprocess
import tensorflow as tf
# from keras.models import load_model
from utils.query_variant import query_rsID, query_values
from utils.graph_modisco import gen_graphs
from utils.predict import predict_main
from utils.load_model import load
import time

def generate_output_values(cell_type, chrom, position, effect_allele, noneffect_allele):
    subprocess.call(['sh' ,'../scripts/reset.sh'])
    model = load(cell_type)
    query_values(chrom, position, effect_allele, noneffect_allele)
    predict_main('data/peaks/app.bed', model)
    gen_graphs('shap/counts_scores.h5')
    subprocess.call(['sh', '../scripts/shap.sh'])
    # vis_shap('../shap/counts_scores.h5', '../shap/profile_scores.h5', [0, 1])

def generate_output_rsID(cell_type, rsID):
    subprocess.call(['sh' ,'../scripts/reset.sh'])
    model = load(cell_type)
    query_rsID(rsID)
    subprocess.call(['sh', '../scripts/shap.sh'])
    predict_main('data/peaks/app.bed', model)
    gen_graphs('shap/counts_scores.h5')

if __name__ == '__main__':
    #generate_output_values('abc', 'chr1', 35641660, 'A', 'G')
    generate_output_rsID('Microglia', 'rs181391313')

import numpy as np
import subprocess
import tensorflow as tf
# from keras.models import load_model
from graph_shap import vis_shap
from query_rsID import query_rsID, query_values

def generate_output_values(cell_type, chrom, position, effect_allele, noneffect_allele):
    query_values(chrom, position, effect_allele, noneffect_allele)
    # subprocess.call(['sh', '../scripts/predict.sh'])
    subprocess.call(['sh', '../scripts/shap.sh'])
    vis_shap('../shap/counts_scores.h5', '../shap/profile_scores.h5', [0, 1])


if __name__ == '__main__':
    generate_output_values('abc', 'chr1', 35641660, 'A', 'G')

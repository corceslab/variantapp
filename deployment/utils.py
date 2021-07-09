import numpy as np
import random
import json
import re
import subprocess
import tensorflow as tf
from keras.models import load_model
from graph_shap import vis_shap

def generate_output(cell_type, chromosome, position, allele1, allele2):
    subprocess.call(['sh', '../scripts/reset.sh'])
    with open('../ENCSR000EGM/data/peaks.bed', 'w') as peaks:
        pos = int(position)
        peaks.write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chromosome, str(pos-120), str(pos+120), '.', '1000', '.', '0.0', '0.0', '0.01', '120'))
    subprocess.call(['sh', '../scripts/predict.sh'])
    subprocess.call(['sh', '../scripts/shap.sh'])
    vis_shap('../ENCSR000EGM/shap/counts_scores.h5', '../ENCSR000EGM/shap/profile_scores.h5', 0)

if __name__ == '__main__':
    generate_output('abc', 'chr1', 35641660, 'A', 'G')

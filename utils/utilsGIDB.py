import numpy as np
import pandas as pd
import subprocess
import tensorflow as tf
import io
import base64

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

from utils.load_model import load_bpnet, load_cbp
from utils.load_seqs import load_compound
from utils.gen_prediction import predict_main, predict_main_mpra_chrombpnet
from utils.gen_shap import shap_scores_main, shap_scores_main_mpra_chrombpnet
from utils.query_motif import run_motif

def generate_output_compound(cell_type, chrom, center, positions, alleleAlt, alleleRef, outputID):
    """ Main method used to evaluate a variant
    """
    
    #load model according to cell type
    model_chrombpnet = load_cbp(cell_type)

    #load sequences associated with the given variant
    X, sequences = load_compound(chrom, center, positions, alleleAlt, alleleRef)
    
    # verify that the variant was inserted at the right location (indexing)
    # not sure if we need to incorporate this meaningfully on the frontend
    print(sequences[0][1007:1107])
    print(sequences[1][1007:1107])

    # generate predictions and importance scores
    predict_main(model_chrombpnet, X, sequences, outputID)
    shap_scores_main(cell_type, X, sequences, outputID)

    # query relevant TFs from Vierstra's database
    run_motif(chrom, center, outputID)
import numpy as np
import pandas as pd
import subprocess
import tensorflow as tf
import time
import io
from PIL import Image
import base64

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
    export = io.StringIO()
    images = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        altpred, refpred, lfcpred, predexport = predict_main(model, g)
        altshap, refshap, delshap, shapexport = shap_scores_main(cell_type, peaks_df, nc)
        export.write(predexport.getvalue())
        export.write(shapexport.getvalue())
        p1 = Image.open(io.BytesIO(base64.b64decode(altpred)))
        p2 = Image.open(io.BytesIO(base64.b64decode(refpred)))
        p3 = Image.open(io.BytesIO(base64.b64decode(lfcpred)))
        s1 = Image.open(io.BytesIO(base64.b64decode(altshap)))
        s2 = Image.open(io.BytesIO(base64.b64decode(refshap)))
        s3 = Image.open(io.BytesIO(base64.b64decode(delshap)))
        variant = Image.new('RGB', (p1.width, p1.height + p2.height + p3.height + s1.height + s2.height + s3.height))
        print("width", p1.width)
        print("heights", p1.height, p2.height, p3.height, s1.height, s2.height, s3.height)
        variant.paste(p1, (0, 0))
        variant.paste(p2, (0, p1.height))
        variant.paste(p3, (0, p1.height + p2.height))
        variant.paste(s1, (0, p1.height + p2.height + p3.height))
        variant.paste(s2, (0, p1.height + p2.height + p3.height + s1.height))
        variant.paste(s3, (0, p1.height + p2.height + p3.height+ s1.height+ s2.height))
        images.append(variant)
    
    table, motifexport = get_motifs(peaks_df.iloc[0]['chrom'], peaks_df.iloc[0]['st'])
    #table, motifexport = get_motifs(peaks_df)
    #export.write(motifexport.getvalue())

    encoded = io.BytesIO()
    for image in images:
        image.save(encoded, format="PNG")
    encoded_img_data = base64.b64encode(encoded.getvalue())


    #subprocess.call(['sh' ,'utils/export.sh'])
    
    # return encoded_img_data, refpred, lfcpred, altshap, refshap, delshap, table, export
    return encoded_img_data, table

if __name__ == '__main__':
    #generate_output_values('abc', 'chr1', 35641660, 'A', 'G')
    generate_output_rsID('C24', 'rs181391313')

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

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

from utils.load_model import load, load_chrombpnet
from utils.query_variant import query_rsID, query_values
from utils.gen_prediction import predict_main, predict_main_chrombpnet
from utils.gen_shap import shap_scores_main, insert_variant, shap_scores_main_chrombpnet
from utils.query_motif import get_motifs, get_motif
from utils.scoring import gen_score, gen_importance


# from load_model import load
# from query_variant import query_rsID, query_values
# from gen_prediction import predict_main
# from gen_shap import shap_scores_main
# from query_motif import get_motifs


def generate_output_values(cell_type, chrom, position, alt_allele, ref_allele):
    subprocess.call(['sh' ,'utils/reset.sh'])
    model = load(cell_type)
    peaks_df = query_values(chrom, position, alt_allele, ref_allele)
    predict_main(model, peaks_df)
    shap_scores_main(cell_type, peaks_df)
    get_motifs(peaks_df.iloc[0]['chrom'], peaks_df.iloc[0]['st'])

def generate_output_rsID(cell_type, rsID, nc):
    subprocess.call(['sh' ,'utils/reset.sh'])
    peaks_df = query_rsID(rsID)
    print(peaks_df)
    export = io.StringIO()
    images = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        model = load(cell_type, nc)
        altpred, refpred, lfcpred, predexport = predict_main(model, g)
        print("prediction: ", i)
        altshap, refshap, delshap, shapexport = shap_scores_main(cell_type, g, nc)
        print("importance scores: ", i)
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
    
    motiftables = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        input_chrom = str(g.iloc[0]['chrom'])
        input_loc = str(g.iloc[0]['st'])
        motiftables.append(get_motif(input_chrom, input_loc))

    rsIDs = rsID.split(", ")
    print(rsIDs)
    export_images = []
    for i in range(len(images)):
        encoded = io.BytesIO()
        print("i:", i, "rsIDs[i]:", rsIDs[i])
        width = 3000
        height = 2230
        export_image = Image.new('RGB', (width, height), (255, 255, 255))
        I1 = ImageDraw.Draw(export_image)
        roboto = ImageFont.truetype('static/ttf/Roboto-BoldItalic.ttf', 50)
        w, h = I1.textsize(rsIDs[i], font=roboto)
        I1.text(((width-w)/2, 50), rsIDs[i], font=roboto, fill =(0, 0, 0))
        export_image.paste(images[i], (-40, 130))
        export_image.save(encoded, format="PNG")
        export_images.append(base64.b64encode(encoded.getvalue()))

    return export_images, motiftables

def generate_output_rsID_chrombpnet(cell_type, rsID):
    subprocess.call(['sh' ,'utils/reset.sh'])
    peaks_df = query_rsID(rsID)
    print(peaks_df)
    export = io.StringIO()
    images = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        model_chrombpnet, model_bias = load_chrombpnet(cell_type)
        altpred, refpred, lfcpred = predict_main_chrombpnet(model_chrombpnet, model_bias, g)
        print("prediction: ", i)
        altshap, refshap, delshap, shapexport = shap_scores_main_chrombpnet(cell_type, g)
        print("importance scores: ", i)

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
        variant.save('output_pdfs/' + rsID + '.pdf')
    
    motiftables = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        input_chrom = str(g.iloc[0]['chrom'])
        input_loc = str(g.iloc[0]['st'])
        motiftables.append(get_motif(input_chrom, input_loc))

    rsIDs = rsID.split(", ")
    print(rsIDs)
    export_images = []
    for i in range(len(images)):
        encoded = io.BytesIO()
        print("i:", i, "rsIDs[i]:", rsIDs[i])
        width = 3000
        height = 2230
        export_image = Image.new('RGB', (width, height), (255, 255, 255))
        I1 = ImageDraw.Draw(export_image)
        roboto = ImageFont.truetype('static/ttf/Roboto-BoldItalic.ttf', 50)
        w, h = I1.textsize(rsIDs[i], font=roboto)
        I1.text(((width-w)/2, 50), rsIDs[i], font=roboto, fill =(0, 0, 0))
        export_image.paste(images[i], (-40, 130))
        export_image.save(encoded, format="PNG")
        export_images.append(base64.b64encode(encoded.getvalue()))

    return export_images, motiftables

def generate_lfc_ranking(cell_type, rsID, nc):
    variant_names = rsID.split(", ")
    print(variant_names)
    peaks_df = query_rsID(rsID)
    print("peaks_df:\n", peaks_df)
    model = load(cell_type, nc)
    lfcscores = gen_score(model, peaks_df, variant_names)
    return lfcscores

def generate_explain_score(cell_type, rsID, nc):
    variant_names = rsID.split(", ")
    peaks_df = query_rsID(rsID)
    explain_score = gen_importance(cell_type, nc, peaks_df, variant_names)
    return explain_score

if __name__ == '__main__':
    #generate_output_values('abc', 'chr1', 35641660, 'A', 'G')
    generate_output_rsID('C24', 'rs181391313, rs636317', '00')

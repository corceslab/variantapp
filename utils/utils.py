import numpy as np
import pandas as pd
import subprocess
import tensorflow as tf
import io
from PIL import Image
import base64

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

from utils.load_model import load, load_chrombpnet, load_mpra_chrombpnet
from utils.load_seqs import load_sequences, load_compound
from utils.query_variant import query_rsID, query_values, query_values_scoring
from utils.gen_prediction import predict_main, predict_main_mpra_chrombpnet
from utils.gen_shap import shap_scores_main, shap_scores_main_mpra_chrombpnet
from utils.query_motif import get_motif
from scoring.scoringV2 import gen_importance

def generate_inputs(rsID):
    subprocess.call(['sh' ,'utils/reset.sh'])
    peaks_df = query_rsID(rsID)
    return peaks_df

def generate_output_values_chrombpnet(cell_type, peaks_df):
    images = []

    # iterates through all variants (2 rows of the peaks_df dataframe at a time)
    # generates all predictions and importances scores, which are saved to the images list
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        model_chrombpnet, model_bias = load_chrombpnet(cell_type)
        
        # load the sequences and one-hot-encode from the peaks_df pandas dataframe
        X, sequences = load_sequences(g)

        # generate predictions and importance scores
        altrefpred, lfcpred = predict_main(model_chrombpnet, X, sequences)
        altshap, refshap, delshap = shap_scores_main(cell_type, X, sequences)
        
        # merge the prediction and shap score outputs and save
        variantgraph = merge_graphs(altrefpred, lfcpred, altshap, refshap, delshap)
        images.append(variantgraph)
        variantgraph.save('output_pdfs/customvalues.pdf')
    
    # query relevant motifs, getting one data table per variant
    motiftables = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        input_chrom = str(g.iloc[0]['chrom'])
        input_loc = str(g.iloc[0]['st'])
        motiftables.append(get_motif(input_chrom, input_loc))

    
    # postprocessing images and adding a title with the rsID information
    export_images = []
    for i in range(len(images)):
        encoded = io.BytesIO()
        # print("i:", i, "rsIDs[i]:", rsIDs[i])
        width = images[i].width
        height = images[i].height + 100
        export_image = Image.new('RGB', (width, height), (255, 255, 255))
        I1 = ImageDraw.Draw(export_image)
        roboto = ImageFont.truetype('static/ttf/Roboto-BoldItalic.ttf', 50)
        # w, h = I1.textsize(rsIDs[i], font=roboto)
        # I1.text(((width-w)/2, 50), rsIDs[i], font=roboto, fill =(0, 0, 0))
        export_image.paste(images[i], (-40, 130))
        export_image.save(encoded, format="PNG")
        export_images.append(base64.b64encode(encoded.getvalue()))

    return export_images, motiftables

def generate_output_rsID_chrombpnet(cell_type, rsID):
    """ Main method used to evaluate a list of variants
    """
    # resets any files written by the previous run
    subprocess.call(['sh' ,'utils/reset.sh'])

    # gets dataframe with relevant information on the variants
    peaks_df = query_rsID(rsID)
    images = []
    
    # iterates through all variants (2 rows of the peaks_df dataframe at a time)
    # generates all predictions and importances scores, which are saved to the images list
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        model_chrombpnet, model_bias = load_chrombpnet(cell_type)
        
        # load the sequences and one-hot-encode from the peaks_df pandas dataframe
        X, sequences = load_sequences(g)

        # generate predictions and importance scores
        altrefpred, lfcpred = predict_main(model_chrombpnet, X, sequences)
        altshap, refshap, delshap = shap_scores_main(cell_type, X, sequences)
        
        # merge the prediction and shap score outputs and save
        variantgraph = merge_graphs(altrefpred, lfcpred, altshap, refshap, delshap)
        images.append(variantgraph)
        variantgraph.save('output_pdfs/' + rsID + '.pdf')
    
    # query relevant motifs, getting one data table per variant
    motiftables = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        input_chrom = str(g.iloc[0]['chrom'])
        input_loc = str(g.iloc[0]['st'])
        motiftables.append(get_motif(input_chrom, input_loc))

    rsIDs = rsID.split(", ")
    print(rsIDs)
    
    # postprocessing images and adding a title with the rsID information
    export_images = []
    for i in range(len(images)):
        encoded = io.BytesIO()
        # print("i:", i, "rsIDs[i]:", rsIDs[i])
        width = images[i].width
        height = images[i].height + 100
        export_image = Image.new('RGB', (width, height), (255, 255, 255))
        I1 = ImageDraw.Draw(export_image)
        roboto = ImageFont.truetype('static/ttf/Roboto-BoldItalic.ttf', 50)
        w, h = I1.textsize(rsIDs[i], font=roboto)
        I1.text(((width-w)/2, 50), rsIDs[i], font=roboto, fill =(0, 0, 0))
        export_image.paste(images[i], (-40, 130))
        export_image.save(encoded, format="PNG")
        export_images.append(base64.b64encode(encoded.getvalue()))

    return export_images, motiftables

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


# def generate_lfc_ranking(cell_type, rsID, nc):
#     variant_names = rsID.split(", ")
#     print(variant_names)
#     peaks_df = query_rsID(rsID)
#     print("peaks_df:\n", peaks_df)
#     model = load(cell_type, nc)
#     lfcscores = gen_score(model, peaks_df, variant_names)
#     return lfcscores


def generate_explain_score_rankings(cell_type, rsID, vars_df):
    """ Utility method in bulk variant scoring that redirects to scoringV2.py
    """
    variant_names = rsID.split(", ")
    peaks_df = query_values_scoring(vars_df)
    explain_score = gen_importance(cell_type, peaks_df, variant_names)
    return explain_score


def generate_output_compound(cell_type, chrom, center, positions, alleleA, alleleB):
    """ Main method used to evaluate a list of variants
    """
    # resets any files written by the previous run
    subprocess.call(['sh' ,'utils/reset.sh'])

    images = []
    
    model_chrombpnet, model_bias = load_chrombpnet(cell_type)
    X, sequences = load_compound(chrom, center, positions, alleleA, alleleB)
    print(sequences[0][1007:1107])
    print(sequences[1][1007:1107])

    # generate predictions and importance scores
    altrefpred, lfcpred = predict_main(model_chrombpnet, X, sequences)
    altshap, refshap, delshap = shap_scores_main(cell_type, X, sequences)
        
    # merge the prediction and shap score outputs and save
    variantgraph = merge_graphs(altrefpred, lfcpred, altshap, refshap, delshap)
    images.append(variantgraph)
    variantgraph.save('output_pdfs/customcompound.pdf')
    
    # query relevant motifs, getting one data table per variant
    motiftables = []
    for position in positions:
        motiftables.append(get_motif(str(chrom), str(center + position)))

    # postprocessing images and adding a title with the rsID information
    export_images = []
    for i in range(len(images)):
        encoded = io.BytesIO()
        # print("i:", i, "rsIDs[i]:", rsIDs[i])
        width = images[i].width
        height = images[i].height + 100
        export_image = Image.new('RGB', (width, height), (255, 255, 255))
        I1 = ImageDraw.Draw(export_image)
        roboto = ImageFont.truetype('static/ttf/Roboto-BoldItalic.ttf', 50)
        export_image.paste(images[i], (-40, 130))
        export_image.save(encoded, format="PNG")
        export_images.append(base64.b64encode(encoded.getvalue()))

    return export_images, motiftables

def generate_output_values_mpra_chrombpnet(cell_type, peaks_df, chrom, pos, alt_allele, ref_allele):
    images = []

    # iterates through all variants (2 rows of the peaks_df dataframe at a time)
    # generates all predictions and importances scores, which are saved to the images list
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        model_chrombpnet = load_mpra_chrombpnet(cell_type)
        
        # load the sequences and one-hot-encode from the peaks_df pandas dataframe
        X, sequences = load_sequences(g)

        # generate predictions and importance scores
        altrefpred, lfcpred = predict_main_mpra_chrombpnet(model_chrombpnet, X, sequences)
        altshap, refshap, delshap = shap_scores_main_mpra_chrombpnet(cell_type, X, sequences)
        
        # merge the prediction and shap score outputs and save
        variantgraph = merge_graphs(altrefpred, lfcpred, altshap, refshap, delshap)
        images.append(variantgraph)
        variantgraph.save('MPRAModelPDFs/'+chrom+str(pos)+alt_allele+ref_allele+'.pdf')
    
    motiftables = []
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        input_chrom = str(g.iloc[0]['chrom'])
        input_loc = str(g.iloc[0]['st'])
        motiftables.append(get_motif(input_chrom, input_loc))

    # postprocessing images and adding a title with the rsID information
    export_images = []
    for i in range(len(images)):
        encoded = io.BytesIO()
        # print("i:", i, "rsIDs[i]:", rsIDs[i])
        width = images[i].width
        height = images[i].height + 100
        export_image = Image.new('RGB', (width, height), (255, 255, 255))
        I1 = ImageDraw.Draw(export_image)
        roboto = ImageFont.truetype('static/ttf/Roboto-BoldItalic.ttf', 50)
        # w, h = I1.textsize(rsIDs[i], font=roboto)
        # I1.text(((width-w)/2, 50), rsIDs[i], font=roboto, fill =(0, 0, 0))
        export_image.paste(images[i], (-40, 130))
        export_image.save(encoded, format="PNG")
        export_images.append(base64.b64encode(encoded.getvalue()))

    return export_images, motiftables
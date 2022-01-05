from types import ModuleType
from flask import Flask, request
from flask import render_template
from PIL import Image
import base64
import io
from markupsafe import Markup
import shutil

from numpy.core.fromnumeric import var

from utils.form import form_values, form_rsID
from utils.utils import generate_lfc_ranking
from utils.load_SNPs import load_variants

import tensorflow as tf
from tensorflow import compat


if __name__ == '__main__':
    vars_df = load_variants()
    print(vars_df.info())
    vars_df = vars_df[vars_df['Has_ML_prediction'] == True]
    # vars_df = vars_df.drop(vars_df.columns[[1, 2, 3]], axis = 1)
    print(vars_df.info())
    print(vars_df['Has_ML_prediction'].head(30))
    vars_df = vars_df[['SNP_rsID', 'Disease', 'hg38_Chromosome', 'hg38_Position', 'Effect_Allele', 'Noneffect_Allele', 'ML_confidence', 'ML_sig_clusters']]
    vars_df = vars_df.sort_values('ML_confidence', ascending=False)
    vars_df = vars_df[vars_df['Disease']=='AD']
    vars_df = vars_df[vars_df['SNP_rsID'].str[:2]=='rs']
    print(vars_df.info())
    print(vars_df.head(30))
    all_rsIDs = vars_df[['SNP_rsID']].values.tolist()
    rsIDs = ''
    for rsID in all_rsIDs:
        if rsID[0][:2] == 'rs':
            rsIDs += rsID[0] + ", "
    rsIDs = rsIDs[:-2]
    print(rsIDs)


    cell_type = 'C24'
    # rsID = 'rs691331, rs691342, rs691346'
    nc = '10'
    lfc_scores = generate_lfc_ranking(cell_type, rsIDs, nc)
    lfc_scores.reset_index(inplace=True)
    vars_df.reset_index(inplace=True)
    print(lfc_scores)
    vars_df['lfc'] = lfc_scores['lfc']
    vars_df['d_lfc'] = lfc_scores['d_lfc']
    vars_df = vars_df.sort_values('d_lfc', ascending=False)
    print(vars_df)
    
    
    
    vars_df.to_csv('data/output/ADPD_SNPs_score.csv')

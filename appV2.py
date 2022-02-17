from types import ModuleType
from flask import Flask, request
from flask import render_template
from PIL import Image
import base64
import io
from markupsafe import Markup
import shutil

from utils.form import form_values, form_rsID, form_scoring
from utils.utils import generate_output_rsID_chrombpnet, generate_output_values, generate_output_rsID
from utils.query_motif import get_motifs

from utils.utils import generate_explain_score
from utils.load_SNPs import load_variants

import tensorflow as tf
from tensorflow import compat


app = Flask(__name__)

@app.route("/", methods=['GET', 'POST'])

def home():
    """Home page of app with form"""
    app.jinja_env.auto_reload = True
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    form = form_rsID(request.form)
    if request.method == 'POST' and form.validate():
        cell_type = request.form['cell_type']
        rsID = request.form['rsID']
        nc = request.form['nc']

        motifs = []
        graphs, tables = generate_output_rsID_chrombpnet(cell_type, rsID)
        for i in range(len(graphs)):
            graphs[i] = graphs[i].decode('utf-8')
        for table in tables:
            motifs.append(Markup(table))
        ind = list(range(len(graphs)))
        print(ind)
        return render_template('outputV3.html', zip=zip(graphs,motifs))
    return render_template('indexV3.html', form=form)

@app.route('/varianteffectvisualizer')
def vev():
    """Home page of app with form"""
    app.jinja_env.auto_reload = True
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    form = form_rsID(request.form)
    if request.method == 'POST' and form.validate():
        cell_type = request.form['cell_type']
        rsID = request.form['rsID']
        nc = request.form['nc']

        motifs = []
        graphs, tables = generate_output_rsID(cell_type, rsID, nc)
        for i in range(len(graphs)):
            graphs[i] = graphs[i].decode('utf-8')
        for table in tables:
            motifs.append(Markup(table))
        ind = list(range(len(graphs)))
        print(ind)
        return render_template('outputV3.html', zip=zip(graphs,motifs))
    return render_template('indexV3.html', form=form)

@app.route('/variantscoring')
def vs():
    app.jinja_env.auto_reload = True
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    form = form_scoring(request.form)
    nc = '20'
    if request.method == 'POST' and form.validate():
        cluster = request.form['cell_type']
        rsID = request.form['rsID']
        vars_df = load_variants()
        print(vars_df.info())
        vars_df = vars_df[vars_df['Has_ML_prediction'] == True]
        vars_df = vars_df[['SNP_rsID', 'Disease', 'hg38_Chromosome', 'hg38_Position', 'Effect_Allele', 'Noneffect_Allele', 'ML_confidence', 'ML_sig_clusters']]
        vars_df = vars_df.sort_values('ML_confidence', ascending=False)
        vars_df = vars_df[vars_df['Disease']=='AD']
        vars_df = vars_df[vars_df['SNP_rsID'].str[:2]=='rs']
        filtered_vars_df = vars_df[vars_df['ML_sig_clusters'].str.find(cluster) != -1]
        print(vars_df.info())
        print(vars_df.head())
        print(filtered_vars_df.head())
        filtered_vars_df.reset_index(drop=True, inplace=True)
        all_rsIDs = filtered_vars_df[['SNP_rsID']].values.tolist()
        rsIDs = ''
        for rsID in all_rsIDs:
            if rsID[0][:2] == 'rs':
                rsIDs += rsID[0] + ", "
        rsIDs = rsIDs[:-2]
        print(rsIDs)
        lfc_scores = generate_explain_score('C' + cluster, rsIDs, nc)
        lfc_scores.reset_index(drop=True, inplace=True)
        print(lfc_scores)
        scores_df = filtered_vars_df
        scores_df['lfc'] = lfc_scores['lfc']
        scores_df['d_lfc'] = lfc_scores['d_lfc']
        scores_df['jsd'] = lfc_scores['jsd']
        scores_df = scores_df.sort_values(by=['d_lfc', 'jsd'], ascending=[False, False])
        print(scores_df)
        # scores_df = scores_df.sort_values('jsd', ascending=False)
        # print(scores_df)
        # generate_explain_score('C' + cluster, rsIDs, nc)
        vars_df.to_csv('data/output/ADPD_SNPs_score.csv')
    return render_template('indexV3.html', form=form)
# def home():
#     """Home page of app with form"""
#     form = form_values(request.form)
#     if request.method == 'POST' and form.validate():
#         cell_type = request.form['cell_type']
#         chromosome = request.form['chromosome']
#         position = request.form['position']
#         allele1 = request.form['allele1']
#         allele2 = request.form['allele2']
#         generate_output_values(cell_type, chromosome, position, allele1, allele2)
#         return render_template('output.html')
#     return render_template('index_values.html', form=form)

if __name__ == '__main__':
    app.run()
    # generate_output_rsID_chrombpnet('C24', 'rs636317')
from types import ModuleType
from flask import Flask, request
from flask import render_template
from PIL import Image
import base64
import io

from utils.form import form_values, form_rsID
from utils.utils import generate_output_values, generate_output_rsID
from utils.query_motif import get_motifs

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
        altpred, refpred, lfcpred, altshap, refshap, delshap = generate_output_rsID(cell_type, rsID, nc)
        return render_template('outputV2.html', altpred=altpred.decode('utf-8'), refpred=refpred.decode('utf-8'), lfcpred=lfcpred.decode('utf-8'), \
            altshap=altshap.decode('utf-8'), refshap=refshap.decode('utf-8'), delshap=delshap.decode('utf-8'))
    return render_template('index_rsID.html', form=form)

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
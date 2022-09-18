from flask import Flask, request
from flask import render_template
from markupsafe import Markup
from werkzeug.datastructures import MultiDict


from utils.form import form_compound, form_rsID, form_values, form_info, form_chrombpnet
from utils.utils import generate_output_rsID_chrombpnet, generate_output_values_chrombpnet, generate_output_compound, generate_output_values_mpra_chrombpnet
from utils.query_variant import gen_var_df, query_rsID

import tensorflow as tf
from tensorflow import compat


app = Flask(__name__)

@app.route("/", methods=['GET', 'POST'])
def home():
    """Variant Effect Prediction Platform Home Page
    """
    # configs
    app.jinja_env.auto_reload = True
    app.config['TEMPLATES_AUTO_RELOAD'] = True

    # use the 'form_rsID' form to query information from the user
    form = form_rsID(request.form)

    # after the user hits 'submit' - what happens?
    if request.method == 'POST' and form.validate():
        
        # collect relevant information
        cell_type = request.form['cell_type']
        rsID = request.form['rsID']

        motifs = []

        # evaluate the variants
        graphs, tables = generate_output_rsID_chrombpnet(cell_type, rsID)
        
        #preprocess results before outputting to user
        for i in range(len(graphs)):
            graphs[i] = graphs[i].decode('utf-8')
        for table in tables:
            motifs.append(Markup(table))
        
        return render_template('outputV3.html', zip=zip(graphs,motifs))
    return render_template('indexV3.html', form=form)

@app.route("/values", methods=['GET', 'POST'])
def valueEntry():
    form = form_values(request.form)

    if request.method == 'POST' and form.validate():
        cell_type = request.form['cell_type']
        chrom = request.form['chromosome']
        pos = int(request.form['position'])
        alt_allele = request.form['allele1']
        ref_allele = request.form['allele2']

        peaks_df = gen_var_df(chrom, pos, alt_allele, ref_allele)
        graphs, tables = generate_output_values_chrombpnet(cell_type, peaks_df)

        motifs = []

        #preprocess results before outputting to user
        for i in range(len(graphs)):
            graphs[i] = graphs[i].decode('utf-8')
        for table in tables:
            motifs.append(Markup(table))
        
        return render_template('outputV3.html', zip=zip(graphs,motifs))
    return render_template('indexV3V.html', form=form)
    # form = form_values(formdata=MultiDict({'cell_type': }))

@app.route("/compound-variants", methods=['GET', 'POST'])
def compoundVariant():
    form = form_compound(request.form, 5)

    if request.method == 'POST':
        cell_type = request.form['cell_type']
        chrom = request.form['chromosome']
        center = int(request.form['position'])
        numSNPs = int(request.form['num_used'])
        positions = []
        alleleA = []
        alleleB = []
        for entry in form.variants.entries[:numSNPs]:
            positions.append(entry.data['pos'])
            alleleA.append(entry.data['alleleA'])
            alleleB.append(entry.data['alleleB'])
        graphs, tables = generate_output_compound(cell_type, chrom, center, positions, alleleA, alleleB)
        
        motifs = []
        for i in range(len(graphs)):
            graphs[i] = graphs[i].decode('utf-8')
        for table in tables:
            motifs.append(Markup(table))
        
        return render_template('outputV3.html', zip=zip(graphs,motifs))
    return render_template('indexV3C.html', form=form)

@app.route("/info", methods=['GET', 'POST'])
def variantInfo():
    form = form_info(request.form)
    
    if request.method == 'POST' and form.validate():
        rsIDs = request.form['rsIDs']
        peaks_df = query_rsID(rsIDs)
        # construct the html string to display to the user
        html_string = '''
        <html>
        <head><title>rsID Query Output</title></head>
        <link rel="stylesheet" type="text/css" href="/static/css/df_style.css"/>
        <body>
            {table}
        </body>
        </html>
        '''
        s = html_string.format(table=peaks_df.to_html(classes='mystyle', escape=False, index=False))
        return render_template('outputInfo.html', table=Markup(s))
    return render_template('indexV3I.html', form=form)

@app.route("/chrombpnet", methods=['GET', 'POST'])
def chrombpnet():
    form = form_chrombpnet(request.form)

    if request.method == 'POST' and form.validate():
        cell_type = request.form['cell_type']
        chrom = request.form['chromosome']
        pos = int(request.form['position'])
        alt_allele = request.form['allele1']
        ref_allele = request.form['allele2']

        peaks_df = gen_var_df(chrom, pos, alt_allele, ref_allele)
        graphs, tables = generate_output_values_mpra_chrombpnet(cell_type, peaks_df, chrom, pos, alt_allele, ref_allele)

        motifs = []

        #preprocess results before outputting to user
        for i in range(len(graphs)):
            graphs[i] = graphs[i].decode('utf-8')
        for table in tables:
            motifs.append(Markup(table))
        
        return render_template('outputV3.html', zip=zip(graphs,motifs))
    return render_template('indexV3V.html', form=form)

# @app.route("/gen_PDF", methods=['GET', 'POST'])
# def gen_PDF():


if __name__ == '__main__':
    app.run()
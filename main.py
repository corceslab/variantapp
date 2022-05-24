from flask import Flask, request
from flask import render_template
from markupsafe import Markup

from utils.form import form_rsID
from utils.utils import generate_output_rsID_chrombpnet

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

if __name__ == '__main__':
    app.run()
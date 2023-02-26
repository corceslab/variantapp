from flask import Flask, request
from flask import render_template
from markupsafe import Markup
from werkzeug.datastructures import MultiDict
from os.path import exists
from time import strftime, gmtime

from utils.form import form_compound
from utils.utilsGIDB import generate_output_compound
from utils.query_motif import create_table
import keras
import tensorflow as tf
from tensorflow import compat
import pandas as pd


app = Flask(__name__)


# Flask app route to the main ChromBPNet webpage
# Queries chromosome, central position, number of variants
# For each variant, asks for variant position, alleleA, alleleB
@app.route("/", methods=['GET', 'POST'])
def home():
    form = form_compound(request.form, 5)

    if request.method == 'POST':
        print("START", strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        cell_type = request.form['cell_type']
        chrom = request.form['chromosome']
        center = int(request.form['position'])
        numSNPs = int(request.form['num_used'])
        positions = []
        alleleA = []
        alleleB = []
        # outputID = 'm' + cell_type.split()[0] + '.' + chrom + '.pos' + str(center)
        outputID = chrom + str(center)
        for entry in form.variants.entries[:numSNPs]:
            positions.append(entry.data['pos'])
            alleleA.append(entry.data['alleleA'])
            alleleB.append(entry.data['alleleB'])
            # outputID = outputID + '.' + str(entry.data['pos']) + entry.data['alleleA'] + entry.data['alleleB']
            outputID = outputID + str(entry.data['alleleA'] + str(entry.data['alleleB']+cell_type.split()[0]))
        
        header = 'Model Output: ' + outputID
        predictions = 'static/GIDBcache/predictions/' + outputID + '.svg'
        preddelta = 'static/GIDBcache/preddelta/' + outputID + '.svg'
        shapref = 'static/GIDBcache/shapref/' + outputID + '.svg'
        shapalt = 'static/GIDBcache/shapalt/' + outputID + '.svg'
        shapdelta = 'static/GIDBcache/shapdelta/' + outputID + '.svg'
        motifs = 'static/GIDBcache/motifs/' + outputID + '.txt'
        print("DONE START PHASE", strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        # if (not(exists(predictions) and exists(preddelta) and exists(shapref) and exists(shapalt) and exists(shapdelta) and exists(motifs))):
        generate_output_compound(cell_type, chrom, center, positions, alleleA, alleleB, outputID)
        motifhtml = create_table(outputID)
        motiftable = Markup(motifhtml)
        print("OUTPUT SENT", strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        return render_template('outputGIDB.html', predictions = predictions, preddelta = preddelta, shapalt = shapalt, shapref = shapref, shapdelta = shapdelta, table=motiftable)
    return render_template('indexGIDB.html', form=form)

if __name__ == '__main__':
    app.run(port='8000')
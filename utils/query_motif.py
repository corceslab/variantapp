import os
import subprocess
import numpy as np
import pandas as pd

import csv
import io

def path_to_image_html(best_model):
    return '<img src="static/images/motifs/'+ best_model + '" style="max-height:100px;"/>'

def get_motifs(peaks_df):
    #print(chrom, loc)
    # utils_dir='/home/ubuntu/variantapp/utils'
    # input_chrom = str(chrom)
    # input_loc = str(loc)
    #print(input_chrom, input_loc)
    # query = "\'" + str(chrom) + '\:' + str(loc) + '\-' + str(loc) + "\'"
    #print('Querying motifs from Vierstra')
    motifs_main = pd.DataFrame(columns = ['variant', 'chromosome', 'start', 'end', 'motif_cluster', 'match_score', 'strand', 'best_model', 'num_models'])
    for i, g in peaks_df.groupby(peaks_df.index // 2):
        input_chrom = str(g.iloc[0]['chrom'])
        input_loc = str(g.iloc[0]['st'])
        subprocess.call(['sh', 'utils/query_motif.sh', '-c', input_chrom, '-l', input_loc])
        #os.system("sh utils/query_motif.sh -c " + input_chrom + " -l " + input_loc)
        # os.system("/home/ubuntu/variantapp/utils/query_motif.sh -a " + query)
        # subprocess.call(['sh query_motif.sh ', input_chrom, input_loc], shell=True, cwd=utils_dir) #../
        #print('Finished querying motifs')
        motifs=pd.read_csv('static/text/motifs.txt', sep='\t', header=None, names=['chromosome', 'start', 'end', 'motif_cluster', 'match_score', 'strand', 'best_model', 'num_models']) #../
        # pd.set_option('colheader_justify', 'center')
        #print(motifs.info())
        motifs = motifs.sort_values(by=['match_score'], ascending=False)
        motifs['image'] = motifs['best_model'] + '.png'
        motifs.insert(0, 'new-col', g)
        motifs_main = motifs_main.append(motifs)
    #motifs.to_html(open('../templates/motifs.html', 'w'))
    html_string = '''
    <html>
    <head><title>HTML Motifs Dataframe with CSS</title></head>
    <link rel="stylesheet" type="text/css" href="/static/css/df_style.css"/>
    <body>
        {table}
    </body>
    </html>
    '''
    export = io.StringIO()
    export.write("Motifs Table:\n")
    csv.writer(export).writerows(motifs_main.values.tolist())
    export.write("\n")
    export.write("\n")

    # with open('static/csv/data.csv', 'a') as f:
    #     motifs.to_csv(f, header=True, index=True)
    # with open('templates/motifs.html', 'w') as f: #../
    #     f.write(html_string.format(table=motifs.to_html(classes='mystyle', escape=False, formatters=dict(image=path_to_image_html), index=False)))
    s = html_string.format(table=motifs_main.to_html(classes='mystyle', escape=False, formatters=dict(image=path_to_image_html), index=False))
    #print(s)
    return s, export

def get_motif(input_chrom, input_loc):
    motifs_main = pd.DataFrame(columns = ['chromosome', 'start', 'end', 'motif_cluster', 'match_score', 'strand', 'best_model', 'num_models'])
    subprocess.call(['sh', 'utils/query_motif.sh', '-c', input_chrom, '-l', input_loc])
    motifs=pd.read_csv('static/text/motifs.txt', sep='\t', header=None, names=['chromosome', 'start', 'end', 'motif_cluster', 'match_score', 'strand', 'best_model', 'num_models']) #../
    motifs = motifs.sort_values(by=['match_score'], ascending=False)
    motifs['image'] = motifs['best_model'] + '.png'
    motifs = motifs[motifs['match_score']>8]
    motifs_main = motifs_main.append(motifs)
    html_string = '''
    <html>
    <head><title>HTML Motifs Dataframe with CSS</title></head>
    <link rel="stylesheet" type="text/css" href="/static/css/df_style.css"/>
    <body>
        {table}
    </body>
    </html>
    '''
    s = html_string.format(table=motifs_main.to_html(classes='mystyle', escape=False, formatters=dict(image=path_to_image_html), index=False))
    return s

if __name__ == '__main__':
    get_motifs('chr3', 52498434)
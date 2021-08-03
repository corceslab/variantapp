import os
import subprocess
import numpy as np
import pandas as pd

def path_to_image_html(best_model):
    return '<img src="static/images/motifs/'+ best_model + '" style="max-height:100px;"/>'

def get_motifs(chrom, loc):
    print(chrom, loc)
    utils_dir='/home/ubuntu/variantapp/utils'
    input_chrom = str(chrom)
    input_loc = str(loc)
    print(input_chrom, input_loc)
    # query = "\'" + str(chrom) + '\:' + str(loc) + '\-' + str(loc) + "\'"
    print('Querying motifs from Vierstra')
    subprocess.call(['sh', 'utils/query_motif.sh', '-c', input_chrom, '-l', input_loc])
    #os.system("sh utils/query_motif.sh -c " + input_chrom + " -l " + input_loc)
    # os.system("/home/ubuntu/variantapp/utils/query_motif.sh -a " + query)
    # subprocess.call(['sh query_motif.sh ', input_chrom, input_loc], shell=True, cwd=utils_dir) #../
    print('Finished querying motifs')
    motifs=pd.read_csv('static/text/motifs.txt', sep='\t', header=None, names=['chromosome', 'start', 'end', 'motif_cluster', 'match_score', 'strand', 'best_model', 'num_models']) #../
    pd.set_option('colheader_justify', 'center')
    print(motifs.info())
    motifs = motifs.sort_values(by=['match_score'], ascending=False)
    motifs['image'] = motifs['best_model'] + '.png'
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
    with open('static/csv/data.csv', 'a') as f:
        motifs.to_csv(f, header=True, index=True)
    with open('templates/motifs.html', 'w') as f: #../
        f.write(html_string.format(table=motifs.to_html(classes='mystyle', escape=False, formatters=dict(image=path_to_image_html), index=False)))

if __name__ == '__main__':
    get_motifs('chr3', 52498434)
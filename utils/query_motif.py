import subprocess
import numpy as np
import pandas as pd

def get_motifs(chrom, loc):
    query = str(chrom) + ":" + str(loc) + "-" + str(loc)
    subprocess.call(['sh', 'query_motif.sh', str(query)])
    motifs=pd.read_csv('static/text/motifs.txt', sep='\t', header=None, names=['chromosome', 'start', 'end', 'motif', 'match_score', 'strand', 'sequence'])
    pd.set_option('colheader_justify', 'center')
    #motifs.to_html(open('../templates/motifs.html', 'w'))
    html_string = '''
    <html>
    <head><title>HTML Motifs Dataframe with CSS</title></head>
    <link rel="stylesheet" type="text/css" href="/static/css/df_style.css"/>
    <body>
        {table}
    </body>
    </html>.
    '''
    with open('templates/motifs.html', 'w') as f:
        f.write(html_string.format(table=motifs.to_html(classes='mystyle')))

if __name__ == '__main__':
    get_motifs('chr3', 52498434)
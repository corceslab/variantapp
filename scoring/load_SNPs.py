import numpy as np
import pandas as pd

def load_variants():
    """ Load variants from csv file containing GWAS AD data
        Method used for variant scoring 
    """
    vars_df = pd.read_csv('data/ADPD_SNPs.csv')
    return vars_df
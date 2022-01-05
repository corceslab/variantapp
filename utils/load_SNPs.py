import numpy as np
import pandas as pd

def load_variants():
    vars_df = pd.read_csv('data/ADPD_SNPs.csv')
    return vars_df
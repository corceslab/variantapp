import pandas as pd
import numpy as np
import tensorflow as tf
import pysam
from basepairmodels.cli.losses import MultichannelMultinomialNLL
from basepairmodels.cli.losses import multinomial_nll
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope
from mseqgen.sequtils import one_hot_encode
from scipy.special import softmax
from matplotlib import pyplot as plt

from gen_prediction import insert_variant

def get_critical(df):
    return df[['chrom', 'st', 'Effect_Allele', 'Noneffect_Allele']]

def load_data(filepath):
    df = pd.read_csv(filepath, \
        names=['SNP_rsID', 'Locus_Number', 'Disease', 'SNP_Type', 'chrom', \
            'st', 'hg19_Position', 'r2_with_LD_Tag', 'LD_Tag_chr', 'LD_Tag_pos_hg38', \
            'Effect_Allele', 'Noneffect_Allele', 'Direction_of_Association', 'GWAS_pvalue', \
            'In_Any_Peak', 'In_scATAC_Cluster_Peak', 'Has_Coloc', 'Only_Coloc', 'Has_eqtl_coloc', \
            'Has_sqtl_coloc', 'Coloc_Effect_Genes', 'Has_HiChIP_Cicero_link', \
            'HiChIP_Cicero_Linked_Genes', 'Has_ML_prediction', 'ML_confidence', 'ML_sig_clusters', \
            'ML_sig_celltypes', 'Has_Allelic_Imbalance', 'RASQUAL_region', 'RASQUAL_effect_size', \
            'RASQUAL_permutation_significant', 'RASQUAL_numSigRegions'])
    print(df.info())
    vals = get_critical(df)
    vals.dropna(inplace=True)

    vals = vals[(vals.Effect_Allele == 'A') | (vals.Effect_Allele == 'C') | (vals.Effect_Allele == 'G') | (vals.Effect_Allele == 'T')]
    vals = vals[(vals.Noneffect_Allele == 'A') | (vals.Noneffect_Allele == 'C') | (vals.Noneffect_Allele == 'G') | (vals.Noneffect_Allele == 'T')]
    
    vals = vals.astype({'chrom': 'str', 'st': 'int', 'Effect_Allele': 'str', 'Noneffect_Allele': 'str'})
    print(vals.info())
    print(vals['chrom'].unique())
    print(vals['Effect_Allele'].unique())
    print(vals['Noneffect_Allele'].unique())

    peaks_df = pd.DataFrame(columns=['chrom', 'st', 'allele'])
    for idx, row in vals.iterrows():
        peaks_df.loc[len(peaks_df.index)] = [row['chrom'], row['st'], row['Effect_Allele']]
        peaks_df.loc[len(peaks_df.index)] = [row['chrom'], row['st'], row['Noneffect_Allele']]
    peaks_df['summit'] = 0
    peaks_df['signalValue'] = 10
    peaks_df['start'] = peaks_df['st'] + peaks_df['summit'] - (2114 // 2)
    peaks_df['end'] = peaks_df['st'] + peaks_df['summit'] + (2114 // 2)
    #peaks_df = peaks_df.reset_index(drop=True)
    print(peaks_df.info())
    print(peaks_df.head())
    return peaks_df

def load_sequences(peaks_df):
    fasta_ref = pysam.FastaFile('../reference/hg38.genome.fa')
    sequences = []
    for idx, row in peaks_df.iterrows():
        start = row['start']
        end = row['end']
        peak_loc = 1057
        allele = row['allele']
        seq = fasta_ref.fetch(row['chrom'], start, end).upper()
        seq = insert_variant(seq, allele, peak_loc)
        if(len(seq)!= 2114):
            continue
        sequences.append(seq)
    return sequences

def run_batch_main(datapath, modelpath):
    with CustomObjectScope({'MultichannelMultinomialNLL': MultichannelMultinomialNLL}):
        model = load_model(modelpath)
    peaks_df = load_data(datapath)
    sequences = load_sequences(peaks_df)
    X = one_hot_encode(sequences, 2114)
    num_entries = X.shape[0]
    control_profile = np.zeros((num_entries, 1000, 2))
    control_logcount = np.zeros(num_entries)
    batch = {
        "sequence": X,
        "control_profile": control_profile,
        "control_logcount": control_logcount
    }
    preds = model.predict(batch, batch_size = num_entries)
    print(preds)


if __name__ == '__main__':
    run_batch_main('../static/csv/ADPD_SNP_List.csv', '../models/C24/model.h5')
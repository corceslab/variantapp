import pandas as pd
import tensorflow as tf

from utils.utils import generate_explain_score_rankings
from utils.load_model import load_chrombpnet
from scoring.scoringV2 import gen_MPRA_preds


def load_variants():
    vars_df = pd.read_csv('data/ADVariants_Ryan.csv')
    return vars_df

def run_cluster(cluster):
    print("Running Cluster", cluster)
    vars_df = load_variants()
    # print(vars_df.info())
    # vars_df = vars_df[vars_df['Has_ML_prediction'] == True]
    vars_df = vars_df[['SNP_rsID', 'Disease', 'hg38_Chromosome', 'hg38_Position', 'Effect_Allele', 'Noneffect_Allele', 'ML_confidence']] # , 'ML_clusters']]
    # vars_df = vars_df.sort_values('ML_confidence', ascending=False)
    vars_df = vars_df[vars_df['Disease']=='AD']
    vars_df = vars_df[vars_df['SNP_rsID'].str[:2]=='rs']
    filtered_vars_df = vars_df #[vars_df['ML_clusters'].str.find(',' + cluster + ',') != -1]
    # print(vars_df.info())
    # print(vars_df.head())
    # print(filtered_vars_df.head())
    filtered_vars_df.reset_index(drop=True, inplace=True)
    peaks_df = pd.DataFrame(columns = ['chrom', 'st', 'effect', 'noneffect'])
    peaks_df['chrom'] = filtered_vars_df['hg38_Chromosome']
    peaks_df['st'] = filtered_vars_df['hg38_Position']
    peaks_df['effect'] = filtered_vars_df['Effect_Allele']
    peaks_df['noneffect'] = filtered_vars_df['Noneffect_Allele']
    print(peaks_df)
    all_rsIDs = filtered_vars_df[['SNP_rsID']].values.tolist()
    rsIDs = ''
    for rsID in all_rsIDs:
        if rsID[0][:2] == 'rs':
            rsIDs += rsID[0] + ", "
    rsIDs = rsIDs[:-2]
    print(rsIDs)
    # lfc_scores = generate_explain_score('C' + cluster, rsIDs)
    lfc_scores = generate_explain_score_rankings('C' + cluster, rsIDs, peaks_df)
    lfc_scores.reset_index(drop=True, inplace=True)
    # print(lfc_scores)
    scores_df = filtered_vars_df
    scores_df['lfc'] = lfc_scores['lfc']
    scores_df['abs_lfc'] = lfc_scores['abs_lfc']
    scores_df['jsd'] = lfc_scores['jsd']
    scores_df['alt_scores'] = lfc_scores['alt_scores']
    scores_df['ref_scores'] = lfc_scores['ref_scores']
    scores_df['max_alleles'] = lfc_scores['max_alleles']
    scores_df['cluster'] = cluster
    # scores_df = scores_df[scores_df['max_alleles']>=60]
    scores_df = scores_df.sort_values(by=['jsd'], ascending=[False])
    # print(scores_df)
    # scores_df = scores_df.sort_values('jsd', ascending=False)
    # print(scores_df)
    # generate_explain_score('C' + cluster, rsIDs, nc)
    return scores_df

def run_MPRAseqs():
    seqA_df = pd.read_csv('data/MPRA/seqAN.csv')
    seqB_df = pd.read_csv('data/MPRA/seqBN.csv')
    seqA = seqA_df['seqA']
    seqB = seqB_df['seqB']
    metadata = seqA_df['metadata']
    clusters = ['C2', 'C5', 'C8', 'C13', 'C19', 'C24']
    results_df = gen_MPRA_preds('C1', seqA, seqB, metadata)
    for cluster in clusters:
        results_df = results_df.append(gen_MPRA_preds(cluster, seqA, seqB, metadata), ignore_index=True)
    results_df.to_csv('MPRA_all_lfc_seqN.csv')

if __name__ == '__main__':
    run_MPRAseqs()

# For GWAS Variant Prioritization
# if __name__ == '__main__':
#     clusters = ['2', '5', '8', '13', '19', '24']
#     results = run_cluster('1')
#     # print("RESULTS AT CLUSTER", '1', results)
#     for cluster in clusters:
#         results = results.append(run_cluster(cluster), ignore_index=True)
#         # print("RESULTS AT CLUSTER", cluster, results)
#     results = results.sort_values(by=['jsd'], ascending=[False])
#     # results = results.drop_duplicates(subset=['SNP_rsID'])
#     results = results.reset_index(drop=True)
#     print('\n\nFINAL RESULTS\n\n', results)
#     results.to_csv('data/output/scored_ADVariants_ISEF.csv')
    
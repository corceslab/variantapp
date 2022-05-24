import pandas as pd

if __name__ == '__main__':
    results = pd.read_csv('../data/output/scored_ADVariants_Ryan.csv')
    results = results.drop_duplicates(subset=['SNP_rsID'])
    results = results.reset_index(drop=True)
    print('\n\nFINAL RESULTS\n\n', results)
    results.to_csv('../data/output/scored_ADVariants_Ryan_deduped.csv')
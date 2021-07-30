import myvariant
import pandas as pd


def add_variant(df, chromosome, position, allele):
    df2 = pd.DataFrame([[chromosome, position, allele]], \
        columns = ['chrom', 'st', 'allele'])
    return df.append(df2)

def query_rsID(rsID):
    mv = myvariant.MyVariantInfo()
    rsID_info = mv.query('dbsnp.rsid:'+rsID, fields='vcf', assembly='hg38')
    hits = rsID_info['hits']
    data = hits[0]
    #print("data _id:", data['_id'])
    chrom = data['_id'].split(':')[0]
    #print("chrom:", chrom)
    vcf = data['vcf']
    position = int(vcf['position'])
    effect_allele = vcf['alt']
    noneffect_allele = vcf['ref']
    var_df = pd.DataFrame(columns = ['chrom', 'st', 'allele'])
    var_df = add_variant(var_df, chrom, position, effect_allele)
    var_df = add_variant(var_df, chrom, position, noneffect_allele)
    var_df['summit'] = 0
    var_df['signalValue'] = 10
    var_df['start'] = var_df['st'] + var_df['summit'] - (2114 // 2)
    var_df['end'] = var_df['st'] + var_df['summit'] + (2114 // 2)
    print(var_df.head())
    return var_df


def query_values(chrom, position, effect_allele, noneffect_allele):
    var_df = pd.DataFrame(columns = ['chrom', 'st', 'allele'])
    var_df = add_variant(var_df, chrom, position, effect_allele)
    var_df = add_variant(var_df, chrom, position, noneffect_allele)
    var_df['summit'] = 0
    var_df['signalValue'] = 10
    var_df['start'] = var_df['st'] + var_df['summit'] - (2114 // 2)
    var_df['end'] = var_df['st'] + var_df['summit'] + (2114 // 2)
    return var_df
    


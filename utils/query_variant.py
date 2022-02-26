import myvariant
import pandas as pd


def add_variant(df, chromosome, position, allele):
    df2 = pd.DataFrame([[chromosome, position, allele]], \
        columns = ['chrom', 'st', 'allele'])
    return df.append(df2)

def query_rsID(rsID):
    mv = myvariant.MyVariantInfo()
    rsIDs = rsID.split(", ")
    print(len(rsIDs))
    var_df = pd.DataFrame(columns = ['chrom', 'st', 'allele'])
    bad_SNPs = ['rs199504614','rs142767245', 'rs113242154']
    counter = 0
    for id in rsIDs:
        counter+=1
        if(id in bad_SNPs):
            continue
        # print("querying: ", id)
        rsID_info = mv.query('dbsnp.rsid:'+id, fields='vcf', assembly='hg38')
        hits = rsID_info['hits']
        data = hits[0]
        chrom = data['_id'].split(':')[0]
        vcf = data['vcf']
        position = int(vcf['position'])
        alt_allele = vcf['alt']
        ref_allele = vcf['ref']
        # print("alt:", alt_allele)
        # print("ref:", ref_allele)
        # alt alleles are in even rows (2k), ref alleles are in odd rows (1+2k)
        var_df = add_variant(var_df, chrom, position, alt_allele)
        var_df = add_variant(var_df, chrom, position, ref_allele)
        if(counter % 100 == 0): print(counter)

    var_df['summit'] = 0
    var_df['signalValue'] = 10
    var_df['start'] = var_df['st'] + var_df['summit'] - (2114 // 2)
    var_df['end'] = var_df['st'] + var_df['summit'] + (2114 // 2)
    var_df = var_df.reset_index(drop=True)
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
    var_df = var_df.reset_index(drop=True)
    return var_df
    


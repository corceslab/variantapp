import myvariant
import pandas as pd


def add_variant(df, chromosome, position, allele):
    """ Inserts a variant (with chromosome, position, and allele information) into a Pandas dataframe
    """
    df2 = pd.DataFrame([[chromosome, position, allele]], \
        columns = ['chrom', 'st', 'allele'])
    return df.append(df2)

def query_rsID(rsID):
    """ For each rsID, query the relevant information from the MyVariant package
        Add each allele of the variant into the 'var_df' pandas dataframe, with
        alt alleles in even rows (added first) and ref alleles in odd rows (added second)
        Returns a pandas dataframe with an allele of a variant in each row
    """
    mv = myvariant.MyVariantInfo()
    rsIDs = rsID.split(", ")
    print(len(rsIDs))
    var_df = pd.DataFrame(columns = ['chrom', 'st', 'allele'])
    bad_SNPs = ['rs199504614','rs142767245', 'rs113242154', 'rs144707726']
    counter = 0
    for id in rsIDs:
        counter+=1
        if(id in bad_SNPs):
            continue
        print("querying: ", id)
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
    """ For each variant, adds two rows to the var_df dataframe, one for
        each allele, following the alt (0) / ref (1) standard
    """
    var_df = pd.DataFrame(columns = ['chrom', 'st', 'allele'])
    var_df = add_variant(var_df, chrom, position, effect_allele)
    var_df = add_variant(var_df, chrom, position, noneffect_allele)
    var_df['summit'] = 0
    var_df['signalValue'] = 10
    var_df['start'] = var_df['st'] + var_df['summit'] - (2114 // 2)
    var_df['end'] = var_df['st'] + var_df['summit'] + (2114 // 2)
    var_df = var_df.reset_index(drop=True)
    return var_df

def query_values_scoring(vars_df):
    """ Creating the peaks_df file for bulk variant scoring
        Compile key information (chrom, st, allele) from the original dataframe of variants
    """
    peaks_df = pd.DataFrame(columns = ['chrom', 'st', 'allele'])
    for i, g in vars_df.groupby(vars_df.index):
        # print("g", g)
        # print("g['chrom']", g['chrom'].iloc[0])
        # print("g['st']", g['st'])
        # print("g['effect']", g['effect'])
        peaks_df = add_variant(peaks_df, g['chrom'].iloc[0], g['st'].iloc[0], g['effect'].iloc[0])
        peaks_df = add_variant(peaks_df, g['chrom'].iloc[0], g['st'].iloc[0], g['noneffect'].iloc[0])
    peaks_df['summit'] = 0
    peaks_df['signalValue'] = 10
    peaks_df['start'] = peaks_df['st'] + peaks_df['summit'] - (2114 // 2)
    peaks_df['end'] = peaks_df['st'] + peaks_df['summit'] + (2114 // 2)
    peaks_df = peaks_df.reset_index(drop=True)
    # print(peaks_df)
    return peaks_df
    


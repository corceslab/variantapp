import pysam
from mseqgen.sequtils import one_hot_encode

def insert_variant(seq, allele, position):
    """ Inserts the specified allele at the given position of the sequence
    """
    left, right = seq[:position-1], seq[position:]
    return left + allele + right

def load_sequences(peaks_df):
    """ Loads the sequences associated with each row of the peaks_df Pandas dataframe
        using chromosome, start position, end position, and allele information
        Returns a list of sequences queried from reference/hg38.genome.fa file
    """
    fasta_ref = pysam.FastaFile('reference/hg38.genome.fa')
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
    X = one_hot_encode(sequences, 2114)
    return X, sequences

def load_compound(chrom, center, positions, alleleA, alleleB):
    fasta_ref = pysam.FastaFile('reference/hg38.genome.fa')
    sequences = []

    #Allele A
    seq = fasta_ref.fetch(chrom, center - 1057, center + 1057).upper()
    seqA = seq
    seqB = seq
    SNP_loc = 1057
    if(len(seq)!= 2114):
        print("QUERIED SEQUENCE LENGTH INCORRECT")
    for i in range(len(positions)):
        seqA = insert_variant(seqA, alleleA[i], SNP_loc + positions[i])
        seqB = insert_variant(seqB, alleleB[i], SNP_loc + positions[i])
    sequences = [seqA, seqB]
    X = one_hot_encode(sequences, 2114)
    return X, sequences

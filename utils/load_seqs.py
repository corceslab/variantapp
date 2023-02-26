import pysam
from utils.sequtils import one_hot_encode

def insert_variant(seq, allele, position):
    """ Inserts the specified allele at the given position of the sequence
    """
    left, right = seq[:position-1], seq[position:]
    return left + allele + right

def load_compound(chrom, center, positions, alleleA, alleleB):
    fasta_ref = pysam.FastaFile('reference/hg38.genome.fa')
    sequences = []

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
